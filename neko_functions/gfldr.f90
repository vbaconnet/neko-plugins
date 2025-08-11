subroutine gfldr(file_name,u,v,w,p,interpolate,tolerance,mesh_file_name, sync_memcpy)
  character(len=*), intent(in) :: file_name
  type(field_t), intent(inout) :: u,v,w,p
  logical, intent(in), optional :: interpolate
  real(kind=rp), intent(in), optional :: tolerance
  character(len=*), intent(inout), optional :: mesh_file_name
  logical, intent(in), optional :: sync_memcpy

  character(len=LOG_SIZE) :: log_buf
  integer :: sample_idx, sample_mesh_idx
  integer :: last_index
  type(fld_file_data_t) :: fld_data
  type(file_t) :: f
  real(kind=rp) :: tolerance_
  logical :: mesh_mismatch, interpolate_, sync_memcpy_
  character(len=1024) :: mesh_fname_

  ! ---- For the mesh to mesh interpolation
  type(global_interpolation_t) :: global_interp
  ! -----
 
  ! ---- For space to space interpolation
  type(space_t) :: prev_Xh
  type(interpolator_t) :: space_interp
  ! ----

  ! --- Setup of default options

  ! Default is to interpolate
  if (present(interpolate)) then
    interpolate_ = interpolate
  else
    interpolate_ = .true.
  end if
  
  ! Default is to not sync
  if (present(sync_memcpy)) then
    sync_memcpy_ = sync_memcpy
  else
    sync_memcpy_ = .false.
  end if

  ! Default is tolerance 1e-6
  if (present(tolerance)) then
    tolerance_ = tolerance
  else
    tolerance_ = 0.000001_rp
  end if
  
  ! Default is mesh_file_name = "none"
  if (present(mesh_file_name)) then
    mesh_fname_ = mesh_file_name
  else
    mesh_fname_ = "none"
  end if
 
  ! -------------------

  ! ---------- EXECUTION ------- 
  call neko_log%section("(GFLDR)", LVL=NEKO_LOG_INFO)
  call neko_log%message("(GFLDR) File name     : " // trim(file_name))
  write (log_buf, '(A,L1)') "(GFLDR) Interpolation : ", interpolate
  call neko_log%message(log_buf, LVL=NEKO_LOG_INFO)
 
  ! Extract sample index from the file name
  sample_idx = extract_fld_file_index(file_name, -1)
 
  if (sample_idx .eq. -1) &
       call neko_error("(GFLDR) Invalid file name. The file format must be e.g. 'mean0.f00001'")
 
  ! Change from "field0.f000*" to "field0.fld" for the fld reader
  call filename_chsuffix(file_name, file_name, 'fld')
 
  call fld_data%init
  f = file_t(trim(file_name))
 
  if (interpolate) then
 
     ! If no mesh file is specified, use the default file name
     if (mesh_file_name .eq. "none") then
        mesh_file_name = trim(file_name)
        sample_mesh_idx = sample_idx
     else
 
        ! Extract sample index from the mesh file name
        sample_mesh_idx = extract_fld_file_index(mesh_file_name, -1)
 
        if (sample_mesh_idx .eq. -1) then
           call neko_error("(GFLDR) Invalid file name. &
             &The file format must be e.g. 'mean0.f00001'")
        end if
 
        write (log_buf, '(A,ES12.6)') "(GFLDR) Tolerance     : ", tolerance
        call neko_log%message(log_buf, LVL=NEKO_LOG_INFO)
        write (log_buf, '(A,A)') "(GFLDR) Mesh file     : ", &
             trim(mesh_file_name)
        call neko_log%message(log_buf, LVL=NEKO_LOG_INFO)
 
     end if ! if mesh_file_name .eq. none
 
     ! Avoid reading the same fld file twice (i.e. if mesh file and fld file
     ! are the same)
     if (sample_mesh_idx .ne. sample_idx) then
        call f%set_counter(sample_mesh_idx)
        call f%read(fld_data)
     end if
 
  end if interp_check
 
  ! Read the field file containing (u,v,w,p)
  call f%set_counter(sample_idx)
  call f%read(fld_data)
 
  !
  ! Check that the data in the fld file matches the current case.
  ! Note that this is a safeguard and there are corner cases where
  ! two different meshes have the same dimension and same # of elements
  ! but this should be enough to cover obvious cases.
  !
  mesh_mismatch = (fld_data%glb_nelv .ne. u%msh%glb_nelv .or. &
       fld_data%gdim .ne. u%msh%gdim)
 
  if (mesh_mismatch .and. .not. interpolate) then
     call neko_error("(GFLDR) The fld file must match the current mesh! &
       &Use 'interpolate': 'true' to enable interpolation.")
  else if (.not. mesh_mismatch .and. interpolate) then
     if (pe_rank .eq. 0) call neko_log%warning("(GFLDR) You have activated interpolation but you might &
       &still be using the same mesh.")
  end if
 
  ! Mesh interpolation if specified
  if (interpolate) then
 
     ! Issue a warning if the mesh is in single precision
     select type (ft => f%file_type)
     type is (fld_file_t)
        if (.not. ft%dp_precision) then
           if (pe_rank .eq. 0) call neko_warning("(GFLDR) The coordinates read from the field file are &
             &in single precision.")
           call neko_log%message("(GFLDR) It is recommended to use a mesh in double &
             &precision for better interpolation results.")
           call neko_log%message("(GFLDR) If the interpolation does not work, you&
             &can try to increase the tolerance.")
        end if
     class default
     end select
 
     ! Generates an interpolator object and performs the point search
     call fld_data%generate_interpolator(global_interp, u%dof, u%msh, &
          tolerance)
 
     ! Evaluate velocities and pressure
     call global_interp%evaluate(u%x, fld_data%u%x)
     call global_interp%evaluate(v%x, fld_data%v%x)
     call global_interp%evaluate(w%x, fld_data%w%x)
     call global_interp%evaluate(p%x, fld_data%p%x)
 
     call global_interp%free
 
  else ! No interpolation, but potentially just from different spaces
 
     ! Build a space_t object from the data in the fld file
     call prev_xh%init(gll, fld_data%lx, fld_data%ly, fld_data%lz)
     call space_interp%init(u%Xh, prev_xh)
 
     ! Do the space-to-space interpolation
     call space_interp%map_host(u%x, fld_data%u%x, fld_data%nelv, u%Xh)
     call space_interp%map_host(v%x, fld_data%v%x, fld_data%nelv, u%Xh)
     call space_interp%map_host(w%x, fld_data%w%x, fld_data%nelv, u%Xh)
     call space_interp%map_host(p%x, fld_data%p%x, fld_data%nelv, u%Xh)
 
     call space_interp%free
 
  end if
 
  ! NOTE: we do not copy (u,v,w) to the GPU since `set_flow_ic_common` does
  ! the copy for us
  if (neko_bcknd_device .eq. 1) then
     call device_memcpy(u%x, u%x_d, u%dof%size(), &
             host_to_device, sync = sync_memcpy_)
     call device_memcpy(v%x, v%x_d, p%dof%size(), &
             host_to_device, sync = sync_memcpy_)
     call device_memcpy(w%x, w%x_d, p%dof%size(), &
             host_to_device, sync = sync_memcpy_)
     call device_memcpy(p%x, p%x_d, p%dof%size(), &
             host_to_device, sync = sync_memcpy_)
  end if
  call fld_data%free
 
end subroutine gfldr
