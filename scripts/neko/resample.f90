program resample
  use neko
  use facet_zone
  use mpi_f08
  use interpolation, only: interpolator_t
  implicit none

  character(len=NEKO_FNAME_LEN) ::  field_fname, out_fname, mesh_fname
  integer :: argc, new_poly_order

  character(len=5) :: suffix
  type(mesh_t) :: mesh

  character(len=LOG_SIZE) :: log_buf
  character(len=NEKO_FNAME_LEN) :: inputchar

  type(file_t) :: field_file, out_file, mesh_file

  type(fld_file_data_t) :: field_data, new_field_data
  type(space_t) :: new_fld_xh, old_fld_xh
  type(interpolator_t) :: space_interp
  type(dofmap_t) :: new_dofmap

  integer :: lx, nelv, n, i

  call neko_init

  !
  ! Read Args
  !
  argc = command_argument_count()

  if ((argc .le. 1) .or. (argc .gt. 3)) then
     call usage()
     stop
  end if

  ! --- Read field name
  call get_command_argument(1, field_fname)
  field_file = file_t(trim(field_fname))
  call field_data%init

  ! Also initialize the new file / fld series
  out_fname = "resampled_"//trim(field_fname)
  out_file = file_t(trim(out_fname))
  call new_field_data%init

  ! -------------------- Mesh
!!$  call get_command_argument(2, mesh_fname)
!!$  call filename_suffix(mesh_fname, suffix)
!!$
!!$  if (suffix .eq. "nmsh") then
!!$     mesh_fname = trim(mesh_fname)
!!$  else
!!$     call neko_error(suffix // " : invalid file type, only .mesh allowed")
!!$     call usage()
!!$     stop
!!$  end if
!!$
!!$  mesh_file = file_t(trim(mesh_fname))
!!$  call mesh_file%read(mesh)
  ! ----------------------


  ! --------- New polynomial order to sample with
  call get_command_argument(2, inputchar)
  read(inputchar, *) new_poly_order
  if (pe_rank .eq. 0) write (*,*) "New poly order: ", new_poly_order
  new_poly_order = new_poly_order + 1
  ! ---------------------

  !
  ! Start processing the field data
  !

  ! ---- Read first fiel file for processing field file data
  call field_file%read(field_data)
  lx = field_data%lx
  call old_fld_xh%init(GLL, lx, lx, lx)
  ! --------

  ! ---- Generate the new field file stuff
  n = new_poly_order*new_poly_order*new_poly_order * field_data%nelv
  if (pe_rank .eq. 0) print *, "init new file data"

  ! Initializes (u,v,w) and copies important data
  call copy_fld_params(new_field_data, field_data, n, new_poly_order)

  ! initialize u,v,w
  call new_field_data%x%init(n)
  call new_field_data%y%init(n)
  call new_field_data%z%init(n)

  ! Need to also build a new dofmap
  call new_fld_xh%init(gll, new_poly_order, new_poly_order, new_poly_order)
  !call new_dofmap%init(mesh, new_fld_xh)

  !do i = 1,n
  !   new_field_data%x%x(i) = new_dofmap%x(i,1,1,1)
  !   new_field_data%y%x(i) = new_dofmap%y(i,1,1,1)
  !   new_field_data%z%x(i) = new_dofmap%z(i,1,1,1)
  !end do
  ! ----------

  ! Generate an interpolator between the two spaces
  call space_interp%init(new_fld_xh, old_fld_xh)

  ! Do the first interpolation
  call interpolate(field_data, new_field_data, new_fld_xh, space_interp, mesh=.true.)

  ! Write the first file
  call out_file%write(new_field_data, field_data%time)

  do i = 1, field_data%meta_nsamples - 1

     if (pe_rank .eq. 0) print *, "reading file", i
     call field_file%read(field_data)
     call interpolate(field_data, new_field_data, new_fld_xh, space_interp, &
          mesh=.false.)
     call out_file%write(new_field_data, field_data%time)

  end do

  !
  ! Free stuff
  !
  call space_interp%free

  call new_field_data%free
  call field_data%free

  call old_fld_xh%free
  call new_fld_xh%free

  call file_free(field_file)
  call file_free(out_file)
  call neko_finalize

contains

  subroutine usage()
     write(*,*) '--------------------------------------------------------------'
     write(*,*) 'Usage:' 
     write(*,*) ' ./resample <field.fld> <polynomial order>'
     write(*,*) ''
     write(*,*) 'Examples:'
     write(*,*) ' ./resample stats0.fld 3'
     write(*,*) '-------------------------------------------------------------'
  end subroutine usage

  subroutine interpolate(old_fld_data, new_fld_data, new_xh, space_interp, mesh)
    type(fld_file_data_t), intent(inout) :: old_fld_data, new_fld_data
    type(space_t), intent(inout) :: new_xh
    type(interpolator_t), intent(inout) :: space_interp
    logical, intent(in) :: mesh

    if (mesh) then
       call space_interp%map_host(new_fld_data%x%x, old_fld_data%x%x, &
            old_fld_data%nelv, new_xh)
       call space_interp%map_host(new_fld_data%y%x, old_fld_data%y%x, &
            old_fld_data%nelv, new_xh)
       call space_interp%map_host(new_fld_data%z%x, old_fld_data%z%x, &
            old_fld_data%nelv, new_xh)
    end if

    call space_interp%map_host(new_fld_data%u%x, old_fld_data%u%x, &
         old_fld_data%nelv, new_xh)
    call space_interp%map_host(new_fld_data%v%x, old_fld_data%v%x, &
         old_fld_data%nelv, new_xh)
    call space_interp%map_host(new_fld_data%w%x, old_fld_data%w%x, &
         old_fld_data%nelv, new_xh)

  end subroutine interpolate

  subroutine copy_fld_params(new_fld, old_fld, n, new_poly_order)
    type(fld_file_data_t), intent(inout) :: old_fld, new_fld
    integer, intent(in) :: n, new_poly_order

    ! This initializes only 3 fields, we still need x,y,z
    call new_field_data%init_n_fields(3, n)

    new_fld%fld_series_fname = trim(old_fld%fld_series_fname)
    new_fld%meta_start_counter = old_fld%meta_start_counter
    new_fld%meta_nsamples = old_fld%meta_nsamples

    new_fld%nelv = old_fld%nelv
    new_fld%lx = new_poly_order
    new_fld%ly = new_poly_order
    new_fld%lz = new_poly_order
    new_fld%gdim = old_fld%gdim
    new_fld%glb_nelv = old_fld%glb_nelv
    new_fld%offset_el = old_fld%offset_el

    allocate(new_fld%idx(new_fld%nelv))
    new_fld%idx = old_fld%idx

  end subroutine copy_fld_params

end program resample
