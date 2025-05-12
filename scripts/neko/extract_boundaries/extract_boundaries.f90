program extract_boundaries
  use neko
  use facet_zone
  use mpi_f08
  implicit none

  character(len=NEKO_FNAME_LEN) :: mesh_fname, field_fname, case_fname, out_fname
  integer :: boundary_label
  character(len=LOG_SIZE) :: log_buf
  character(len=5) :: suffix
  type(file_t) :: field_file, mesh_file, out_file
  type(fld_file_data_t) :: field_data
  type(mesh_t) :: mesh
  type(json_file) :: params
  character(len=:), allocatable :: string_val
  integer :: argc
  type(matrix_t) :: out, global_out, global_out_trsp
  integer :: ierr, i, lx, labeled_zone_index = -1, npts_global
  integer, allocatable :: npts_local(:), offsets(:)
  character(len=NEKO_FNAME_LEN) :: inputchar
  type(facet_zone_t) :: mesh_zone
character(len=NEKO_MSH_MAX_ZLBL_LEN), allocatable :: bc_labels(:)
  type(coef_t) :: coef
  type(space_t) :: Xh
  type(dofmap_t) :: dof
  type(gs_t) :: gs_h

  integer :: NROWS = 10

  type(dirichlet_t) :: target_bc

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
  field_file = file_t(trim(field_fname), precision = dp)
  call field_data%init
  ! --------------------

  ! --- Read case file
  call get_command_argument(2, case_fname)
  call filename_suffix(case_fname, suffix)

  if (suffix .eq. "case") then
     call params%load_file(filename = trim(case_fname))
     call json_get(params, 'case.mesh_file', string_val)

     mesh_fname = trim(string_val)
  else if (suffix .eq. "nmsh") then
     mesh_fname = trim(case_fname)
  else
     call neko_error(suffix // " : invalid file type, only .case allowed")
     call usage()
     stop
  end if

  mesh_file = file_t(trim(mesh_fname))
  call mesh_file%read(mesh)
  ! ---------------------------

  ! Read name of the boundary to extract (inlet, outlet...)
  call get_command_argument(3, inputchar)
  read(inputchar, *) boundary_label
  write (out_fname, '(A,I2.2,A)') "extracted", boundary_label,".csv"
  out_file = file_t(trim(out_fname))

  !
  ! Start processing the field data
  !

  ! READ first fiel file for processing field file data
  call field_file%read(field_data)
  
  !
  ! Generate coef xh etc -------
  !
  call Xh%init(GLL, field_data%lx, field_data%ly, field_data%lz)
  call dof%init(mesh, Xh)
  call gs_h%init(dof)
  call coef%init(gs_h) 

  call target_bc%init_base(coef)
  call target_bc%mark_zone(mesh%labeled_zones(boundary_label))
  call target_bc%finalize()
  ! ----------------------------

  print *, pe_rank, "has", target_bc%msk(0), "points"

  allocate(npts_local(0:(pe_size-1)))
  allocate(offsets(0:(pe_size-1)))
  offsets = 0
  npts_local = 0

  ! Get number of points per rank for everyone
  npts_local(pe_rank) = target_bc%msk(0)
  call MPI_Allreduce(MPI_IN_PLACE, npts_local, pe_size, MPI_INTEGER, MPI_SUM, &
          NEKO_COMM, ierr)

  ! Calculate the local offset based on number of points from other ranks
  do i = 1, pe_size-1
     offsets(i) = offsets(i-1) + npts_local(i-1)
  end do

  npts_global = sum(npts_local)

  if (npts_global .eq. 0) call neko_error("npts global is 0")
  call out%init(NROWS, npts_local(pe_rank))

  print *, npts_global

  if (pe_rank .eq. 0) then
        call global_out%init(NROWS, npts_global) !x,y,z,u,v,w,p,t,s1,s2
        call global_out_trsp%init(npts_global, NROWS) !x,y,z,u,v,w,p
        print *, global_out%ncols
        call out_file%set_header("t,x,y,z,u,v,w,p,t,s1,s2")
  end if

  ! extract and write inlet values
  call extract_values(field_data, target_bc%msk, target_bc%msk(0), out)
  
  call MPI_Gatherv(out%x, NROWS*npts_local(pe_rank), MPI_DOUBLE_PRECISION, &
          global_out%x, NROWS*npts_local, NROWS*offsets, MPI_DOUBLE_PRECISION, &
          0, NEKO_COMM, ierr)

  if (pe_rank .eq. 0) then
     call trsp(global_out_trsp%x, npts_global, global_out%x, NROWS)
     call out_file%write(global_out_trsp, field_data%time)
  end if

  do i = 1, field_data%meta_nsamples - 1
  
     call field_file%read(field_data)
  
     ! extract and write inlet values
     call extract_values(field_data, target_bc%msk, target_bc%msk(0), out)

     call MPI_Gatherv(out%x, NROWS*npts_local(pe_rank), MPI_DOUBLE_PRECISION, &
          global_out%x, NROWS*npts_local, NROWS*offsets, MPI_DOUBLE_PRECISION, &
          0, NEKO_COMM, ierr)

     if (pe_rank .eq. 0) then
        call trsp(global_out_trsp%x, npts_global, global_out%x, NROWS)
        call out_file%write(global_out_trsp, field_data%time)
     end if
     
  end do

  !
  ! Free stuff
  !
  deallocate(offsets)
  deallocate(npts_local)
  call out%free

  print *, global_out%ncols
  print *, global_out%ncols
  !call global_out%free
  !call global_out_trsp%free

  call file_free(mesh_file)
  call file_free(out_file)
  call file_free(field_file)
  call neko_finalize

contains

  subroutine usage()
     write(*,*) '--------------------------------------------------------------'
     write(*,*) 'Usage:' 
     write(*,*) ' ./extract_boundaries <field.fld> <case file> <boundary>'
     write(*,*) ''
     write(*,*) 'Examples:'
     write(*,*) ' ./extract_boundaries stats0.fld mycase.case "w"'
     write(*,*) ''
     write(*,*) 'Valid boundaries:'
     write(*,*) ' Basically any boundary labels supported in Neko'
     write(*,*) ' "w", "o", "on", "sym", "d_vel_*", "d_pres", etc...'
     write(*,*) '-------------------------------------------------------------'
  end subroutine usage

  !> Extract u,v,w,p at all GLL points that are in a specified zone.
  !! @param field_data Field data to read from
  !! @param mesh_zone Zone that defines which points to extract (e.g.
  !! mesh%inlet, mesh%outlet...)
  !! @param out_values matrix of x,y,z,u,v,w,p for each point
  !! @note There is no need to initialize `out_values`
  subroutine extract_values(fdata, bc_mask, n, out_values)
    type(fld_file_data_t), intent(in) :: fdata
    integer, intent(in) :: bc_mask(0:n)
    integer, intent(in) :: n
    type(matrix_t), intent(inout) :: out_values
    character(len=LOG_SIZE) :: log_buf

    integer :: i, idx

    write(log_buf,'(A,F10.6)') "Extracting values at t = ", fdata%time
    call neko_log%message(log_buf)

    if (n .eq. 0) return

    do idx = 1, n
       i = bc_mask(idx)

       out_values%x(1, idx) = fdata%x%x(i)
       out_values%x(2, idx) = fdata%y%x(i)
       out_values%x(3, idx) = fdata%z%x(i)
       out_values%x(4, idx) = fdata%u%x(i)
       out_values%x(5, idx) = fdata%v%x(i)
       out_values%x(6, idx) = fdata%w%x(i)
       out_values%x(7, idx) = fdata%p%x(i)
       !out_values%x(8, idx) = fdata%t%x(i)
       !out_values%x(9, idx) = fdata%s(1)%x(i)
       !out_values%x(10, idx) = fdata%s(2)%x(i)

    end do

  end subroutine extract_values

end program extract_boundaries
