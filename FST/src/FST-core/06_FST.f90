! * Free stream turbulence implementation for neko
! * Author: Victor Baconnet (baconnet@kth.se)
! * Based in original implementation by Elektra Kluesberg, Prabal Negi
! 	and Philipp Schlatter. 
! * For FST theory see Schlatter (2001)

module FST

  !use global_params
  use fst_operator, only: fst_bc_compute
  use turbu, only: make_turbu, glb_uinf
  use utils, only: neko_error
  use field, only: field_t
  use coefs, only: coef_t
  use logger, only: neko_log, LOG_SIZE
  use point_zone, only: point_zone_t
  use math, only: masked_gather_copy
  use device_math, only: device_masked_gather_copy_0
  use num_types, only: rp, xp
  use comm, only: pe_rank, MPI_EXTRA_PRECISION, NEKO_COMM
  use neko_config, only: NEKO_BCKND_DEVICE
  use mpi_f08
  use device, only: device_map, device_memcpy, HOST_TO_DEVICE, device_get_ptr
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none


  type, public :: FST_t

     real(kind=xp) :: Uinf

     ! periodic directions
     logical :: periodic_x
     logical :: periodic_y
     logical :: periodic_z

     ! x fringe
     real(kind=rp) :: xmin
     real(kind=rp) :: xmax
     real(kind=rp) :: xstart
     real(kind=rp) :: xend
     real(kind=rp) :: x_delta_rise
     real(kind=rp) :: x_delta_fall

     ! y fringe
     real(kind=rp) :: ymin
     real(kind=rp) :: ymax
     real(kind=rp) :: ystart
     real(kind=rp) :: yend
     real(kind=rp) :: y_delta_rise
     real(kind=rp) :: y_delta_fall

     !> Total fringe amplitude
     real(kind=xp) :: fringe_max

     !> Final ramp time
     real(kind=xp) :: t_end
     real(kind=xp) :: t_start

     logical :: is_forcing
     logical :: is_bc

     integer :: nshells
     integer :: n_modes_total ! = k_length

     !> Fringe in space
     real(kind=xp), allocatable :: fringe_space(:)
     type(c_ptr) :: fringe_space_d = C_NULL_PTR

     !> Baseflows, if applying on a non-uniform inflow
     real(kind=rp), allocatable :: u_baseflow(:)
     real(kind=rp), allocatable :: v_baseflow(:)
     real(kind=rp), allocatable :: w_baseflow(:)
     type(c_ptr) :: u_baseflow_d = C_NULL_PTR
     type(c_ptr) :: v_baseflow_d = C_NULL_PTR
     type(c_ptr) :: w_baseflow_d = C_NULL_PTR

     !> Variable that is precomputed to save some time
     real(kind=xp), allocatable :: phi_0(:,:)
     type(c_ptr) :: phi_0_d = C_NULL_PTR

     real(kind=xp), allocatable :: random_vectors(:,:) ! u_hat_pn but reshaped
     type(c_ptr) :: random_vectors_d = C_NULL_PTR
     integer, allocatable :: shell(:)
     type(c_ptr) :: shell_d = C_NULL_PTR
     real(kind=xp), allocatable :: shell_amp(:)
     type(c_ptr) :: shell_amp_d = C_NULL_PTR

     real(kind=xp), allocatable :: k_y(:)
     real(kind=xp), allocatable :: k_z(:)
     real(kind=xp), allocatable :: k_x(:)
     type(c_ptr) :: k_x_d = C_NULL_PTR
     real(kind=xp), allocatable :: phase_shifts(:)

   contains

     ! ======== Init/Free procedures
     procedure, pass(this) :: init_common => FST_init_common
     procedure, pass(this) :: init_bc => FST_init_bc
     procedure, pass(this) :: read_from_files => FST_read_from_files
     procedure, pass(this) :: init_forcing => FST_init_forcing
     procedure, pass(this) :: free => FST_free_params
     ! =========================================================================
     procedure, pass(this) :: apply_baseflow => FST_apply_baseflow
     ! =========================================================================
     ! ======== Generate FST
     procedure, pass(this) :: generate_common => FST_generate_common
     procedure, pass(this) :: generate_forcing => FST_generate_forcing
     procedure, pass(this) :: generate_bc => FST_generate_bc
     ! =========================================================================
     ! ======= Apply FST forcing/BC
     procedure, pass(this) :: apply_forcing => FST_forcing_zone
     procedure, pass(this) :: apply_BC => FST_apply_BC
     ! ========================================================================
     procedure, pass(this) :: print => FST_print_params
  end type FST_t

contains

  !> Initialize all parameters
  subroutine FST_init_common(this, &
       xmin, xmax, xstart, xend, x_delta_rise, x_delta_fall, &
       ymin, ymax, ystart, yend, y_delta_rise, y_delta_fall, &
       fringe_max, &
       t_start, t_end, &
       periodic_x, periodic_y, periodic_z)
    class(FST_t), intent(inout) :: this
    real(kind=rp), intent(in) :: xstart, xend, xmin, xmax
    real(kind=rp), intent(in) :: ystart, yend, ymin, ymax
    real(kind=rp), intent(in) :: x_delta_rise
    real(kind=rp), intent(in) :: x_delta_fall
    real(kind=rp), intent(in) :: y_delta_rise
    real(kind=rp), intent(in) :: y_delta_fall
    real(kind=xp), intent(in) :: fringe_max
    real(kind=xp), intent(in) :: t_start
    real(kind=xp), intent(in) :: t_end
    logical, intent(in) :: periodic_x, periodic_y, periodic_z

    this%periodic_x = periodic_x
    this%periodic_y = periodic_y
    this%periodic_z = periodic_z

    this%xstart = xstart
    this%xend   = xend
    this%ystart = ystart
    this%yend   = yend

    this%xmin   = xmin
    this%xmax   = xmax
    this%ymin   = ymin
    this%ymax   = ymax

    this%fringe_max = fringe_max
    this%x_delta_rise = x_delta_rise!0.002
    this%x_delta_fall = x_delta_fall!0.002
    this%y_delta_rise = y_delta_rise!0.002
    this%y_delta_fall = y_delta_fall!0.002
    this%t_end = t_end
    this%t_start = t_start

    call this%print() ! show parameters
  
  end subroutine FST_init_common


  !> Initialize the FST to use with forcing.
  subroutine FST_init_forcing(this, &
       xstart, xend, x_delta_rise, x_delta_fall, &
       ystart, yend, y_delta_rise, y_delta_fall, &
       fringe_max, &
       t_start, t_end, &
       periodic_x, periodic_y, periodic_z)

    class(FST_t), intent(inout) :: this
    real(kind=rp), intent(in) :: xstart
    real(kind=rp), intent(in) :: xend
    real(kind=rp), intent(in) :: ystart
    real(kind=rp), intent(in) :: yend
    real(kind=rp), intent(in) :: x_delta_rise
    real(kind=rp), intent(in) :: x_delta_fall
    real(kind=rp), intent(in) :: y_delta_rise
    real(kind=rp), intent(in) :: y_delta_fall
    real(kind=xp), intent(in) :: fringe_max
    real(kind=xp), intent(in) :: t_start
    real(kind=xp), intent(in) :: t_end
    logical, intent(in) :: periodic_x, periodic_y, periodic_z

    call neko_log%section('Initializing FST (Forcing)')

    call this%init_common(xstart, xend,xstart,xend,x_delta_rise, x_delta_fall, ystart, &
         yend, ystart, yend, y_delta_rise, y_delta_fall, fringe_max, t_start, t_end, &
         periodic_x, periodic_y, periodic_z)

    call neko_log%end_section('')

    this%is_forcing = .true.
    this%is_bc = .false.

  end subroutine FST_init_forcing

  !> Initialize the FST to use as a boundary condition
  subroutine FST_init_bc(this, &
       xmin, xmax, xstart, xend, x_delta_rise, x_delta_fall, &
       ymin, ymax, ystart, yend, y_delta_rise, y_delta_fall, &
       t_start, t_end, &
       periodic_x, periodic_y, periodic_z)

    class(FST_t), intent(inout) :: this
    real(kind=rp), intent(in) :: xstart, xend, xmin, xmax
    real(kind=rp), intent(in) :: ystart, yend, ymin, ymax
    real(kind=rp), intent(in) :: x_delta_rise
    real(kind=rp), intent(in) :: x_delta_fall
    real(kind=rp), intent(in) :: y_delta_rise
    real(kind=rp), intent(in) :: y_delta_fall
    real(kind=xp), intent(in) :: t_start
    real(kind=xp), intent(in) :: t_end
    logical, intent(in) :: periodic_x, periodic_y, periodic_z

    call neko_log%section('Initializing FST (BC)')

    call this%init_common(xmin, xmax, xstart, xend, x_delta_rise, &
            x_delta_fall, ymin, ymax, ystart, yend, y_delta_rise, &
            y_delta_fall, 1.0_xp, t_start, t_end, &
            periodic_x, periodic_y, periodic_z)

    call neko_log%end_section('')

    this%is_forcing = .false.
    this%is_bc = .true.

  end subroutine FST_init_bc
  
  !> Initialize all the stuff from files generated by the FST code:
  !! - fst_spectrum.csv
  !! - sphere.dat
  !! - bb.txt
  !!
  !! NOTE: The structure of fst_spectrum.csv should be:
  !!  shellno, kx, ky, kz, amp, u_hat_pn(1), u_hat_pn(2), u_hat_pn(3)
  subroutine FST_read_from_files(this)
    class(FST_t), intent(inout) :: this

    integer :: unit, ios, num_columns, num_lines, n_modes_total, i, np_eff, &
         ierr, prev_shell, idx_shell_amp
    character(len=1) :: delimiter
    character(len=1024) :: line
    character(len=20) :: keyword
    real(kind=xp) :: tmp

    delimiter = ','

    if (pe_rank .eq. 0) then

       !
       ! Read sphere.dat to get number of spheres
       !
       call neko_log%message("Reading sphere.dat")
       open(file="sphere.dat", unit=unit, status="old", action="read", &
            iostat=ios)
       if (ios /= 0) then
          call neko_error("Error opening sphere.dat")
       end if

       read(unit,*) line
       read(unit,*) keyword, this%nshells
       print *, this%nshells

       close(unit)

       !
       ! Read FST spectrum, count # of lines to allocate all the arrays
       !
       open(file="fst_spectrum.csv", unit=unit, status="old", action="read", &
            iostat=ios)
       write (*,*) "Reading fst_spectrum"
       if (ios /= 0) then
          call neko_error("Error opening fst_spectrum.csv")
       end if

       num_columns = 1
       num_lines = 0

       ! Read the file line by line
       do
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit

          ! If it's the first line, count the columns
          if (num_columns .eq. 1) then

             ! Count the number of delimiters in the line
             do i = 1, len_trim(line)
                if (line(i:i) == delimiter) then
                   num_columns = num_columns + 1
                end if
             end do

          end if ! if num_columns .eq. 1

          num_lines = num_lines + 1
       end do
       close(unit)

       ! NOTE: this requirement on the number of columns is a bit hardcoded..
       !        AND different from the original implementation. Hence the
       !        extra verbose error message.
       if (num_columns .ne. 8) then
          if (pe_rank .eq. 0) write (*,*) " ***ERROR READING fst_spectrum.csv***"
          if (pe_rank .eq. 0) write (*,*) "fst_spectrum should have 8 cols."
          if (pe_rank .eq. 0) write (*,*) "shell,kx,ky,kz,amp,u_hat_pn1,u_hat_pn2,u_hat_pn3"
          call neko_error("fst_spectrum.csv should have 8 columns")
       end if

       ! NOTE: The total number of modes in the file can sometimes be
       ! different from the theoretical value 2*N_shells*N_points_per_shell.
       ! This is because some modes are removed in the process.
       this%n_modes_total = num_lines - 1 ! Remove the header
       np_eff = this%n_modes_total / this%nshells
       print *, this%n_modes_total, np_eff

    end if ! pe_rank .eq. 0

    call MPI_Bcast(this%n_modes_total, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
    call MPI_Bcast(this%nshells, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

    call MPI_Barrier(NEKO_COMM, ierr)
    print *, pe_rank, this%n_modes_total, this%nshells

    !
    ! Allocate all the relevant arrays:
    ! shell,  kx, ky, kz, amplitudes, u_hat_pn, v_hat_pn, w_hat_pn
    !
    allocate(this%shell(this%n_modes_total))
    allocate(this%k_x(this%n_modes_total))
    allocate(this%k_y(this%n_modes_total))
    allocate(this%k_z(this%n_modes_total))
    allocate(this%shell_amp(this%nshells))
    allocate(this%random_vectors(this%n_modes_total,3))
    allocate(this%phase_shifts(this%n_modes_total))

    if (pe_rank .eq. 0) then
       
       !
       ! Now read fst_spectrum again and populate all the arrays
       !
       open(file="fst_spectrum.csv", unit=unit, status="old", action="read", &
            iostat=ios)
       if (ios /= 0) then
          call neko_error("Error opening fst_spectrum.csv")
       end if

       read(unit,*) line! read the header

       ! Read the file line by line
       do i = 1, this%n_modes_total

          !if (i//this%nshells .eq. 0)
          read(unit,*) this%shell(i), this%k_x(i), &
               this%k_y(i), this%k_z(i), this%shell_amp((i-1)/np_eff+1), &
               this%random_vectors(i,1), this%random_vectors(i,2), &
               this%random_vectors(i,3)
          !print *, "<<0>>",this%shell(i), this%k_x(i), this%k_y(i), this%k_z(i), this%shell_amp((i-1)/np_eff+1), &
          !     this%random_vectors(i,1), this%random_vectors(i,2), this%random_vectors(i,3)
       end do
       close(unit)

       !
       ! Read the phase shifts in bb.txt
       !
       open(file="bb.txt", unit=unit, status="old", action="read", iostat=ios)
       do i = 1, this%n_modes_total
          read(unit,*) this%phase_shifts(i), tmp
          !print *, this%phase_shifts(i)
       end do

    end if ! pe_rank .eq. 0

    call MPI_Bcast(this%k_x, this%n_modes_total, MPI_EXTRA_PRECISION, 0, &
         NEKO_COMM, ierr)
    call MPI_Bcast(this%k_y, this%n_modes_total, MPI_EXTRA_PRECISION, 0, &
         NEKO_COMM, ierr)
    call MPI_Bcast(this%k_z, this%n_modes_total, MPI_EXTRA_PRECISION, 0, &
         NEKO_COMM, ierr)
    call MPI_Bcast(this%shell, this%n_modes_total, MPI_INTEGER, 0, NEKO_COMM, &
         ierr)
    call MPI_Bcast(this%shell_amp, this%nshells, MPI_EXTRA_PRECISION, 0, &
         NEKO_COMM, ierr)
    call MPI_Bcast(this%random_vectors, this%n_modes_total*3, &
         MPI_EXTRA_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(this%phase_shifts, this%n_modes_total, MPI_EXTRA_PRECISION, &
         0, NEKO_COMM, ierr)
    call MPI_Barrier(NEKO_COMM, ierr)

    !if (pe_rank .eq. 1) then
    !   do i = 1, this%n_modes_total
    !      !print *, "<<1>>", this%shell(i), this%k_x(i), this%k_y(i), this%k_z(i), &
    !      !this%shell_amp((i-1)/78 +1), &
    !      !     this%random_vectors(i,1), this%random_vectors(i,2), this%random_vectors(i,3)
    !   end do
    !end if

  end subroutine FST_read_from_files


  !! Free parameters in global params
  subroutine FST_free_params(this)
    class(FST_t), intent(inout) :: this

    if(allocated(this%fringe_space)) deallocate(this%fringe_space)
    if(allocated(this%phi_0)) deallocate(this%phi_0)
    if(allocated(this%k_x)) deallocate(this%k_x)
    if(allocated(this%k_y)) deallocate(this%k_y)
    if(allocated(this%k_z)) deallocate(this%k_z)
    if(allocated(this%shell)) deallocate(this%shell)
    if(allocated(this%shell_amp)) deallocate(this%shell_amp)
    if(allocated(this%phi_0)) deallocate(this%phi_0)

    if (allocated(this%u_baseflow)) deallocate(this%u_baseflow)
    if (allocated(this%v_baseflow)) deallocate(this%v_baseflow)
    if (allocated(this%w_baseflow)) deallocate(this%w_baseflow)

  end subroutine FST_free_params

  subroutine FST_print_params(this)
    class(FST_t) :: this

    character(len=LOG_SIZE) :: log_buf
    
    if (pe_rank .ne. 0) return

    write(log_buf, '(A ,F15.7)') "xstart      ", this%xstart
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "xend        ", this%xend
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "xmin        ", this%xmin
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "xmax        ", this%xmax
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "ystart      ", this%ystart
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "yend        ", this%yend
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "ymin        ", this%ymin
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "ymax        ", this%ymax
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "fringe_max  ", this%fringe_max
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "x_delta_rise", this%x_delta_rise
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "x_delta_fall", this%x_delta_fall
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "y_delta_rise", this%y_delta_rise
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "y_delta_fall", this%y_delta_fall
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "t_start     ", this%t_start
    call neko_log%message(log_buf)
    write(log_buf, '(A ,F15.7)') "t_end       ", this%t_end
    call neko_log%message(log_buf)

  end subroutine FST_print_params

  !> Apply a specific baseflow in our region, from a boundary mask.
  !! If we specify v or w it takes the norm.
  subroutine FST_apply_baseflow(this, mask, n, u, v, w)
    class(FST_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w
    integer, intent(in) :: mask(0:n)
    integer, intent(in) :: n

    type(c_ptr) :: mask_d
    integer :: i, idx

    if (allocated(this%u_baseflow)) deallocate(this%u_baseflow)
    if (allocated(this%v_baseflow)) deallocate(this%v_baseflow)
    if (allocated(this%w_baseflow)) deallocate(this%w_baseflow)

    allocate(this%u_baseflow(n))
    allocate(this%v_baseflow(n))
    allocate(this%w_baseflow(n))

    if (neko_bcknd_device .eq. 1) then
       call device_map(this%u_baseflow, this%u_baseflow_d, n)
       call device_map(this%v_baseflow, this%v_baseflow_d, n)
       call device_map(this%w_baseflow, this%w_baseflow_d, n)

       mask_d = device_get_ptr(mask)

       call device_masked_gather_copy_0(this%u_baseflow_d, u%x_d, mask_d, u%dof%size(), n)
       call device_masked_gather_copy_0(this%v_baseflow_d, v%x_d, mask_d, u%dof%size(), n)
       call device_masked_gather_copy_0(this%w_baseflow_d, w%x_d, mask_d, u%dof%size(), n)
    else
       call masked_gather_copy(this%u_baseflow, u%x, mask, u%dof%size(), n)
       call masked_gather_copy(this%v_baseflow, v%x, mask, u%dof%size(), n)
       call masked_gather_copy(this%w_baseflow, w%x, mask, u%dof%size(), n)
    end if
         
  end subroutine FST_apply_baseflow

  !> Common function for generation
  subroutine FST_generate_common(this, coef)
    class(FST_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    integer :: ierr

    !
    ! First, generate everything as usual
    !
    call neko_log%section ('Generating FST')
    call make_turbu(coef, this%periodic_x, this%periodic_y, this%periodic_z)
    call neko_log%end_section('Done --> Generating FST')

    !
    ! Next, use the arrays defined

  end subroutine FST_generate_common

  !> Generate FST for forcing
  subroutine FST_generate_forcing(this, coef, zone, u, v, w)
    class(FST_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    class(point_zone_t), intent(in) :: zone
    type(field_t), intent(in) :: u, v, w

    real(kind=xp) :: x, y, z
    integer :: ierr, i, idx

    !! Do the general generation
    !call this%generate_common(coef)

    !!
    !! Copy the baseflow in the zone
    !!
    !call this%apply_baseflow(zone%mask, zone%size, u, v, w)

    !! Generate the fringe in space
    !allocate(this%fringe_space(zone%size))

    !! Initialize the fringe in space
    !do idx = 1, zone%size
    !   i = zone%mask(idx)
    !   x = coef%dof%x(i,1,1,1)
    !   y = coef%dof%y(i,1,1,1)
    !   z = coef%dof%z(i,1,1,1)

    !   if ( x .lt. this%xmin .or. x .gt. this%xmax .or. &
    !        y .lt. this%ymin .or. y .gt. this%ymax ) then
    !      this%fringe_space(idx) = 0.0_xp
    !   else
    !      this%fringe_space(idx) = fringe(x, y, this)
    !   end if

    !end do

  end subroutine FST_generate_forcing

  !> Do the generation for BC.
  subroutine FST_generate_bc(this, coef, bc_mask, n, u, v, w, regen, Uinf)
    class(FST_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: bc_mask(0:n)
    integer, intent(in) :: n
    type(field_t), intent(in) :: u, v, w
    logical, intent(in) :: regen
    real(kind=xp), optional :: Uinf

    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: x, y, z
    integer :: ierr, i, idx, m, j
    
    !
    ! Generate everything from scratch. This will create the files
    ! bb.txt, sphere.dat and fst_spectrum.csv
    !
    if (regen) then

       call this%generate_common(coef)
       this%Uinf = glb_uinf

    else

       if (.not. present(Uinf)) then
          call neko_error("Provide a value for Uinf if you init from file!")
       else 
          this%Uinf = Uinf
       end if

    end if

    !
    ! Read data from files (will not use anything from global_params!)
    !
    call this%read_from_files()

    !
    ! Print some diagnostics just in case
    !
    write (log_buf, '(A,F15.7)') "(FST) Uinf set to ", this%Uinf
    call neko_log%message(log_buf)



    !
    ! Apply baseflow in the bc zone
    !
    call this%apply_baseflow(bc_mask, bc_mask(0), u, v, w)

    !
    ! Initialize the fringe in space
    !
    allocate(this%fringe_space(bc_mask(0)))
    do idx = 1, size(this%fringe_space)

       i = bc_mask(idx)
       x = coef%dof%x(i,1,1,1)
       y = coef%dof%y(i,1,1,1)
       z = coef%dof%z(i,1,1,1)

       if ( z .lt. this%xmin .or. z .gt. this%xmax .or. &
            y .lt. this%ymin .or. y .gt. this%ymax ) then
          this%fringe_space(idx) = 0.0_xp
       else
          this%fringe_space(idx) = fringe(z, y, this)
       end if

    end do

    !
    ! Precompute time-independent term
    ! 
    allocate(this%phi_0(this%n_modes_total, bc_mask(0)))

    do j = 1, bc_mask(0)
      x = coef%dof%x(bc_mask(j), 1,1,1)
      y = coef%dof%y(bc_mask(j), 1,1,1)
      z = coef%dof%z(bc_mask(j), 1,1,1)
      do m = 1, this%n_modes_total
        this%phi_0(m,j) = this%k_x(m)*x + this%k_y(m)*y + this%k_z(m)*z + &
             this%phase_shifts(m)
      end do
    end do

    !
    ! Copy the wavenumbers in x direction
    !
!!$    allocate(this%k_x(fst_modes))
!!$
!!$    do m = 1, fst_modes
!!$       this%k_x(m) = k_num_all(m,1)
!!$    end do
!!$
    !
    ! Copy a proper version of u_hat_pn with correct size
    !
!!$    allocate(this%random_vectors(this%n_modes_total, 3))
!!$
!!$    do m = 1, this%n_modes_total
!!$       do j =  1, 3
!!$          this%random_vectors(m,j) = u_hat_pn(m,j)
!!$       end do
!!$    end do
!!$
    ! Copy everything to device and map with relevant device array pointers
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%fringe_space, this%fringe_space_d, bc_mask(0))
       call device_memcpy(this%fringe_space, this%fringe_space_d, bc_mask(0), &
            HOST_TO_DEVICE, .false.)

       call device_map(this%phi_0, this%phi_0_d, this%n_modes_total*bc_mask(0))
       call device_memcpy(this%phi_0, this%phi_0_d, &
            this%n_modes_total*bc_mask(0), HOST_TO_DEVICE, .false.)

       call device_map(this%random_vectors, this%random_vectors_d, &
            this%n_modes_total*3)
       call device_memcpy(this%random_vectors, this%random_vectors_d, &
            this%n_modes_total*3, HOST_TO_DEVICE, .false.)

       call device_map(this%shell, this%shell_d, this%n_modes_total)
       call device_memcpy(this%shell, this%shell_d, this%n_modes_total, &
            HOST_TO_DEVICE, .false.)

       call device_map(this%shell_amp, this%shell_amp_d, this%nshells)
       call device_memcpy(this%shell_amp, this%shell_amp_d, this%nshells, &
            HOST_TO_DEVICE, .false.)

       call device_map(this%k_x, this%k_x_d, this%n_modes_total)
       call device_memcpy(this%k_x, this%k_x_d, this%n_modes_total, &
            HOST_TO_DEVICE, .true.)
    end if

  end subroutine FST_generate_bc

  ! Forcing to be performed on entire domain, on a local element ix
  ! Final values of the forcing are to be applied
  ! on f%u, f%v and f%w
  subroutine FST_forcing_zone(this, zone, &
       x, y, z, t, &
       u, v, w, &
       fu, fv, fw)

    class(FST_t), intent(in) :: this
    class(point_zone_t), intent(in) :: zone
    real(kind=rp), intent(in), dimension(:,:,:,:) :: x, y, z, u, v, w
    real(kind=rp), intent(in) :: t
    real(kind=rp), intent(inout), dimension(:,:,:,:) :: fu, fv, fw

    integer :: idx, l, m, i, shellno
    integer, parameter :: gdim = 3
    real(kind=rp) :: phase_shft, phi, amp, pert, vel_mag
    real(kind=rp) :: rand_vec(gdim)
    real(kind=rp) :: fringe_time

!!$    fringe_time = time_ramp(t, this%t_end, this%t_start)
!!$
!!$    ! Loop on all points in the point zone
!!$    do idx = 1, zone%size
!!$       i = zone%mask(idx)
!!$
!!$       !> This vector will contain the sum of all Fourier modes
!!$       rand_vec = 0.0_rp
!!$
!!$       ! Sum all sin modes for each gll point
!!$       do m = 1, this%n_modes_total
!!$          phase_shft = bb(m,1)
!!$
!!$          ! k_x(x - U*t) + ky*y + kz*z - random_phase[-pi, pi]
!!$          ! = k_x*x + ky*y + kz*z - k_x*U*t - random_phase
!!$          phi = k_num_all(m,1) * (x(i,1,1,1) - glb_uinf*t) + &
!!$               k_num_all(m,2) *  y(i,1,1,1) + &
!!$               k_num_all(m,3) *  z(i,1,1,1) + &
!!$               phase_shft
!!$
!!$          shellno = shell(m)
!!$          pert = shell_amp(shellno)*sin(phi)
!!$
!!$          rand_vec(1) = rand_vec(1) + u_hat_pn(m,1)*pert
!!$          rand_vec(2) = rand_vec(2) + u_hat_pn(m,2)*pert
!!$          rand_vec(3) = rand_vec(3) + u_hat_pn(m,3)*pert
!!$       enddo
!!$
!!$       fu(i,1,1,1) = fringe_time*this%fringe_space(idx)*( &
!!$            this%u_baseflow(idx) + rand_vec(1) - u(i,1,1,1))
!!$
!!$       fv(i,1,1,1) = fringe_time*this%fringe_space(idx)*( &
!!$            this%v_baseflow(idx) + rand_vec(2) - v(i,1,1,1))
!!$
!!$       fw(i,1,1,1) = fringe_time*this%fringe_space(idx)*( &
!!$            this%w_baseflow(idx) + rand_vec(3) - w(i,1,1,1))
!!$    end do

  end subroutine FST_forcing_zone

  ! Apply FST as a boundary condition based on the bc mask
  ! Assumes that u,v,w have the same bc masks
  subroutine FST_apply_BC(this, bc_mask, n, &
       t, &
       u_bc, v_bc, w_bc, angleXY, on_host)

    class(FST_t), intent(in) :: this
    integer, intent(in) :: n ! size of the bc mask
    integer, intent(in) :: bc_mask(0:n)
    real(kind=rp), intent(in) :: t
    real(kind=rp), intent(inout), dimension(:,:,:,:) :: u_bc, v_bc, w_bc
    real(kind=xp), intent(in) :: angleXY
    logical, intent(in) :: on_host

!!$    integer :: idx, l, m, i, shellno
!!$    integer, parameter :: gdim = 3
!!$    real(kind=rp) :: phase_shft, phi, amp, pert, urand, vrand, wrand
!!$    real(kind=rp) :: rand_vec(gdim), vel_mag, phi_t

    real(kind=xp) :: fringe_time, cosa, sina

    fringe_time = time_ramp(t, this%t_end, this%t_start)
    cosa = cos(angleXY)
    sina = sin(angleXY)

    call fst_bc_compute(t, this%Uinf, u_bc, v_bc, w_bc, bc_mask, n, &
       this%u_baseflow, this%v_baseflow, this%w_baseflow, &
       this%k_x, this%n_modes_total, this%phi_0, this%shell, this%shell_amp, &
       this%random_vectors, angleXY, fringe_time, this%fringe_space, on_host)

!!$    phi_t = glb_uinf*t
!!$    ! Loop on all points in the point zone
!!$    do idx = 1, bc_mask(0)
!!$
!!$       i = bc_mask(idx)
!!$
!!$       !> This vector will contain the sum of all Fourier modes
!!$       rand_vec = 0.0_rp
!!$
!!$       ! Sum all sin modes for each gll point
!!$       do m = 1, this%n_modes_total
!!$
!!$          ! Random phase shifts
!!$          !phase_shft = bb(m,1)
!!$
!!$          ! k_x(x - U*t) + ky*y + kz*z + phi
!!$          ! = k_x*x + ky*y + kz*z - k_x*U*t + phi
!!$          !phi = k_num_all(m,1) * (x(i,1,1,1) - glb_uinf*t) + &
!!$          !     k_num_all(m,2) *  y(i,1,1,1) + &
!!$          !     k_num_all(m,3) *  z(i,1,1,1) + &
!!$          !     phase_shft
!!$
!!$          phi = this%phi_0(m,idx) - this%k_x(m)*phi_t
!!$
!!$          shellno = shell(m)
!!$
!!$          pert = shell_amp(shellno)*sin(phi)
!!$
!!$          rand_vec(1) = rand_vec(1) + u_hat_pn(m,1)*pert
!!$          rand_vec(2) = rand_vec(2) + u_hat_pn(m,2)*pert
!!$          rand_vec(3) = rand_vec(3) + u_hat_pn(m,3)*pert
!!$
!!$       enddo
!!$
!!$       ! Project the rand_vec if there is a rotation
!!$       urand = rand_vec(1)*cosa - rand_vec(2)*sina
!!$       vrand = rand_vec(1)*sina + rand_vec(2)*cosa
!!$       wrand = rand_vec(3)
!!$
!!$       u_bc(i,1,1,1) = this%u_baseflow(idx) + &
!!$            fringe_time*this%fringe_space(idx)*urand
!!$       v_bc(i,1,1,1) = this%v_baseflow(idx) + &
!!$            fringe_time*this%fringe_space(idx)*vrand
!!$       w_bc(i,1,1,1) = this%w_baseflow(idx) + &
!!$            fringe_time*this%fringe_space(idx)*wrand
!!$
!!$    end do

  end subroutine FST_apply_BC

  ! ! Fringe function as defined in original FST.
  ! function fringe(x, f) result(y)
  !   real(kind=rp), intent(in) :: x
  !   type(FST_t) :: f
  !   real(kind=rp) :: y, S

  !   if(x.le.f%xstart) then
  !      S=0
  !   else if (x.ge.f%xmax) then
  !      S=0
  !   else if (x.ge.f%xend) then
  !      S=1
  !   else
  !      S=1/(1+exp(1/(x-f%xend)+1/(x-f%xstart)))
  !   endif

  !   y = f%fringe_max * (S*(x-f%xstart)/(f%delta_rise)-S*((x-f%xend)/f%delta_fall+1))

  ! end function fringe

  !
  ! Linear ramp in time
  function time_ramp(t, t_end, t_start) result(ramp)
    real(kind=rp), intent(in) :: t
    real(kind=xp), intent(in) :: t_end
    real(kind=xp), intent(in) :: t_start

    real(kind=xp) :: ramp

    if (t .le. t_start) then
       ramp = 0.0_xp
    else if (t .lt. t_end) then 
       ramp = (t - t_start)/(t_end - t_start)
    else
       ramp = 1.0_xp
    end if

  end function time_ramp

  !
  ! Fringe function as described in Schlatter (2001), extended to take y bounds into account too
  !
  !   Here is what the fringe looks like, except the ramp-up is not linear
  !   but exponential (see function S below)
  !
  ! fringe_max      ________
  !                /        \
  !               /          \
  ! 0.0 _________/            \_______
  !
  !   xmin   xstart         xend   xmax
  !
  ! The ramp-up 
  ! after xstart is of length delta_rise. The same applies for the ramp down
  ! and the distance between xmax and xend. If y is specified then it computes
  ! a 2D fringe
  function fringe(x, y, f) result(fr)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in), optional :: y
    type(FST_t), intent(in) :: f
    real(kind=xp) :: fr
    integer :: i
    character :: a

    fr = f%fringe_max * ( S((x-f%xstart)/f%x_delta_rise) - S((x-f%xend)/f%x_delta_fall + 1.0_rp) )

    if (present(y)) then
       fr = fr * ( S((y-f%ystart)/f%y_delta_rise) - S((y-f%yend)/f%y_delta_fall + 1._rp) )
    end if

  end function fringe

  ! Smooth step function, 0 if x <= 0, 1 if x >= 1, 1/exp(1/(x-1) + 1/x) between 0 and 1
  function S(x) result(y)
    real(kind=rp), intent(in) :: x
    real(kind=rp)             :: y

    if ( x.le.0._xp ) then
       y = 0._xp
    else if ( x.ge.1._xp ) then
       y = 1._xp
    else
       y = 1._xp / (1._xp + exp( 1._xp/(x-1._xp) + 1._xp/x))
    end if

  end function S


end module FST
