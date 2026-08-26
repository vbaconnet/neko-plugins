! * Free stream turbulence implementation for neko
! * Author: Victor Baconnet (baconnet@kth.se)
! * Based in original implementation by Elektra Kluesberg, Prabal Negi
! 	and Philipp Schlatter.
! * For FST theory see Schlatter (2001)

module FST

  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  !use global_params
  use fst_operator, only: fst_bc_compute
  use turbu, only : make_turbu
  use fst_utils, only : print_param
  use logger, only: LOG_SIZE, neko_log
  use field, only: field_t
  use coefs, only: coef_t
  use num_types, only: rp
  use utils, only: neko_error
  use point_zone, only: point_zone_t
  use comm, only: pe_rank
  use math, only: masked_gather_copy_0, pi, glmin, glmax
  use device_math, only: device_masked_gather_copy_0
  use device, only: device_map, device_memcpy, HOST_TO_DEVICE, device_unmap
  use mpi_f08, only: MPI_IN_PLACE, MPI_MAX, MPI_MIN, MPI_INTEGER, MPI_Bcast, &
      MPI_Allreduce
  use comm, only: MPI_REAL_PRECISION, NEKO_COMM
  use device, only: device_get_ptr
  use neko_config, only: NEKO_BCKND_DEVICE
  implicit none


  type, public :: FST_t

     ! Seed for random number generator
     integer :: seed = -143

     ! periodic directions
     logical :: periodic_x
     real(kind=rp) :: Lx
     logical :: periodic_y
     real(kind=rp) :: Ly
     logical :: periodic_z
     real(kind=rp) :: Lz

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
     real(kind=rp) :: fringe_max

     !> Final ramp time
     real(kind=rp) :: t_end
     real(kind=rp) :: t_start

     logical :: is_forcing
     logical :: is_bc

     ! ------- 
     ! FST generation parameters
     !> Free-stream velocity
     real(kind=rp) :: Uinf
     !> Turbulence intensity
     real(kind=rp) :: Tu
     !> Turbulent length scale
     real(kind=rp) :: L
     !> Number of shells
     integer :: n_shells
     !> Max number of points per shell
     integer :: n_max_pts_per_shell
     !> Start and end bounds of the total wavenumber range
     real(kind=rp) :: k_start, k_end
     !> Effective number of points per shell, usually equal to but sometimes 
     !! lower than n_max_pts_per_shell
     integer :: n_eff_pts_per_shell = 0
     !> Total number of modes (2*n_shells*n_max_eff_per_shell)
     integer :: n_modes = 0
     
     !> Random, divergence-free unit vectors
     real(kind=rp), allocatable :: random_vectors(:,:) ! u_hat_pn but reshaped
     type(c_ptr) :: random_vectors_d = C_NULL_PTR
     !> Map to the total wavenumber id/shell
     integer, allocatable :: shell(:)
     type(c_ptr) :: shell_d = C_NULL_PTR
     !> Amplitude of mode on a given shell
     real(kind=rp), allocatable :: shell_amp(:)
     type(c_ptr) :: shell_amp_d = C_NULL_PTR

     !> Wavenumbers in the x,y,z direction
     real(kind=rp), allocatable :: k_x(:)
     real(kind=rp), allocatable :: k_y(:)
     real(kind=rp), allocatable :: k_z(:)
     type(c_ptr) :: k_x_d = C_NULL_PTR

     real(kind=rp), allocatable :: phase_shifts(:)
     ! -------

     !> Fringe in space
     real(kind=rp), allocatable :: fringe_space(:)
     type(c_ptr) :: fringe_space_d = C_NULL_PTR

     !> Baseflows, if applying on a non-uniform inflow
     real(kind=rp), allocatable :: u_baseflow(:)
     real(kind=rp), allocatable :: v_baseflow(:)
     real(kind=rp), allocatable :: w_baseflow(:)
     type(c_ptr) :: u_baseflow_d = C_NULL_PTR
     type(c_ptr) :: v_baseflow_d = C_NULL_PTR
     type(c_ptr) :: w_baseflow_d = C_NULL_PTR

     !> Variable that is precomputed to save some time
     real(kind=rp), allocatable :: phi_0(:,:)
     type(c_ptr) :: phi_0_d = C_NULL_PTR

   contains

     ! ======== Init/Free procedures
     procedure, pass(this) :: init_common => FST_init_common
     procedure, pass(this) :: init_bc => FST_init_bc
     !procedure, pass(this) :: init_forcing => FST_init_forcing
     procedure, pass(this) :: unmap => FST_unmap
     procedure, pass(this) :: free => FST_free
     ! =========================================================================
     procedure, pass(this) :: apply_baseflow => FST_apply_baseflow
     procedure, pass(this) :: apply_baseflow_0 => FST_apply_baseflow_0
     ! =========================================================================
     ! ======== Generate FST
     procedure, pass(this) :: generate_common => FST_generate_common
     !procedure, pass(this) :: generate_forcing => FST_generate_forcing
     procedure, pass(this) :: generate_bc => FST_generate_bc
     ! =========================================================================
     ! ======= Apply FST forcing/BC
     !   procedure, pass(this) :: apply_forcing => FST_forcing_zone
     procedure, pass(this) :: apply_BC => FST_apply_BC
     ! ========================================================================
     procedure, pass(this) :: validate => FST_validate
     procedure, pass(this) :: print => FST_print_params
  end type FST_t

contains

  !> Initialize all parameters
  subroutine FST_init_common(this, &
       Uinf, Tu, L, k_start, k_end, Nshells, Npmax, &
       periodic_x, periodic_y, periodic_z, &
       xmin, xmax, xstart, xend, x_delta_rise, x_delta_fall, &
       ymin, ymax, ystart, yend, y_delta_rise, y_delta_fall, &
       fringe_max, &
       t_start, t_ramp, &
       seed)
    class(FST_t), intent(inout) :: this
    real(kind=rp), intent(in) :: Uinf, Tu, L, k_start, k_end
    integer, intent(in) :: Nshells, Npmax
    real(kind=rp), intent(in) :: xmin, xmax, xstart, xend
    real(kind=rp), intent(in) :: ymin, ymax, ystart, yend
    real(kind=rp), intent(in) :: x_delta_rise, x_delta_fall
    real(kind=rp), intent(in) :: y_delta_rise, y_delta_fall
    real(kind=rp), intent(in) :: fringe_max
    real(kind=rp), intent(in) :: t_start
    real(kind=rp), intent(in) :: t_ramp
    logical, intent(in) :: periodic_x, periodic_y, periodic_z
    integer, intent(inout), optional :: seed

    integer :: seed_ = -143

    call this%free()

    if (present(seed)) seed_ = seed
    this%seed = seed_

    this%Uinf = Uinf
    this%Tu = Tu
    this%L = L
    this%k_start = k_start
    this%k_end = k_end

    this%n_shells = Nshells
    allocate(this%shell_amp(Nshells))

    this%n_max_pts_per_shell = Npmax

    this%periodic_x = periodic_x
    this%periodic_y = periodic_y
    this%periodic_z = periodic_z

    this%xstart = xstart
    this%xend = xend
    this%ystart = ystart
    this%yend = yend

    this%xmin = xmin
    this%xmax = xmax
    this%ymin = ymin
    this%ymax = ymax

    this%fringe_max = fringe_max
    this%x_delta_rise = x_delta_rise!0.002
    this%x_delta_fall = x_delta_fall!0.002
    this%y_delta_rise = y_delta_rise!0.002
    this%y_delta_fall = y_delta_fall!0.002
    this%t_start = t_start
    this%t_end = t_ramp
    if (t_start .lt. 0.0_rp .or. t_ramp .lt. 0.0_rp .or. t_ramp .lt. t_start) then
      call neko_error("t_start or t_ramp is invalid!")
    end if

  end subroutine FST_init_common


  !> Initialize the FST to use with forcing.
!   subroutine FST_init_forcing(this, &
!        xstart, xend, x_delta_rise, x_delta_fall, &
!        ystart, yend, y_delta_rise, y_delta_fall, &
!        fringe_max, &
!        t_start, t_end, &
!        periodic_x, periodic_y, periodic_z)

!     class(FST_t), intent(inout) :: this
!     real(kind=rp), intent(in) :: xstart
!     real(kind=rp), intent(in) :: xend
!     real(kind=rp), intent(in) :: ystart
!     real(kind=rp), intent(in) :: yend
!     real(kind=rp), intent(in) :: x_delta_rise
!     real(kind=rp), intent(in) :: x_delta_fall
!     real(kind=rp), intent(in) :: y_delta_rise
!     real(kind=rp), intent(in) :: y_delta_fall
!     real(kind=rp), intent(in) :: fringe_max
!     real(kind=rp), intent(in) :: t_start
!     real(kind=rp), intent(in) :: t_end
!     logical, intent(in) :: periodic_x, periodic_y, periodic_z

!     call neko_log%section('Initializing FST')

!     call this%init_common(xstart, xend,xstart,xend,x_delta_rise, x_delta_fall, ystart, &
!          yend, ystart, yend, y_delta_rise, y_delta_fall, fringe_max, t_start, t_end, &
!          periodic_x, periodic_y, periodic_z)

!     call this%print() ! show parameters
!     call neko_log%end_section('Done --> Intializing FST')

!     this%is_forcing = .true.
!     this%is_bc = .false.

!   end subroutine FST_init_forcing

  !> Initialize the FST to use as a boundary condition
  subroutine FST_init_bc(this, &
       Uinf, Tu, L, k_start, k_end, Nshells, Npmax, &
       periodic_x, periodic_y, periodic_z, &
       xmin, xmax, xstart, xend, x_delta_rise, x_delta_fall, &
       ymin, ymax, ystart, yend, y_delta_rise, y_delta_fall, &
       t_start, t_end, &
       seed)

    class(FST_t), intent(inout) :: this
    real(kind=rp), intent(in) :: Uinf, Tu, L, k_start, k_end
    integer, intent(in) :: Nshells, Npmax
    real(kind=rp), intent(in) :: xmin, xmax, xstart, xend
    real(kind=rp), intent(in) :: ymin, ymax, ystart, yend
    real(kind=rp), intent(in) :: x_delta_rise, x_delta_fall
    real(kind=rp), intent(in) :: y_delta_rise, y_delta_fall
    real(kind=rp), intent(in) :: t_start
    real(kind=rp), intent(in) :: t_end
    logical, intent(in) :: periodic_x, periodic_y, periodic_z
    integer, intent(inout), optional :: seed

    call neko_log%section('Initializing FST')

    call this%init_common(Uinf, Tu, L, k_start, k_end, Nshells, Npmax, &
       periodic_x, periodic_y, periodic_z, &
       xmin, xmax, xstart, xend, x_delta_rise, x_delta_fall, &
       ymin, ymax, ystart, yend, y_delta_rise, y_delta_fall, &
       1.0_rp, &
       t_start, t_end, &
       seed)

    call this%print() ! show parameters

    call neko_log%end_section()

    this%is_forcing = .false.
    this%is_bc = .true.

  end subroutine FST_init_bc

  !> Unmap device arrays
  subroutine FST_unmap(this)
    class(FST_t), intent(inout) :: this

    if (allocated(this%random_vectors)) call device_unmap(this%random_vectors, this%random_vectors_d)
    if (allocated(this%k_x)) call device_unmap(this%k_x, this%k_x_d)
    if (allocated(this%fringe_space)) call device_unmap(this%fringe_space, this%fringe_space_d)
    if (allocated(this%shell)) call device_unmap(this%shell, this%shell_d)
    if (allocated(this%shell_amp)) call device_unmap(this%shell_amp, this%shell_amp_d)
    if (allocated(this%u_baseflow)) call device_unmap(this%u_baseflow, this%u_baseflow_d)
    if (allocated(this%v_baseflow)) call device_unmap(this%v_baseflow, this%v_baseflow_d)
    if (allocated(this%w_baseflow)) call device_unmap(this%w_baseflow, this%w_baseflow_d)

  end subroutine FST_unmap

  !! Free parameters in global params
  subroutine FST_free(this)
    class(FST_t), intent(inout) :: this

    if (NEKO_BCKND_DEVICE .eq. 1) call this%unmap()

    if (allocated(this%fringe_space)) deallocate(this%fringe_space)
    if (allocated(this%phi_0)) deallocate(this%phi_0)
    if (allocated(this%u_baseflow)) deallocate(this%u_baseflow)
    if (allocated(this%v_baseflow)) deallocate(this%v_baseflow)
    if (allocated(this%w_baseflow)) deallocate(this%w_baseflow)
    if (allocated(this%random_vectors)) deallocate(this%random_vectors)
    if (allocated(this%phase_shifts)) deallocate(this%phase_shifts)
    if (allocated(this%k_x)) deallocate(this%k_x)
    if (allocated(this%k_y)) deallocate(this%k_y)
    if (allocated(this%k_z)) deallocate(this%k_z)
    if (allocated(this%shell)) deallocate(this%shell)
    if (allocated(this%shell_amp)) deallocate(this%shell_amp)

  end subroutine FST_free

  subroutine FST_print_params(this)
    class(FST_t) :: this

    call neko_log%section("--- FST Generation ---")
    call print_param("Free-stream vel.   ", this%Uinf, fmt='F6.2')
    call print_param("Turb. intensity (%)", this%Tu * 100_rp, fmt='F6.2')
    call print_param("Turb. length scale ", this%L, fmt='F10.6')
    call print_param("k_start            ", this%k_start)
    call print_param("k_end              ", this%k_end)
    call print_param("# of shells        ", this%n_shells)
    call print_param("max # pts per shell", this%n_max_pts_per_shell)
    call print_param("Periodic in x      ", this%periodic_x)
    call print_param("Periodic in y      ", this%periodic_y)
    call print_param("Periodic in z      ", this%periodic_z)
    call print_param("seed               ", this%seed)
    call neko_log%end_section()

    call neko_log%section("--- Fringe ---")
    call neko_log%section("Space")
    if (.not. this%periodic_z) then
      call print_param("xmin           ", this%xmin)
      call print_param("xmax           ", this%xmax)
      call print_param("xstart         ", this%xstart)
      call print_param("xend           ", this%xend)
      call print_param("x_delta_rise   ", this%x_delta_rise)
      call print_param("x_delta_fall   ", this%x_delta_fall)
    else
      call neko_log%message("(periodic in y, no fringe required)")
    end if

    if (.not. this%periodic_y) then
      call print_param("ymin           ", this%ymin)
      call print_param("ymax           ", this%ymax)
      call print_param("ystart         ", this%ystart)
      call print_param("yend           ", this%yend)
      call print_param("y_delta_rise   ", this%y_delta_rise)
      call print_param("y_delta_fall   ", this%y_delta_fall)
    else
      call neko_log%message("(periodic in z, no fringe required)")
    end if
   
    call print_param("Fringe amplitude ", this%fringe_max, fmt='F3.1')
    call neko_log%end_section()

    call neko_log%section("Time")
    call print_param("t_start          ", this%t_start)
    call print_param("t_end            ", this%t_end)
    call neko_log%end_section()

    call neko_log%end_section()

  end subroutine FST_print_params

  !> Apply a specific baseflow in our region, from a boundary mask.
  !! If we specify v or w it takes the norm.
  subroutine FST_apply_baseflow_0(this, mask, n, u, v, w)
    class(FST_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w
    integer, intent(in) :: mask(0:n)
    integer, intent(in) :: n

    type(c_ptr) :: mask_d
    integer :: i, idx

    if (allocated(this%u_baseflow)) call neko_error("Baseflow already allocated!")
    if (allocated(this%v_baseflow)) call neko_error("Baseflow already allocated!")
    if (allocated(this%w_baseflow)) call neko_error("Baseflow already allocated!")

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
       call masked_gather_copy_0(this%u_baseflow, u%x, mask, u%dof%size(), n)
       call masked_gather_copy_0(this%v_baseflow, v%x, mask, u%dof%size(), n)
       call masked_gather_copy_0(this%w_baseflow, w%x, mask, u%dof%size(), n)
    end if

  end subroutine FST_apply_baseflow_0

  !> Apply a specific baseflow in our region, from a boundary mask.
  !! If we specify v or w it takes the norm.
  subroutine FST_apply_baseflow(this, mask, n, u, v, w)
    class(FST_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w
    integer, intent(in) :: mask(n)
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
       call masked_gather_copy_0(this%u_baseflow, u%x, mask, u%dof%size(), n)
       call masked_gather_copy_0(this%v_baseflow, v%x, mask, u%dof%size(), n)
       call masked_gather_copy_0(this%w_baseflow, w%x, mask, u%dof%size(), n)
    end if

  end subroutine FST_apply_baseflow

  !> Common function for generation
  subroutine FST_generate_common(this, path, gdim, Lx, Ly, Lz)
    class(FST_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    integer, intent(in) :: gdim
    real(kind=rp), intent(in), optional :: Lx, Ly, Lz

    integer :: ierr

    call neko_log%section ('Generating FST')

    ! NOTE: the generation is done on rank 0 and the data is broadcast
    ! to al ranks
    call make_turbu(this%phase_shifts, this%random_vectors, &
      this%n_eff_pts_per_shell, this%L, this%Tu, this%Uinf, &
      this%n_max_pts_per_shell, this%n_shells, this%k_start, this%k_end, &
      this%k_x, this%k_y, this%k_z, this%n_modes, this%shell, &
      this%shell_amp, this%periodic_x, this%periodic_y, this%periodic_z, &
      this%seed, path, .true., gdim, Lx, Ly, Lz)
    
    call neko_log%end_section('')

    call this%validate()

  end subroutine FST_generate_common

  subroutine FST_validate(this)
   class(FST_t), intent(in) :: this

   real(kind=rp) :: kmin, kmax
   character(len=LOG_SIZE) :: log_buf
   integer :: ierr

   if (pe_rank .ne. 0) return

   call neko_log%section("FST diagnostics")

   !
   ! Compute min/max wavenumbers
   !
   call neko_log%section("Wavenumbers")

   ! x-direction
   kmin = glmin(abs(this%k_x), this%n_modes)
   kmax = glmax(abs(this%k_x), this%n_modes)
   call print_param("(x) min wavelength", 2.0_rp*pi/kmax, fmt='F10.4')
   call print_param("(x) max wavelength", 2.0_rp*pi/kmin, fmt='F10.4')

   ! y-direction
   kmin = glmin(abs(this%k_y), this%n_modes)
   kmax = glmax(abs(this%k_y), this%n_modes)
   call print_param("(y) min wavelength", 2.0_rp*pi/kmax, fmt='F10.4')
   call print_param("(y) max wavelength", 2.0_rp*pi/kmin, fmt='F10.4')

   ! z-direction
   kmin = glmin(abs(this%k_z), this%n_modes)
   kmax = glmax(abs(this%k_z), this%n_modes)
   call print_param("(z) min wavelength", 2.0_rp*pi/kmax, fmt='F10.4')
   call print_param("(z) max wavelength", 2.0_rp*pi/kmin, fmt='F10.4')
   call neko_log%end_section()

   !
   ! From the random vectors and amplitudes, compute estimations of Tu and 
   ! TKE
   !
   call neko_log%section("Amplitudes & rand. vectors")

   block
    integer :: shellno, i
    real(kind=rp) :: amp, uamp, vamp, wamp, ue, ve, we

    do i=1, this%n_modes
      shellno = this%shell(i)
      amp = this%shell_amp(shellno)

      uamp = this%random_vectors(i,1)*amp
      vamp = this%random_vectors(i,2)*amp
      wamp = this%random_vectors(i,3)*amp

      ue = ue + ((uamp)**2.0_rp)/2.0_rp
      ve = ve + ((vamp)**2.0_rp)/2.0_rp
      we = we + ((wamp)**2.0_rp)/2.0_rp
    enddo

    call print_param('Energy in u  ', ue, fmt='E12.6')
    call print_param('Energy in v  ', ve, fmt='E12.6')
    call print_param('Energy in w  ', we, fmt='E12.6')
    call print_param('Target TKE   ', 1.5_rp * this%Uinf**2 * this%Tu**2, fmt='E12.6')
    call print_param('Estimated TKE', (ue+ve+we)/2.0_rp, fmt='E12.6')
    call print_param('Target Tu    ', this%Tu, fmt='E12.6')
    call print_param('Estimated Tu ', sqrt((ue+ve+we)/3.0_rp / this%Uinf), fmt='E12.6')

   end block

   call neko_log%end_section() ! end amplitudes section


   call neko_log%end_section() ! end diagnostics section

  end subroutine FST_validate

  !> Generate FST for forcing
  subroutine FST_generate_forcing(this, coef, zone, u, v, w, path)
    class(FST_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    class(point_zone_t), intent(in) :: zone
    type(field_t), intent(in) :: u, v, w
    character(len=*), intent(in) :: path

    real(kind=rp) :: x, y, z
    integer :: ierr, i, idx

    integer, pointer :: mask(:)

    ! Do the general generation
    call this%generate_common(path, coef%msh%gdim)

    !
    ! Copy the baseflow in the zone
    !
    call this%apply_baseflow(mask, zone%size, u, v, w)

    ! Generate the fringe in space
    allocate(this%fringe_space(zone%size))

    ! Initialize the fringe in space
    do idx = 1, zone%size
       i = mask(idx)
       x = coef%dof%x(i,1,1,1)
       y = coef%dof%y(i,1,1,1)
       z = coef%dof%z(i,1,1,1)

       if ( x .lt. this%xmin .or. x .gt. this%xmax .or. &
            y .lt. this%ymin .or. y .gt. this%ymax ) then
          this%fringe_space(idx) = 0.0_rp
       else
          this%fringe_space(idx) = fringe(x, y, this)
       end if

    end do

  end subroutine FST_generate_forcing

  !> Do the generation for BC.
  subroutine FST_generate_bc(this, x_dof, y_dof, z_dof, bc_mask, n, u, v, w, &
      path, gdim)
    class(FST_t), intent(inout) :: this
    real(kind=rp), intent(in) :: x_dof(:,:,:,:), y_dof(:,:,:,:), z_dof(:,:,:,:)
    integer, intent(in) :: bc_mask(0:n)
    integer, intent(in) :: n
    type(field_t), intent(in) :: u, v, w
    character(len=*), intent(in) :: path
    integer, intent(in) :: gdim

    real(kind=rp) :: x, y, z, ymin, ymax, zmin, zmax, Ly, Lz
    integer :: ierr, i, idx, m, j

    !
    ! Compute the spatial bounds in the mask
    !
    ymin = huge(1.0_rp); zmin = huge(1.0_rp)
    ymax = -huge(1.0_rp); zmax = -huge(1.0_rp)
    do idx = 1, n
      i = bc_mask(idx)
      ymin = min(ymin, y_dof(i,1,1,1))
      zmin = min(zmin, z_dof(i,1,1,1))
      ymax = max(ymax, y_dof(i,1,1,1))
      zmax = max(zmax, z_dof(i,1,1,1))
    end do

    call mpi_allreduce(MPI_IN_PLACE, ymin, 1, &
         mpi_real_precision, mpi_min, neko_comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, zmin, 1, &
         mpi_real_precision, mpi_min, neko_comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, ymax, 1, &
         mpi_real_precision, mpi_max, neko_comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, zmax, 1, &
         mpi_real_precision, mpi_max, neko_comm, ierr)

    Ly = ymax - ymin
    Lz = zmax - zmin

    !
    ! Update the fringe parameters with new bounds for periodic directions
    ! NOTE the notations here are a bit confusing with (x,y) but our directions
    ! are (y,z). 
    ! The point of updating the parameters in periodic directions is to
    ! have a fringe of constant value 1.0 across the entire boundary, therefore
    ! to stretch the bounds xmin/max etc so the fringe region is far outside the
    ! domain
    !
    if (this%periodic_z) then
      this%xmin = zmin - 10.0_rp * Lz 
      this%xmax = zmax + 10.0_rp * Lz
      this%xstart = this%xmin
      this%xend = this%xmax
      this%x_delta_rise = 0.001_rp * Lz
      this%x_delta_fall = 0.001_rp * Lz
    end if

    if (this%periodic_y) then
      this%ymin = ymin - 10.0_rp * Ly 
      this%ymax = ymax + 10.0_rp * Ly
      this%ystart = this%ymin
      this%yend = this%ymax
      this%y_delta_rise = 0.1_rp * Ly
      this%y_delta_fall = 0.1_rp * Ly
    end if

    !
    ! Do the common generation (not passing Lx since it is not supported yet)
    !
    call this%generate_common(path, gdim, Ly = Ly, Lz = Lz)

    !
    ! Apply baseflow in the bc zone
    !
    call this%apply_baseflow_0(bc_mask, n, u, v, w)

    !
    ! Initialize the fringe in space
    !
    allocate(this%fringe_space(n))
    do idx = 1, size(this%fringe_space)

       i = bc_mask(idx)
       y = y_dof(i,1,1,1)
       z = z_dof(i,1,1,1)

       if ( z .lt. this%xmin .or. z .gt. this%xmax .or. &
            y .lt. this%ymin .or. y .gt. this%ymax ) then
          this%fringe_space(idx) = 0.0_rp
       else
          this%fringe_space(idx) = fringe(z, y, this)
       end if

    end do

    !
    ! Precompute spatial, time-independent term
    !
    allocate(this%phi_0(this%n_modes, n))

    do j = 1, n
       x = x_dof(bc_mask(j), 1,1,1)
       y = y_dof(bc_mask(j), 1,1,1)
       z = z_dof(bc_mask(j), 1,1,1)
       do m = 1, this%n_modes
          this%phi_0(m,j) = this%k_x(m)*x + this%k_y(m)*y + this%k_z(m)*z &
               + this%phase_shifts(m)
       end do
    end do

    ! Copy everything to device and map with relevant device array pointers
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%fringe_space, this%fringe_space_d, n)
       call device_memcpy(this%fringe_space, this%fringe_space_d, n, &
            HOST_TO_DEVICE, .false.)

       call device_map(this%phi_0, this%phi_0_d, this%n_modes*n)
       call device_memcpy(this%phi_0, this%phi_0_d, this%n_modes*n, &
            HOST_TO_DEVICE, .false.)

       call device_map(this%random_vectors, this%random_vectors_d, this%n_modes*3)
       call device_memcpy(this%random_vectors, this%random_vectors_d, this%n_modes*3, &
            HOST_TO_DEVICE, .false.)

       call device_map(this%shell, this%shell_d, this%n_modes)
       call device_memcpy(this%shell, this%shell_d, this%n_modes, &
            HOST_TO_DEVICE, .false.)

       call device_map(this%shell_amp, this%shell_amp_d, this%n_shells)
       call device_memcpy(this%shell_amp, this%shell_amp_d, this%n_shells, &
            HOST_TO_DEVICE, .false.)

       call device_map(this%k_x, this%k_x_d, this%n_modes)
       call device_memcpy(this%k_x, this%k_x_d, this%n_modes, &
            HOST_TO_DEVICE, .true.)
    end if

  end subroutine FST_generate_bc

  ! Forcing to be performed on entire domain, on a local element ix
  ! Final values of the forcing are to be applied
  ! on f%u, f%v and f%w
!   subroutine FST_forcing_zone(this, zone, &
!        x, y, z, t, &
!        u, v, w, &
!        fu, fv, fw)

!     class(FST_t), intent(in) :: this
!     class(point_zone_t), intent(in) :: zone
!     real(kind=rp), intent(in), dimension(:,:,:,:) :: x, y, z, u, v, w
!     real(kind=rp), intent(in) :: t
!     real(kind=rp), intent(inout), dimension(:,:,:,:) :: fu, fv, fw

!     integer :: idx, l, m, i, shellno
!     integer, parameter :: gdim = 3
!     real(kind=rp) :: phase_shft, phi, amp, pert, vel_mag
!     real(kind=rp) :: rand_vec(gdim)
!     real(kind=rp) :: fringe_time

!     integer, pointer :: mask(:)
!     mask => zone%mask%get()

!     fringe_time = time_ramp(t, this%t_end, this%t_start)

!     ! Loop on all points in the point zone
!     do idx = 1, zone%size
!        i = mask(idx)

!        !> This vector will contain the sum of all Fourier modes
!        rand_vec = 0.0_rp

!        ! Sum all sin modes for each gll point
!        do m = 1, this%n_modes
!           phase_shft = bb(m,1)

!           ! kx(x - U*t) + ky*y + kz*z - random_phase[-pi, pi]
!           ! = kx*x + ky*y + kz*z - kx*U*t - random_phase
!           phi = k_num_all(m,1) * (x(i,1,1,1) - glb_uinf*t) + &
!                k_num_all(m,2) * y(i,1,1,1) + &
!                k_num_all(m,3) * z(i,1,1,1) + &
!                phase_shft

!           shellno = shell(m)
!           pert = shell_amp(shellno)*sin(phi)

!           rand_vec(1) = rand_vec(1) + u_hat_pn(m,1)*pert
!           rand_vec(2) = rand_vec(2) + u_hat_pn(m,2)*pert
!           rand_vec(3) = rand_vec(3) + u_hat_pn(m,3)*pert
!        enddo

!        fu(i,1,1,1) = fringe_time*this%fringe_space(idx)*( &
!             this%u_baseflow(idx) + rand_vec(1) - u(i,1,1,1))

!        fv(i,1,1,1) = fringe_time*this%fringe_space(idx)*( &
!             this%v_baseflow(idx) + rand_vec(2) - v(i,1,1,1))

!        fw(i,1,1,1) = fringe_time*this%fringe_space(idx)*( &
!             this%w_baseflow(idx) + rand_vec(3) - w(i,1,1,1))
!     end do

!   end subroutine FST_forcing_zone

  ! Apply FST as a boundary condition based on the bc mask
  ! Assumes that u,v,w have the same bc masks
  subroutine FST_apply_BC(this, bc_mask, n, &
       x, y, z, t, &
       u_bc, v_bc, w_bc, angleXY, on_host)

    class(FST_t), intent(in) :: this
    integer, intent(in) :: n ! size of the bc mask
    integer, intent(in) :: bc_mask(0:n)
    real(kind=rp), intent(in), dimension(:,:,:,:) :: x, y, z
    real(kind=rp), intent(in) :: t
    real(kind=rp), intent(inout), dimension(:,:,:,:) :: u_bc, v_bc, w_bc
    real(kind=rp), intent(in) :: angleXY
    logical, intent(in) :: on_host

    real(kind=rp) :: fringe_time

    fringe_time = time_ramp(t, this%t_end, this%t_start)

    call fst_bc_compute(t, this%Uinf, u_bc, v_bc, w_bc, bc_mask, n, &
         this%u_baseflow, this%v_baseflow, this%w_baseflow, &
         this%k_x, this%n_modes, this%phi_0, this%shell, this%shell_amp, &
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
!!$       do m = 1, this%n_modes
!!$
!!$          ! Random phase shifts
!!$          !phase_shft = bb(m,1)
!!$
!!$          ! kx(x - U*t) + ky*y + kz*z + phi
!!$          ! = kx*x + ky*y + kz*z - kx*U*t + phi
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
  !      S=1/(1+erp(1/(x-f%xend)+1/(x-f%xstart)))
  !   endif

  !   y = f%fringe_max * (S*(x-f%xstart)/(f%delta_rise)-S*((x-f%xend)/f%delta_fall+1))

  ! end function fringe

  !
  ! Linear ramp in time
  function time_ramp(t, t_end, t_start) result(ramp)
    real(kind=rp), intent(in) :: t
    real(kind=rp), intent(in) :: t_end
    real(kind=rp), intent(in) :: t_start

    real(kind=rp) :: ramp

    if (t .le. t_start) then
       ramp = 0.0_rp
    else if (t .lt. t_end) then
       ramp = (t - t_start)/(t_end - t_start)
    else
       ramp = 1.0_rp
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
    real(kind=rp) :: fr
    integer :: i
    character :: a

    fr = f%fringe_max * ( S((x-f%xstart)/f%x_delta_rise) - S((x-f%xend)/f%x_delta_fall + 1.0_rp) )

    if (present(y)) then
       fr = fr * ( S((y-f%ystart)/f%y_delta_rise) - S((y-f%yend)/f%y_delta_fall + 1._rp) )
    end if

  end function fringe

  ! Smooth step function, 0 if x <= 0, 1 if x >= 1, 1/erp(1/(x-1) + 1/x) between 0 and 1
  function S(x) result(y)
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: y

    if ( x.le.0._rp ) then
       y = 0._rp
    else if ( x.ge.1._rp ) then
       y = 1._rp
    else
       y = 1._rp / (1._rp + exp( 1._rp/(x-1._rp) + 1._rp/x))
    end if

  end function S


end module FST
