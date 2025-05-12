! * Free stream turbulence implementation for neko
! * Author: Victor Baconnet (baconnet@kth.se)
! * Based in original implementation by Elektra Kluesberg, Prabal Negi
! 	and Philipp Schlatter. 
! * For FST theory see Schlatter (2001)

module FST

  use global_params
  use turbu
  use utils, only: neko_error
  use point_zone, only: point_zone_t
  use device, only: HOST_TO_DEVICE, device_free, device_memcpy, &
          device_map
  use neko_config, only: NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding
  implicit none


  type, public :: FST_t

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
     real(kind=rp) :: fringe_max

     !> Final ramp time
     real(kind=rp) :: t_end
     real(kind=rp) :: t_start

     logical :: is_forcing
     logical :: is_bc

     !> Fringe in space
     type(vector_t) :: fringe_space

     !> Baseflows, if applying on a non-uniform inflow
     type(vector_t) :: u_baseflow
     type(vector_t) :: v_baseflow
     type(vector_t) :: w_baseflow

   contains

     ! ======== Init/Free procedures
     procedure, pass(this) :: init_common => FST_init_common
     procedure, pass(this) :: init_bc => FST_init_bc
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
    real(kind=rp), intent(in) :: fringe_max
    real(kind=rp), intent(in) :: t_start
    real(kind=rp), intent(in) :: t_end
    logical, intent(in) :: periodic_x, periodic_y, periodic_z

    this%periodic_x = periodic_x
    this%periodic_y = periodic_y
    this%periodic_z = periodic_z

    this%xstart = xstart
    this%xend   = xend
    this%ystart = ystart
    this%yend   = yend

    this%xmin   = this%xstart + x_delta_rise/3.0_rp
    this%xmax   = this%xend   - x_delta_fall/3.0_rp
    this%ymin   = this%ystart + y_delta_rise/3.0_rp
    this%ymax   = this%yend   - y_delta_fall/3.0_rp

    this%fringe_max = fringe_max
    this%x_delta_rise = x_delta_rise!0.002
    this%x_delta_fall = x_delta_fall!0.002
    this%y_delta_rise = y_delta_rise!0.002
    this%y_delta_fall = y_delta_fall!0.002
    this%t_end = t_end
    this%t_start = t_start

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
    real(kind=rp), intent(in) :: fringe_max
    real(kind=rp), intent(in) :: t_start
    real(kind=rp), intent(in) :: t_end
    logical, intent(in) :: periodic_x, periodic_y, periodic_z

    call neko_log%section('Initializing FST')

    call this%init_common(xstart, xend, x_delta_rise, x_delta_fall, ystart, &
         yend, y_delta_rise, y_delta_fall, fringe_max, t_start, t_end, &
         periodic_x, periodic_y, periodic_z)

    call this%print() ! show parameters
    call neko_log%end_section('Done --> Intializing FST')

    this%is_forcing = .true.
    this%is_bc = .false.

  end subroutine FST_init_forcing

  !> Initialize the FST to use as a boundary condition
  subroutine FST_init_bc(this, &
       xstart, xend, x_delta_rise, x_delta_fall, &
       ystart, yend, y_delta_rise, y_delta_fall, &
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
    real(kind=rp), intent(in) :: t_start
    real(kind=rp), intent(in) :: t_end
    logical, intent(in) :: periodic_x, periodic_y, periodic_z

    call neko_log%section('Initializing FST')

    call this%init_common(xstart, xend, x_delta_rise, x_delta_fall, ystart, &
         yend, y_delta_rise, y_delta_fall, 1.0_rp, t_start, t_end, &
         periodic_x, periodic_y, periodic_z)

    call this%print() ! show parameters
    call neko_log%end_section('Done --> Intializing FST')

    this%is_forcing = .false.
    this%is_bc = .true.

  end subroutine FST_init_bc

  !! Free parameters in global params
  subroutine FST_free_params(this)
    class(FST_t), intent(inout) :: this

    call this%u_baseflow%free()
    call this%v_baseflow%free()
    call this%w_baseflow%free()
    call this%fringe_space%free()

  end subroutine FST_free_params

  subroutine FST_print_params(this)
    class(FST_t) :: this

    call print_param("nshells", real(nshells, kind=rp))
    call print_param("Npmax", real(Npmax, kind=rp))
    call print_param("fst_ti", fst_ti)
    call print_param("fst_il", fst_il)
    call print_param("kstart", kstart)
    call print_param("kend", kend)
    call print_param("xstart", this%xstart)
    call print_param("xend", this%xend)
    call print_param("ystart", this%ystart)
    call print_param("yend", this%yend)
    call print_param("fringe_max", this%fringe_max)
    call print_param("x_delta_rise", this%x_delta_rise)
    call print_param("x_delta_fall", this%x_delta_fall)
    call print_param("y_delta_rise", this%y_delta_rise)
    call print_param("y_delta_fall", this%y_delta_fall)
    call print_param("t_start", this%t_start)
    call print_param("t_end", this%t_end)

  end subroutine FST_print_params

  !> Apply a specific baseflow in our region, from a boundary mask.
  !! If we specify v or w it takes the norm.
  subroutine FST_apply_baseflow(this, mask, n, u, v, w)
    class(FST_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w
    integer, intent(in) :: mask(n)
    integer, intent(in) :: n

    integer :: i, idx

    call this%u_baseflow%init(n)
    call this%v_baseflow%init(n)
    call this%w_baseflow%init(n)

    ! Could use a gather/scatter copy here?
    do i = 1, n
       idx = mask(i)
       this%u_baseflow%x(i) = u%x(idx,1,1,1)
       this%v_baseflow%x(i) = v%x(idx,1,1,1)
       this%w_baseflow%x(i) = w%x(idx,1,1,1)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%u_baseflow%x, this%u_baseflow%x_d, &
               this%u_baseflow%size(), HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(this%v_baseflow%x, this%v_baseflow%x_d, &
               this%v_baseflow%size(), HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(this%w_baseflow%x, this%w_baseflow%x_d, &
               this%w_baseflow%size(), HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine FST_apply_baseflow

  !> Common function for generation
  subroutine FST_generate_common(this, coef)
    class(FST_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    integer :: ierr

    call neko_log%section ('Generating FST')
    call make_turbu(coef, this%periodic_x, this%periodic_y, this%periodic_z)

    call MPI_Bcast(k_length , 1                   , &
         MPI_INTEGER         , 0, NEKO_COMM, ierr)
    call MPI_Bcast(k_num_all, fst_modes*coef%msh%gdim, &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(k_num    , fst_modes*coef%msh%gdim, &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(u_hat_pn , fst_modes*coef%msh%gdim, &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(bb       , fst_modes*coef%msh%gdim, &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(shell    , fst_modes           , &
         MPI_INTEGER         , 0, NEKO_COMM, ierr)
    call MPI_Bcast(shell_amp, nshells             , &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)

    call neko_log%end_section('Done --> Generating FST')

  end subroutine FST_generate_common

  !> Generate FST for forcing
  subroutine FST_generate_forcing(this, coef, zone, u, v, w)
    class(FST_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    class(point_zone_t), intent(in) :: zone
    type(field_t), intent(in) :: u, v, w

    real(kind=rp) :: x, y, z
    integer :: ierr, i, idx

    ! Do the general generation
    call this%generate_common(coef)

    !
    ! Copy the baseflow in the zone
    !
    call this%apply_baseflow(zone%mask, zone%size, u, v, w)

    ! Generate the fringe in space
    allocate(this%fringe_space(zone%size))

    ! Initialize the fringe in space
    do idx = 1, zone%size
       i = zone%mask(idx)
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
  subroutine FST_generate_bc(this, coef, bc_mask, n, u, v, w)
    class(FST_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: bc_mask(0:n)
    integer, intent(in) :: n
    type(field_t), intent(in) :: u, v, w

    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: x, y, z
    integer :: ierr, i, idx

    ! Do the general generation
    call this%generate_common(coef)

    !
    ! Apply baseflow in the bc zone
    !
    call this%apply_baseflow(bc_mask(1:n), bc_mask(0), u, v, w)

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
          this%fringe_space(idx) = 0.0_rp
       else
          this%fringe_space(idx) = fringe(z, y, this)
       end if

    end do

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

    fringe_time = time_ramp(t, this%t_end, this%t_start)

    ! Loop on all points in the point zone 
    do idx = 1, zone%size
       i = zone%mask(idx)

       !> This vector will contain the sum of all Fourier modes
       rand_vec = 0.0_rp

       ! Sum all sin modes for each gll point
       do m = 1, k_length
          phase_shft = bb(m,1)

          ! kx(x - U*t) + ky*y + kz*z - random_phase[-pi, pi]
          ! = kx*x + ky*y + kz*z - kx*U*t - random_phase
          phi = k_num_all(m,1) * (x(i,1,1,1) - glb_uinf*t) + &
               k_num_all(m,2) *  y(i,1,1,1) + &
               k_num_all(m,3) *  z(i,1,1,1) + &
               phase_shft

          shellno = shell(m)
          pert = shell_amp(shellno)*sin(phi)

          rand_vec(1) = rand_vec(1) + u_hat_pn(m,1)*pert
          rand_vec(2) = rand_vec(2) + u_hat_pn(m,2)*pert
          rand_vec(3) = rand_vec(3) + u_hat_pn(m,3)*pert
       enddo

       fu(i,1,1,1) = fringe_time*this%fringe_space(idx)*( &
            this%u_baseflow(idx) + rand_vec(1) - u(i,1,1,1))

       fv(i,1,1,1) = fringe_time*this%fringe_space(idx)*( &
            this%v_baseflow(idx) + rand_vec(2) - v(i,1,1,1))

       fw(i,1,1,1) = fringe_time*this%fringe_space(idx)*( &
            this%w_baseflow(idx) + rand_vec(3) - w(i,1,1,1))
    end do

  end subroutine FST_forcing_zone

  ! Apply FST as a boundary condition based on the bc mask
  ! Assumes that u,v,w have the same bc masks
  subroutine FST_apply_BC(this, bc_mask, n, &
       x, y, z, t, &
       u_bc, v_bc, w_bc, angleXY)

    class(FST_t), intent(in) :: this
    integer, intent(in) :: n ! size of the bc mask
    integer, intent(in) :: bc_mask(0:n)
    real(kind=rp), intent(in), dimension(:,:,:,:) :: x, y, z
    real(kind=rp), intent(in) :: t
    real(kind=rp), intent(inout), dimension(:,:,:,:) :: u_bc, v_bc, w_bc
    real(kind=rp), intent(in) :: angleXY

    integer :: idx, l, m, i, shellno
    integer, parameter :: gdim = 3
    real(kind=rp) :: phase_shft, phi, amp, pert, urand, vrand, wrand
    real(kind=rp) :: rand_vec(gdim), fringe_time, vel_mag

    fringe_time = time_ramp(t, this%t_end, this%t_start)

    ! Loop on all points in the point zone
    do idx = 1, bc_mask(0)

       i = bc_mask(idx)

       !> This vector will contain the sum of all Fourier modes
       rand_vec = 0.0_rp
       
       ! Sum all sin modes for each gll point
       do m = 1, k_length

          ! Random phase shifts
          phase_shft = bb(m,1)

          ! kx(x - U*t) + ky*y + kz*z + phi
          ! = kx*x + ky*y + kz*z - kx*U*t + phi
          phi = k_num_all(m,1) * (x(i,1,1,1) - glb_uinf*t) + &
               k_num_all(m,2) *  y(i,1,1,1) + &
               k_num_all(m,3) *  z(i,1,1,1) + &
               phase_shft

          shellno = shell(m)

          pert = shell_amp(shellno)*sin(phi)

          rand_vec(1) = rand_vec(1) + u_hat_pn(m,1)*pert
          rand_vec(2) = rand_vec(2) + u_hat_pn(m,2)*pert
          rand_vec(3) = rand_vec(3) + u_hat_pn(m,3)*pert

       enddo

       ! Project the rand_vec if there is a rotation
       urand = rand_vec(1)*cos(angleXY) - rand_vec(2)*sin(angleXY)
       vrand = rand_vec(1)*sin(angleXY) + rand_vec(2)*cos(angleXY)
       wrand = rand_vec(3)

       u_bc(i,1,1,1) = this%u_baseflow(idx) + &
            fringe_time*this%fringe_space(idx)*urand
       v_bc(i,1,1,1) = this%v_baseflow(idx) + &
            fringe_time*this%fringe_space(idx)*vrand
       w_bc(i,1,1,1) = this%w_baseflow(idx) + &
            fringe_time*this%fringe_space(idx)*wrand

    end do

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
  !   xstart   xmin         xmax   xend
  !
  ! The distance between xstart and xmin is delta_rise/3, and the ramp-up
  ! after xmin is of length delta_rise. The same applies for the ramp down
  ! and the distance between xmax and xend. If y is specified then it computes
  ! a 2D fringe
  function fringe(x, y, f) result(fr)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in), optional :: y
    type(FST_t), intent(in) :: f
    real(kind=rp) :: fr
    integer :: i
    character :: a

    fr = f%fringe_max * ( S((x-f%xstart)/f%x_delta_rise) - S((x-f%xend)/f%x_delta_fall + 1d0) )

    if (present(y)) then
       fr = fr * ( S((y-f%ystart)/f%y_delta_rise) - S((y-f%yend)/f%y_delta_fall + 1d0) )
    end if

  end function fringe

  ! Smooth step function, 0 if x <= 0, 1 if x >= 1, 1/exp(1/(x-1) + 1/x) between 0 and 1
  function S(x) result(y)
    real(kind=rp), intent(in) :: x
    real(kind=rp)             :: y

    if ( x.le.0d0 ) then
       y = 0d0
    else if ( x.ge.1d0 ) then
       y = 1d0
    else
       y = 1d0 / (1d0 + exp( 1d0/(x-1d0) + 1d0/x))
    end if

  end function S


end module FST
