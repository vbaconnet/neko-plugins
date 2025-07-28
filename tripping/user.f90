module user
  use neko
  use sponge
  use fst_bc_driver
  use trip
  implicit none

  type(matrix_t) :: INFLOW_PROFILE
  type(trip_t) :: tripper
  logical :: enable_tripping
  !type(sponge_t) :: spng

  !Things for tripping
  ! Parameters taken from Schlatter & Örlu 2012 doi:10.1017/jfm.2012.324
  ! NOTE: Tripping only works for 1 line parallel to the y-direction
  real(kind=rp), parameter :: delta_star = 0.07_rp ! displacement thickness of our inflow BL

  real(kind=rp), parameter :: TIAMP   = 0.5_rp       ! Time independent amplitude
  real(kind=rp), parameter :: TDAMP   = 1.0_rp       ! Time dependent amplitude

  real(kind=rp), parameter :: SPOSX01 = -0.80_rp     ! Starting point (x)
  real(kind=rp), parameter :: SPOSY01 = 0.0_rp         ! Starting point (y)
  real(kind=rp), parameter :: SPOSZ01 = 0.002_rp         ! Starting point (z)

  real(kind=rp), parameter :: EPOSX01 = SPOSX01    ! Ending point (x)
  real(kind=rp), parameter :: EPOSY01 = 1.0_rp       ! Ending point (y)
  real(kind=rp), parameter :: EPOSZ01 = SPOSZ01      ! Ending point (z)

  real(kind=rp), parameter :: SMTHX01 = 3.0_rp * delta_star !0.15_rp       ! Smoothing length (x)
  real(kind=rp), parameter :: SMTHY01 = 0.05_rp      ! Smoothing length (y)
  real(kind=rp), parameter :: SMTHZ01 = delta_star !0.002_rp       ! Smoothing length (z)

  real(kind=rp), parameter :: ROTA01  = 0.0*3.1415926_rp/180.0_rp       ! Rotation angle (radians)
  integer, parameter :: NMODE01 = 20             ! Number of Fourier modes
  real(kind=rp), parameter :: TDT01 = 4.0_rp * delta_star / 1.0_rp !0.001_rp    ! Time step for tripping
  integer, parameter :: nline = 1                 ! Number of tripping lines
  logical, parameter :: LEXT01 = .true.          ! Line extension (??)

  integer :: counter = 1

  integer, parameter :: BL_REFERENCE = 2, BL_THINNER_OUTER = 3, &
       BL_THINNER_INNER = 4

contains


  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_dirichlet_update => set_boundary_conditions
    u%fluid_user_f_vector => userf
    u%fluid_user_f => tripline_org
    u%user_init_modules => initialize
    u%user_finalize_modules => finalize
    u%user_check => usercheck
    u%fluid_user_ic => user_ic
  end subroutine user_setup

  ! Finalize user variables or external objects
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    call INFLOW_PROFILE%free() 
    call fst_bc_driver_finalize(t, params)
    !call spng%free()

  end subroutine finalize

  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j, ierr
    real(kind=rp) :: z, vel
 
    type(file_t) :: f
    !call init_file() 

    !
    ! Read the inflow profile
    !
    f = file_t("blasius_custom_half.csv")

    ! There are 22 lines in the input file with the header
    call INFLOW_PROFILE%init(203, 2)
    call f%read(INFLOW_PROFILE)

    select type(ft => f%file_type)
    class is (csv_file_t)
    call neko_log%message("HEADER " // ft%header)
    class default
    end select

    call MPI_Bcast(INFLOW_PROFILE%x , INFLOW_PROFILE%n , MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)
    call file_free(f)

    w = 0.0_rp

    do i = 1, u%dof%size()
       z = u%dof%z(i,1,1,1)
   
       vel = local_interpolation(z, INFLOW_PROFILE%x(:,2),&
            INFLOW_PROFILE%x(:,1), INFLOW_PROFILE%nrows)

       ! Project the BL profile u = cos(alpha), v = sin(alpha)
       u%x(i,1,1,1) = vel*cos(40.0_rp*PI/180.0_rp)
       v%x(i,1,1,1) = vel*sin(40.0_rp*PI/180.0_rp)
    end do
  
    call neko_log%message("useric")

     if (neko_bcknd_device .eq. 1) then
       call device_memcpy(u%x, u%x_d, u%dof%size(), &
                          host_to_device, sync=.false.)
       call device_memcpy(v%x, v%x_d, v%dof%size(), &
                          host_to_device, sync=.false.)
       call device_memcpy(w%x, w%x_d, w%dof%size(), &
                          host_to_device, sync=.false.)
    end if
   
  end subroutine user_ic

  subroutine initialize(t, u, v, w, p, c_Xh, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: c_Xh
    type(json_file), intent(inout) :: params
    integer :: nmode(trip_nline_max), ierr
    real(kind=rp), dimension(3,trip_nline_max) :: spos, epos, smth
    real(kind=rp) :: rota(trip_nline_max), tdt(trip_nline_max)
    logical :: lext(trip_nline_max)


    type(field_t), pointer :: fx, fy, fz

    !
    ! Tripping
    !

    call json_get(params, "case.tripping.enabled", enable_tripping)

    if (.not. enable_tripping) call neko_log%message("<<< Tripping disabled >>>")
    if (enable_tripping) call neko_log%message("<<< Tripping enabled >>>")

    spos(1,1) = SPOSX01
    spos(2,1) = SPOSY01
    spos(3,1) = SPOSZ01
    epos(1,1)= EPOSX01
    epos(2,1)= EPOSY01
    epos(3,1)= EPOSZ01
    smth(1,1)= SMTHX01
    smth(2,1)= SMTHY01
    smth(3,1)= SMTHZ01
    rota(1) = ROTA01
    nmode(1) = NMODE01
    tdt(1) = TDT01
    lext(1) = LEXT01

    call tripper%init( u%dof, nline, nmode, tiamp, tdamp, &
                       spos, epos, smth, lext, rota, tdt, t)

    !
    ! Sponge
    !
    !call spng%init_from_attributes(0.1_rp, 100.0_rp)
    !call spng%assign_baseflow_field(u,v,w)
    !call spng%compute_fringe((/"outflow"/),1)

    !
    ! FST
    !
    call fst_bc_driver_initialize(t, u, v, w, p, c_Xh, params)
  
  end subroutine initialize

  subroutine init_file
    !character(len=*), intent(inout) :: fname
    !integer, intent(in) :: nrows, ncols

    type(file_t) :: f
    integer :: ierr

    if (allocated(INFLOW_PROFILE%x)) return

    f = file_t(".csv")

    ! There are 22 lines in the input file with the header
    call INFLOW_PROFILE%init(21, 4)
    call f%read(INFLOW_PROFILE)

    select type(ft => f%file_type)
    class is (csv_file_t)
    call neko_log%message("HEADER " // ft%header)
    class default
    end select

    call MPI_Bcast(INFLOW_PROFILE%x , INFLOW_PROFILE%n , MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)
    call file_free(f)

  end subroutine init_file

  ! user-defined boundary condition
  subroutine set_boundary_conditions(field_bc_list, bc, coef, t, tstep)
    type(field_list_t), intent(inout) :: field_bc_list
    type(field_dirichlet_t), intent(in) :: bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    integer :: i, idx_u, idx_v
    real(kind=rp) :: z, vel, ts, te

    ! Only do this at the first time step since our BCs are constants.
    !if (tstep .ne. 1) return

    ! Check that we are being called by `fluid`
    if (field_bc_list%items(1)%ptr%name .eq. "u") then

       associate(u => field_bc_list%items(1)%ptr, &
            v => field_bc_list%items(2)%ptr, &
            w => field_bc_list%items(3)%ptr)

         if (tstep .eq. 1) then
            !
            ! Perform operations on u%x, v%x, w%x and p%x here
            ! Note that we are checking if fields are allocated. If the
            ! boundary types only contains e.g. "d_vel_u/d_pres", the fields
            ! v%x and w%x will not be allocated.
            !

            ! The masks should be the same, but we check just in case

            do i = 1, bc%msk(0)

               idx_u = bc%msk(i)

               z = coef%dof%z(idx_u, 1, 1, 1)

               ! Interpolated BL (multiply y by 0.001 to get in mm)
               vel = local_interpolation(z, INFLOW_PROFILE%x(:,2),&
                    INFLOW_PROFILE%x(:,1), INFLOW_PROFILE%nrows)

               ! Project the BL profile u = cos(alpha), v = sin(alpha)
               u%x(idx_u,1,1,1) = vel*cos(40.0_rp*PI/180.0_rp)
               v%x(idx_u,1,1,1) = vel*sin(40.0_rp*PI/180.0_rp)
            end do

            w = 0.0_rp

         end if ! end if tstep .eq. 1

         ts = MPI_Wtime()
         call fst_bc_driver_apply(u, v, w, bc, coef, t, tstep, 40.0_rp * PI/180_rp)
         te = MPI_Wtime()
         if (pe_rank .eq. 0) write (*,*) ">>> FST apply: ", te - ts

         ! No need to copy v and w since the v = 0 does it on the GPU
         if (neko_bcknd_device .eq. 1) then
            ts = MPI_Wtime()
            call device_memcpy(u%x, u%x_d, u%dof%size(), &
                 host_to_device, sync=.false.)
            call device_memcpy(v%x, v%x_d, v%dof%size(), &
                 host_to_device, sync=.false.)
            call device_memcpy(w%x, w%x_d, w%dof%size(), &
                 host_to_device, sync=.false.)
            te = MPI_Wtime()
            if (pe_rank .eq. 0) write (*,*) ">>> FST copy: ", te - ts
         end if

       end associate

    end if

  end subroutine set_boundary_conditions

  !> Tripline forcing
  subroutine tripline_org(u, v, w, j, k, l, e , t)
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    real(kind=rp), intent(in) :: t

    u = 0.0
    v = 0.0
    w = 0.0

!!    write(*,*) 'inside tripline_org'
    call tripper%apply(u,v,w,j,k,l,e)

  end subroutine tripline_org

  subroutine userf(f,t)
!!!!    class(source_t) :: f
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i, j, k, e
    real(kind=rp) :: u, v, w

    type(field_t), pointer :: fx, fy, fz

    ! Compute sponge
    !call spng%compute(f, t, 0)

    ! Compute tripping
    !fx => neko_field_registry%get_field('fx')
    !fy => neko_field_registry%get_field('fy')
    !fz => neko_field_registry%get_field('fz')

    if (.not. enable_tripping) return

    if (t .le. 3.0_rp) return

    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_rzero(f%u_d,tripper%dof%size())
       call device_rzero(f%v_d,tripper%dof%size())
       call device_rzero(f%w_d,tripper%dof%size())
       call tripper%apply_device(f%u_d, f%v_d, f%w_d)
    else
       do e = 1, tripper%msh%nelv
          do k = 1, tripper%Xh%lz
             do j = 1, tripper%Xh%ly
                do i = 1, tripper%Xh%lx
                   call tripline_org(u,v,w,i,j,k,e,t)
                   f%u(i,j,k,e) = u
                   f%v(i,j,k,e) = v
                   f%w(i,j,k,e) = w
                end do
             end do
          end do
       end do
    end if
    
    !call copy(fx%x, f%u, fx%dof%size())
    !call copy(fy%x, f%v, fy%dof%size())
    !call copy(fz%x, f%w, fz%dof%size())

  end subroutine userf

  !> Linearly interpolates a point from an array of x and u(x).
  !! The function finds the interval [x1, x2] that contains x_GLL,
  !! and interpolates linearly between [u1, u2]:
  !!
  !! u_GLL = (u2 - u1)/(x2 - x1)*(x_GLL - x1) + u1
  !!
  !! @param x_GLL Coordinate of the point at which to interpolate
  !! @param u Array of field from which to interpolate (must be in the same order as x)
  !! @param x Array of coordinates from which to interpolate
  !! @param n Size of u and x
  function local_interpolation(x_GLL, u, x, n) result(u_GLL)
    real(kind=rp), intent(in) :: x_GLL
    real(kind=rp), intent(in) :: u(n), x(n)
    integer, intent(in) :: n

    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: u_GLL, x1, x2
    integer :: i

    do i = 1, n-1
       x1 = x(i)
       x2 = x(i+1)

       ! Check if we are in the interval ]x1, x2[ or if x = x1 or x = x2
       if ( ((x_GLL - x1)*(x_GLL - x2) .lt. 0.0_rp) .or. (abs(x_GLL - x1) .lt. 1d-7) .or. (abs(x_GLL - x2) .lt. 1d-7)) then
          u_GLL = interpolate_1D(x_GLL, x1, u(i), x2, u(i+1))
          !if (x_GLL .lt. 5d-3) print *, x_GLL,x1, x2,u_GLL, u(i), u(i+1)
          return
       end if
    end do


    ! If we reach this point it means we haven't found an interval in which
    !to interpolate
    if (x_GLL .lt. maxval(x) .and. x_GLL .gt. minval(x)) then
       write (log_buf, '(A, F10.6)') "Could not interpolate point ", x_GLL
       call neko_error(log_buf)
    else
       u_GLL = 1.0_rp
    end if

  end function local_interpolation

  ! Linear interpolation between two points (xa,ya) and (xb,yb)
  function interpolate_1D(x, xa, ya, xb, yb) result(y)
    real(kind=rp), intent(in) :: x, xa, ya, xb, yb
    real(kind=rp) :: y

    y = (yb - ya)/(xb - xa)*(x-xa) + ya

  end function interpolate_1D



!=============================================================================!
  subroutine usercheck(t, tstep, u, v, w, p, coef, param )
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: param

    call tripper%update(t)


  end subroutine usercheck
!=============================================================================!

end module user
