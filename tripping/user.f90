module user
  use neko
  use trip
  implicit none

  type(trip_t) :: tripper
  logical :: enable_tripping

  !Things for tripping
  ! Parameters taken from Schlatter & Örlu 2012 doi:10.1017/jfm.2012.324
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

contains


  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_f_vector => userf
    u%fluid_user_f => tripline_org
    u%user_init_modules => initialize
    u%user_finalize_modules => finalize
    u%user_check => usercheck
  end subroutine user_setup

  ! Finalize user variables or external objects
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

  end subroutine finalize

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

  end subroutine initialize

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

    call tripper%apply(u,v,w,j,k,l,e)

  end subroutine tripline_org

  !> User source term
  subroutine userf(f,t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i, j, k, e
    real(kind=rp) :: u, v, w

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
    
  end subroutine userf

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
