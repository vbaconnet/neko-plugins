module global_params
  use logger, only: neko_log, LOG_SIZE
  use num_types, only: rp, dp

  implicit none

  real(kind=rp) :: glb_uinf = 6.0_rp
  real(kind=rp) :: fst_ti = 3.7d-2 ! turbulence intensity
  real(kind=rp) :: fst_il = 11.53d-3 ! integral length scale		

  integer, parameter :: nshells = 80! No of spherical shells
  integer, parameter :: Npmax = 40  ! Npmax  -  Number of points in a shell
  real(kind=rp), parameter :: kstart = 63.77_rp ! smallest wavenumber
  real(kind=rp), parameter :: kend = 1660.0_rp ! largest wavenumber

  integer, parameter :: fst_modes = 2 * nshells * Npmax ! No of freestream modes  
  integer :: shell_modes(nshells) ! Modes saved per shell    

  ! --- These two parameters don't seem to be used
  real(kind=rp) :: kmin = kstart*0.5
  real(kind=rp) :: kmax = kend*2.0
  ! ----------------------------------------------

  real(kind=rp) :: bb(fst_modes, 3)  ! RENAME THESE!!! So confusing

  integer :: k_length
  integer :: shell(fst_modes)
  integer :: seed

  real(kind=rp) :: k_num(fst_modes, 3)
  real(kind=rp) :: k_num_all(fst_modes, 3)
  real(kind=rp) :: u_hat_pn(fst_modes, 3)
  real(kind=rp) :: shell_amp(nshells)

  logical :: write_files = .true.

contains

  ! Use to print one parameter
  subroutine print_param(name, value)
    implicit none

    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: value
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, *) name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_param

  ! Use to print one parameter
  subroutine print_int(name, value)
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, *) name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_int

  subroutine gparams_print_params()
    character(len=LOG_SIZE) :: log_buf

    call print_param("nshells", real(nshells, kind=rp))
    call print_param("Npmax", real(Npmax, kind=rp))
    call print_param("fst_ti", fst_ti)
    call print_param("fst_il", fst_il)

  end subroutine gparams_print_params

  real(kind=rp) FUNCTION ran2(idum)
!
!     A simple portable random number generator
!
!     Requires 32-bit integer arithmetic
!     Taken from Numerical Recipes, William Press et al.
!     gives correlation free random numbers but does not have a very large
!     dynamic range, i.e only generates 714025 different numbers
!     for other use consult the above
!     Set idum negative for initialization
!
      implicit none

      integer idum,ir(97),m,ia,ic,iff,iy,j
      real rm
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      save iff,ir,iy
      data iff /0/

      if (idum.lt.0.or.iff.eq.0) then

!     Initialize
!
         iff=1
         idum=mod(ic-idum,m)
         do j=1,97
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         end do
         idum=mod(ia*idum+ic,m)
         iy=idum
      end if
!
!     Generate random number
!
      j=1+(97*iy)/m
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum

      return
  end function ran2
end module global_params
