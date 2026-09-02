module fst_utils
    use num_types, only : rp, dp, xp, sp
    use logger, only : LOG_SIZE, neko_log


    interface print_param
      module procedure print_param_int, print_param_sp, print_param_dp, print_param_bool
    end interface print_param

contains

    !> A simple portable random number generator
    !! Requires 32-bit integer arithmetic
    !! Taken from Numerical Recipes, William Press et al.
    !! gives correlation free random numbers but does not have a very large
    !! dynamic range, i.e only generates 714025 different numbers
    !! for other use consult the above
    !! Set idum negative for initialization
    function ran2(idum) result(res)

    implicit none

    integer idum,ir(97),m,ia,ic,iff,iy,j
    real rm
    real(kind=xp) :: res
    parameter (m=714025,ia=1366,ic=150889,rm=1.0_xp/real(m, kind=xp))
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
    res=iy*rm
    idum=mod(ia*idum+ic,m)
    ir(j)=idum

    return
  end function ran2

  ! Use to print one parameter
  subroutine print_param_sp(name, value, fmt)
    implicit none

    character(len=*), intent(in) :: name
    real(kind=sp), intent(in) :: value
    character(len=*), intent(in), optional :: fmt

    character(len=LOG_SIZE) :: log_buf
    character(len=128) :: fmt_

    fmt_ = 'g0'
    if (present(fmt)) fmt_ = trim(fmt)

    write(log_buf, '(A,A,' // trim(fmt_) // ')') name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_param_sp

  ! Use to print one parameter
  subroutine print_param_dp(name, value, fmt)
    implicit none

    character(len=*), intent(in) :: name
    real(kind=dp), intent(in) :: value
    character(len=*), intent(in), optional :: fmt

    character(len=LOG_SIZE) :: log_buf
    character(len=128) :: fmt_

    fmt_ = 'g0'
    if (present(fmt)) fmt_ = trim(fmt)

    write(log_buf, '(A,A,' // trim(fmt_) // ')') name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_param_dp

  ! Use to print one parameter
  subroutine print_param_int(name, value)
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A,A,I5)') name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_param_int

    ! Use to print one parameter
  subroutine print_param_bool(name, value)
    implicit none

    character(len=*), intent(in) :: name
    logical, intent(in) :: value
    character(len=LOG_SIZE) :: log_buf
    character(len=3) :: yes_or_no

    yes_or_no = "no"
    if (value) yes_or_no = "yes"

    write(log_buf, '(A,A,A)') name, ": ", trim(yes_or_no)
    call neko_log%message(log_buf)

  end subroutine print_param_bool
end module fst_utils