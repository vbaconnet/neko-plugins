module fst_utils
    use num_types, only : rp
    use logger, only : LOG_SIZE, neko_log
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
    real(kind=rp) :: res
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
    res=iy*rm
    idum=mod(ia*idum+ic,m)
    ir(j)=idum

    return
  end function ran2

  ! Use to print one parameter
  subroutine print_param(name, value)
    implicit none

    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: value
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A,A,g0)') name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_param

end module fst_utils