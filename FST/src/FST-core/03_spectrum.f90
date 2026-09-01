module spectrum
  use num_types, only: xp
  implicit none

contains

  !> Computes the Von Karman spectrum
  !! @param k Wave number.
  !! @param L Integral length scale.
  !! @param q Turbulence intensity.
  function ek(k,L,q) result(E)
    real(kind=xp), intent(in) :: k, L, q
    real(kind=xp) :: E

    E = 2._xp/3._xp*q*1.606_xp * (k*L)**4._xp * L / &
         (1.350_xp+(k*L)**2._xp)**(17._xp/6._xp)
  end function ek

end module spectrum
