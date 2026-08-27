!
!     Philipp Schlatter's routine to generate wavenumbers on a sphere.
!
!----------------------------------------------------------------------

module spec
  use math
  use sphere, only: compute_sphere
  use spectrum, only: ek
  !use global_params
  use num_types, only: rp
  use fst_utils, only : ran2
  use utils, only : neko_error
  use fst_utils, only : print_param
  use logger, only : LOG_SIZE, neko_log, NEKO_LOG_INFO

  implicit none

contains

  subroutine spec_s(Npeff, IL, Tu, U_inf, Npmax, Nshells, k_start, k_end, &
       k_x, k_y, k_z, shell, shell_amp, dlx, dly, dlz, periodic_x, periodic_y, &
       periodic_z, seed, write_file_path, write_files)

    integer, intent(out) :: Npeff
    real(kind=rp), intent(in) :: IL, Tu, U_inf
    integer, intent(in) :: Npmax
    integer, intent(in) :: Nshells
    real(kind=rp), intent(in) :: k_start, k_end
    real(kind=rp), allocatable, intent(inout) :: k_x(:), k_y(:), k_z(:)
    integer, allocatable, intent(inout) :: shell(:)
    real(kind=rp), intent(inout) :: shell_amp(nshells)
    real(kind=rp), intent(in) :: dlx, dly, dlz
    logical, intent(in) :: periodic_x, periodic_y, periodic_z
    integer, intent(inout) :: seed
    character(len=*), intent(in) :: write_file_path
    logical, intent(in) :: write_files

    character(len=LOG_SIZE) :: log_buf
    integer :: shell_modes(nshells) ! Modes saved per shell

    ! integer :: Nsmax
    ! Nsmax = nshells

    !     Npmax  -  Number of points in a shell
    !     Nsmax  -  Number of shells

    real(kind=rp) :: k2

    integer start, Ndk

    integer i,j,l,k

    integer n, Np
    integer z1,z2
    !integer :: seed

    real(kind=rp) :: dk, dkint

    real(kind=rp) :: co(2*Npmax,nshells,3)
    real(kind=rp) :: kk(0:nshells)
    real(kind=rp) :: q_truncated, q_continuous
    real(kind=rp) :: q(Nshells) ! tke in each shell

    real(kind=rp) :: q_theoretical ! Theoretical TKE
    
    Np = Npmax
    dk = (k_end - k_start)/real(nshells-1, kind=rp)
    q_theoretical = (3.0_rp/2.0_rp*(Tu*U_inf)**2.0_rp) ! = 3/2 * Tu**2 * Uinf**2

    !
    ! Generate wavenumbers, this will also give us a definitive value for
    ! Np
    !
    call neko_log%section("Wavenumbers")    
    if (periodic_y .or. periodic_x .or. periodic_z) &
      call neko_log%message("Enforcing periodicity on wavenumbers", &
         lvl=NEKO_LOG_INFO)

    kk(0) = 0.0_rp
    do i=1,nshells
       ! Fill the total wavenumber vector
       kk(i) = k_start + (i-1)*dk ! kk = k_start, k_start+dk, k_start+2dk + ... + k_end
       ! Fill            co(1:Np,i,1), co(1:Np,i,2), co(1:Np,i,3)
       call gen_dodeca_k(co(:,i,1), co(:,i,2), co(:,i,3), &
            kk(i),Np,seed)
       Npeff = Np

       call periodicity_chk(co(:,i,1),co(:,i,2),co(:,i,3), &
           Np,kk(i), dlx, dly, dlz, periodic_x, periodic_y, periodic_z, seed)

    end do

    !
    ! Allocate the arrays here because they will be filled below
    !
    call neko_log%message("Allocating arrays kx,ky,kz", lvl=NEKO_LOG_INFO)
    allocate(k_x(Np*2*nshells))
    allocate(k_y(Np*2*nshells))
    allocate(k_z(Np*2*nshells))
    allocate(shell(Np*2*nshells))

    !     Write wavenumbers to ffst_ile
    if (write_files) then
       call neko_log%message("Creating file " // trim(write_file_path) // &
            '/sphere.dat', lvl=NEKO_LOG_INFO)

       open(file=trim(write_file_path) // '/sphere.dat', unit=10)
       write(10,*) 'energy shell parameters'
       write(10,'(a20,i18)') 'Nshells',nshells
       write(10,'(a20,f18.9)') 'k_start',k_start
       write(10,'(a20,f18.9)') 'k_end',k_end
       write(10,'(a20,i18)') 'Np',Np
       write(10,*) 'isotropic coordinates'
       write(10,'(2a5,3a18)') 'i','j','x','y','z'
    endif

    !
    ! Remove mode (0,0,0) and mirror modes in x axis
    !

    do i=1,nshells

       ! add second dodecaeder mirrored at (x)-axis
       do j=Np+1,2*Np
          co(j,i,1) = co(j-Np,i,1)
          co(j,i,2) = -co(j-Np,i,2)
          co(j,i,3) = -co(j-Np,i,3)
       end do

       if (write_files) then
          do j=1,2*Np
             write(10,'(2i5,3e18.9)') i,j,co(j,i,1),co(j,i,2),co(j,i,3)
          end do
       endif
    end do

    if (write_files) close(10) ! sphere.dat
    
    !
    ! Remove mode (0,0,0) if it exists, and assign the result to our
    ! arrays kx, ky, kz
    !
    call neko_log%message("Removing mode (0,0,0) if it exists", &
      lvl=NEKO_LOG_INFO)
    shell_modes = 0
    z1=0
    z2=0
    l=0

    do i=1, nshells
       do j=1, 2*Np
          !         If some modes need to be removed.
          !         Currently all modes except (0,0,0) are preserved.
          if (co(j,i,1).eq.0.0.and.co(j,i,2).eq.0.0.and.co(j,i,3).eq.0.0) then
             z2=z2+1
             continue
          else
             z1=z1+1
             
             shell_modes(i)=shell_modes(i)+1
             
             l=l+1

             k_x(l) = co(j,i,1)
             k_y(l) = co(j,i,2)
             k_z(l) = co(j,i,3)

             shell(l) = i

          endif ! if (.not.(0,0,0))
       end do ! j=1,2*Np
    end do ! i=1,nshells

    if (z2 .ne. 0) &
      call neko_log%message('(0,0,0) wavenumber removed', lvl=NEKO_LOG_INFO)

    write(log_buf, '(A,I0,A,I0,A)') 'Saved ',z1,' of ',z1+z2, ' fst modes.'
    call neko_log%message(log_buf, lvl=NEKO_LOG_INFO)

    call neko_log%end_section()

    !
    ! Generate amplitudes
    !
    call neko_log%section('Amplitudes')

    write (log_buf, '(A,F10.6)') "Theoretical TKE, (q) : ", q_theoretical
    call neko_log%message(log_buf, lvl=NEKO_LOG_INFO)

    !  ------ integrate the energy spectrum (mimics continuous integral) ---
    Ndk = 5000 ! just a large no of points on the spectrum
    dkint = (k_end-k_start)/float(Ndk)

    ! Include the bounds first
    q_continuous = ek(k_start, IL, 1.0_rp) + ek(k_end, IL, 1.0_rp)

    do i=1,Ndk-1
       q_continuous = q_continuous + ek(k_start + i*dkint, IL, 1._rp)
    end do
    q_continuous = q_continuous*dkint
    write (log_buf, '(A,F10.6)') 'Truncated integral of spectrum :', q_continuous
    call neko_log%message(log_buf, lvl=NEKO_LOG_INFO)
    ! ------------------------------------------------------------------------

    ! ----- integrate the energy spectrum with nshells points ----------------
    ! This is the "discretized" energy.
    q_truncated = 0.0_rp
    do i=1,nshells
       q_truncated = q_truncated + ek(k_start + (i-1)*dk, IL, 1._rp)
    end do
    q_truncated = q_truncated*dk
    write (log_buf, '(A,I0,A,F10.6)') 'Truncated, discrete integral on ', &
      nshells, ' shells, (q_hat) :' , q_truncated
    call neko_log%message(log_buf, lvl=NEKO_LOG_INFO)
    ! -------------------------------------------------------------------------
    
    call print_param("Ratio q / q_hat", q_theoretical/q_truncated, fmt='F10.6')

    !
    ! Generate amplitudes
    !
    do i = 1, Nshells

       ! Generate local TKE
       q(i) = ek(kk(i), IL, q_theoretical/q_truncated)

       shell_amp(i) = sqrt(2.0_rp * q(i)*dk * 2.0_rp / &
            (real(shell_modes(i), kind=rp)))
       
    end do

    call neko_log%end_section()


    return
  end subroutine spec_s

  !----------------------------------------------------------------------

  !> Based on 3d wavenumbers, makes the wavenumbers in the direction
  !! kp periodic, based on the length Lp.
  !! kp is the array of wavenumbers in the periodic direction. kp is filled
  !! with wavenumbers that are multiple of 2pi/Lp, so kp = n*2pi/Lp
  subroutine make_periodic_1D(k1, k2, kp, np, K_total, Lp, seed)
    real(kind=rp), intent(inout) :: k1(:), k2(:), kp(:)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: K_total
    real(kind=rp), intent(in) :: Lp
    integer, intent(inout) :: seed

    integer :: nmax, nmin, n_j, n_j_signed, j
    real(kind=rp) :: twopi_over_L, rtmp, flip, K_total_sq

    twopi_over_L = 2.0_rp * pi / Lp
    K_total_sq = K_total**2.0_rp

    !
    ! First, check if we can fix at least one wavenumber in the direction Lp
    ! nmax gives the number of wavenumbers we can fit in the length Lp
    ! based on the total radius K_total.

    ! n * 2 pi / Lp gives the discrete wavenumbers.
    nmax = floor(K_total/twopi_over_L)
    nmin = 1

    ! Can at least one multiple of twopi_over_L fit inside this total wavenumber?
    if (nmax .lt. nmin) then
       call print_param('K_total:', K_total)
       call print_param('2 pi / L:', twopi_over_L)
       call neko_log%message("2 pi / L must be smaller than K_total!")
       call neko_error('Increase minimum total wave number!')
    endif

    do j = 1, np

       ! how many multiples of 2pi/L can fit in this direction
       n_j = floor( abs(kp(j)) / twopi_over_L)
       n_j_signed = sign(1.0_rp, kp(j)) * floor( abs(kp(j)) / twopi_over_L)

       if (n_j .gt. nmax) then
          n_j_signed = n_j_signed - sign(1.0_rp, kp(j))

       elseif (n_j .eq. 0) then
          ! Force to not be zero
          n_j_signed = n_j_signed + sign(1.0_rp, kp(j))
       endif

       ! Set the discrete wavenumber
       kp(j) = real(n_j_signed, kind=rp) * twopi_over_L

       !
       ! Now we need to adjust the other wavenumbers to make sure we still get
       ! K_tot^2 = k_1^2 + k_2^2 + k_p^2.
       ! The way we do it is to fix the value of k1 or k2 (decided randomly),
       ! while the other value is recomputed as e.g.
       ! k1^2 = K_tot^2 - k_2^2 - k_p^2
       !
       flip = ran2(seed) ! coin toss

       ! > 0.5 means we fix k1 and recompute k2
       if (flip .gt. 0.5_rp) then

          !k1(j) = k1(j)       ! k1 stays the same
          rtmp = K_total_sq - k1(j)**2.0_rp - kp(j)**2.0_rp

          if (rtmp .gt. 1) then
             k2(j) = sign(1.0_rp, k2(j))*sqrt(rtmp)
          else
             rtmp = sqrt((K_total_sq - kp(j)**2.0_rp)/2.0_rp)
             k1(j) = sign(1.0_rp,k1(j))*rtmp
             k2(j) = sign(1.0_rp,k2(j))*rtmp
          endif

          ! < 0.5 means we fix k2 and recompute k1
       else

          !k2(j) = k2(j)       ! k2 stays the same
          rtmp = K_total_sq - kp(j)**2.0_rp - k2(j)**2.0_rp

          if (rtmp .gt. 1) then
             k1(j) = sign(1.0_rp, k1(j)) * sqrt(rtmp)
          else
             rtmp = sqrt((K_total_sq - kp(j)**2.0_rp)/2.0_rp)
             k1(j) = sign(1.0_rp,k1(j))*rtmp
             k2(j) = sign(1.0_rp,k2(j))*rtmp
          endif

       endif ! flip
    enddo ! j=1,Np

  end subroutine make_periodic_1D

  !> Based on 3d wavenumbers, makes the wavenumbers in the directions
  !! kp1 and kp2 periodic, based on the lengths Lp1 and Lp2.
  !! See make_periodic_1D for more details.
  subroutine make_periodic_2D(k1, kp1, kp2, np, K_total, L1, L2)
    real(kind=rp), intent(inout) :: k1(:), kp1(:), kp2(:)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: K_total
    real(kind=rp), intent(in) :: L1, L2

    integer :: nmax, nmin, n_j1, n_j1_signed, j
    integer :: n_j2, n_j2_signed, n_signed_cand
    real(kind=rp) :: twopi_over_L1, twopi_over_L2, rtmp, K_total_sq, flip
    logical :: valid_config

    twopi_over_L1 = 2.0_rp * pi / L1
    twopi_over_L2 = 2.0_rp * pi / L2
    K_total_sq = K_total**2.0_rp

    !
    ! First, check if we can fix at least one wavenumber in the direction Lp
    ! nmax gives the number of wavenumbers we can fit in the length Lp
    ! based on the total radius K_total.

    ! ---------------------------------------------------------
    ! Check in the direction 1

    ! n * 2 pi / Lp gives the discrete wavenumbers.
    nmax = floor(K_total/twopi_over_L1)
    nmin = 1

    ! Can at least one multiple of twopi_over_L fit inside this total wavenumber?
    if (nmax .lt. nmin) then
       call print_param('K_total:', K_total)
       call print_param('2 pi / L1:', twopi_over_L1)
       call neko_log%message("2 pi / L must be smaller than K_total!")
       call neko_error('Increase minimum total wave number!')
    endif

    ! ---------------------------------------------------------
    ! Check in the direction 2

    ! n * 2 pi / Lp gives the discrete wavenumbers.
    nmax = floor(K_total/twopi_over_L2)
    nmin = 1

    ! Can at least one multiple of twopi_over_L fit inside this total wavenumber?
    if (nmax .lt. nmin) then
       call print_param('K_total:', K_total)
       call print_param('2 pi / L2:', twopi_over_L2)
       call neko_log%message("2 pi / L must be smaller than K_total!")
       call neko_error('Increase minimum total wave number!')
    endif

    !
    ! Now, compute the discrete wavenumbers
    !
    do j = 1, np

       ! ---- Set the discrete wavenumber in direction 1
       n_j1 = floor( abs(kp1(j)) / twopi_over_L1)
       n_j1_signed = sign(1.0_rp, kp1(j)) * floor( abs(kp1(j)) / twopi_over_L1)

       if (n_j1 .gt. nmax) then
          n_j1_signed = n_j1_signed - sign(1.0_rp, kp1(j))

       elseif (n_j1 .eq. 0) then
          ! Force to not be zero
          n_j1_signed = n_j1_signed + sign(1.0_rp, kp1(j))
       endif

       kp1(j) = real(n_j1_signed, kind=rp) * twopi_over_L1

       ! ---- Set the discrete wavenumber in direction 2
       n_j2 = floor( abs(kp2(j)) / twopi_over_L2)
       n_j2_signed = sign(1.0_rp, kp2(j)) * floor( abs(kp2(j)) / twopi_over_L2)

       if (n_j2 .gt. nmax) then
          n_j2_signed = n_j2_signed - sign(1.0_rp, kp2(j))

       elseif (n_j2 .eq. 0) then
          ! Force to not be zero
          n_j2_signed = n_j2_signed + sign(1.0_rp, kp2(j))
       endif

       kp2(j) = real(n_j2_signed, kind=rp) * twopi_over_L2

       !
       ! Now we need to adjust the other wavenumbers to make sure we still get
       ! K_tot^2 = k_1^2 + k_p1^2 + k_p2^2, so
       ! k_1**2 = K_tot**2 - k_p1**2 - k_p2**2
       ! and we also make sure to keep the sign

       ! warning: rtmp can be negative! In that case we need to adjust the values
       ! of nj1 and/or nj2 until we can get something.
       rtmp = K_total_sq - kp1(j)**2.0_rp - kp2(j)**2.0_rp
       valid_config = (rtmp .gt. 1.0_rp)

       do while (.not. valid_config)
         
          ! Always reduce the component that is highest
          if (kp1(j) .gt. kp2(j)) then
             n_j1_signed = n_j1_signed - sign(1.0_rp, kp1(j))
             kp1(j) = real(n_j1_signed, kind=rp) * twopi_over_L1
          else
             n_j2_signed = n_j2_signed - sign(1.0_rp, kp2(j))
             kp2(j) = real(n_j2_signed, kind=rp) * twopi_over_L2
          end if

          rtmp = K_total_sq - kp1(j)**2.0_rp - kp2(j)**2.0_rp
          valid_config = (rtmp .gt. 0.0_rp)

       end do

       k1(j) = sign(1.0_rp, k1(j)) * sqrt( rtmp )

    enddo

  end subroutine make_periodic_2D


  subroutine periodicity_chk(kx, ky, kz, np, kk, dlx, dly, dlz, ifxp, ifyp, &
       ifzp, seed)
    real(kind=rp), intent(inout) :: kx(:), ky(:), kz(:)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: kk
    real(kind=rp), intent(in) :: dlx,dly,dlz
    logical, intent(in) :: ifxp,ifyp,ifzp
    integer, intent(inout) :: seed

    ! integer :: i,j
    ! integer :: nmax,nmin,kn

    ! real(kind=rp) :: flip, k2
    ! real(kind=rp) :: rtmp
    ! k2 = kk**2

    logical :: periodic_1d, periodic_2d

    if (ifxp .and. ifyp .and. ifzp) &
         call neko_error("Only two periodic directions supported")

    ! Whether we have to apply correction for periodicity in two directions
    periodic_2d = (ifxp .and. ifyp) .or. &
         (ifxp .and. ifzp) .or. &
         (ifyp .and. ifzp)

    if (periodic_2d) then

       if (ifxp .and. ifyp) then
          call make_periodic_2D(kz, kx, ky, np, kk, dlx, dly)
       elseif (ifxp .and. ifzp) then
          call make_periodic_2D(ky, kz, kx, np, kk, dlz, dlx)
       elseif (ifyp .and. ifzp) then
          call make_periodic_2D(kx, ky, kz, np, kk, dly, dlz)
       end if

    else ! Apply 1D correction

       if (ifxp) then
          call make_periodic_1D(ky, kz, kx, np, kk, dlx, seed)
       elseif (ifyp) then
          call make_periodic_1D(kz, kx, ky, np, kk, dly, seed)
       elseif (ifzp) then
          call make_periodic_1D(kx, ky, kz, np, kk, dlz, seed)
       end if

    end if

    ! OLD WAY OF DOING THE PERIODICITY CHECKS BELOW
    !     Add periodicity check
    ! if (ifxp) then
    !   call neko_log%message('Checking periodicity in x')
    !   nmax = floor(kk*dlx/(2.0*pi))
    !   nmin = 1
    !   if (nmax.lt.nmin) then
    !     call neko_log%message('Check allowed wavenumbers in FST')
    !     call print_param('nmax:', real(nmax, kind=rp))
    !     call print_param('nmin:', real(nmin, kind=rp))
    !     call print_param('k   :', kk)
    !     call exit
    !   endif

    !   do j=1,np
    !     !          kn = nint(kx(j)*dlx/(2.0*pi))
    !     kn = sign(1.0_rp, kx(j)) * floor(abs(kx(j))*dlx/(2.0_rp*pi))  ! always
    !     ! make k
    !     ! smaller

    !     if (abs(kn).gt.nmax) then
    !       kn=kn-sign(1.0_rp,kx(j))
    !     elseif (abs(kn).eq.0) then
    !       kn=kn+sign(1.0_rp,kx(j))
    !     endif
    !     kx(j)=real(kn)*2.0_rp*pi/dlx

    !     flip = ran2(seed)            ! coin toss
    !     if (flip.gt.0.5_rp) then
    !       ky(j) = ky(j)       ! ky stays the same
    !       rtmp = k2-ky(j)**2.0_rp-kx(j)**2.0_rp
    !       if (rtmp.gt.1) then
    !         kz(j) = sign(1.0_rp,kz(j))*sqrt(rtmp)
    !       else
    !         rtmp = sqrt((k2-kx(j)**2.0_rp)/2.0_rp)
    !         ky(j) = sign(1.0_rp,ky(j))*rtmp
    !         kz(j) = sign(1.0_rp,kz(j))*rtmp
    !       endif
    !     else
    !       kz(j) = kz(j)       ! kz stays the same
    !       rtmp = k2-kx(j)**2.0_rp-kz(j)**2.0_rp
    !       if (rtmp.gt.1) then
    !         ky(j) = sign(1.0_rp,ky(j))*sqrt(rtmp)
    !       else
    !         rtmp = sqrt((k2-kx(j)**2.0_rp)/2.0_rp)
    !         ky(j) = sign(1.0_rp,ky(j))*rtmp
    !         kz(j) = sign(1.0_rp,kz(j))*rtmp
    !       endif
    !     endif       ! flip
    !   enddo         ! j=1,Np

    ! end if

    ! if (ifyp) then
    !   call neko_log%message('Checking periodicity in y')
    !   nmax = floor(kk*dly/(2.0_rp*pi))
    !   nmin = 1
    !   if (nmax.lt.nmin) then
    !     call neko_log%message('Check allowed wavenumbers in FST')
    !     call print_param('nmax:', real(nmax, kind=rp))
    !     call print_param('nmin:', real(nmin, kind=rp))
    !     call print_param('k   :', kk)
    !     call exit
    !   endif

    !   do j=1,np
    !     !          kn = nint(ky(j)*dly/(2.0*pi))
    !     kn = sign(1.0_rp,ky(j))*floor(abs(ky(j))*dly/(2.0*pi))  ! always
    !     ! make k
    !     ! smaller


    !     if (abs(kn).gt.nmax) then
    !       kn=kn-sign(1.0_rp,ky(j))
    !       elseif (abs(kn).eq.0) then
    !         kn=kn+sign(1.0_rp,ky(j))
    !       endif
    !       ky(j)=real(kn)*2.0*pi/dly

    !       flip = ran2(seed)            ! coin toss
    !       if (flip.gt.0.5) then
    !         kz(j) = kz(j)       ! kz stays the same
    !         rtmp = k2-ky(j)**2.-kz(j)**2.
    !         if (rtmp.gt.1) then
    !           kx(j) = sign(1.0_rp,kx(j))*sqrt(rtmp)
    !         else
    !           rtmp = sqrt((k2-ky(j)**2.)/2.)
    !           kx(j) = sign(1.0_rp,kx(j))*rtmp
    !           kz(j) = sign(1.0_rp,kz(j))*rtmp
    !         endif
    !       else
    !         kx(j) = kx(j)       ! kx stays the same
    !         rtmp = k2-ky(j)**2.-kx(j)**2.
    !         if (rtmp.gt.1) then
    !           kz(j) = sign(1.0_rp,kz(j))*sqrt(rtmp)
    !         else
    !           rtmp = sqrt((k2-ky(j)**2.)/2.)
    !           kx(j) = sign(1.0_rp,kx(j))*rtmp
    !           kz(j) = sign(1.0_rp,kz(j))*rtmp
    !         endif
    !       endif       ! flip
    !     enddo         ! j=1,Np

    ! end if

    ! if (ifzp) then
    !       call neko_log%message('Checking periodicity in z')
    !       nmax = floor(kk*dlz/(2.0*pi))
    !       nmin = 1
    !       if (nmax.lt.nmin) then
    !         call neko_log%message('Check allowed wavenumbers in FST')
    !         call print_param('nmax:', real(nmax, kind=rp))
    !         call print_param('nmin:', real(nmin, kind=rp))
    !         call print_param('k   :', kk)
    !         call exit
    !       endif

    !       do j=1,np
    !         kn = sign(1.0_rp,kz(j))*floor(abs(kz(j))*dlz/(2.0*pi))  ! always
    !         ! make k
    !         ! smaller
    !         if (abs(kn).gt.nmax) then
    !           kn=kn-sign(1.0_rp,kz(j))
    !           elseif (abs(kn).eq.0) then
    !             kn=kn+sign(1.0_rp,kz(j))
    !           endif
    !           kz(j)=real(kn)*2.0*pi/dlz

    !           flip = ran2(seed)            ! coin toss
    !           if (flip.gt.0.5) then
    !             kx(j) = kx(j)       ! kx stays the same
    !             rtmp = k2-kx(j)**2.0_rp-kz(j)**2.0_rp
    !             if (rtmp.gt.1) then
    !               ky(j) = sign(1.0_rp,ky(j))*sqrt(rtmp)
    !             else
    !               rtmp = sqrt((k2-kz(j)**2.)/2.)
    !               kx(j) = sign(1.0_rp,kx(j))*rtmp
    !               ky(j) = sign(1.0_rp,ky(j))*rtmp
    !             endif
    !           else
    !             ky(j) = ky(j)       ! ky stays the same
    !             rtmp = k2-ky(j)**2.0_rp-kz(j)**2.0_rp
    !             if (rtmp.gt.1) then
    !               kx(j) = sign(1.0_rp,kx(j))*sqrt(rtmp)
    !             else
    !               rtmp = sqrt((k2-kz(j)**2.0_rp)/2.0_rp)
    !               kx(j) = sign(1.0_rp,kx(j))*rtmp
    !               ky(j) = sign(1.0_rp,ky(j))*rtmp
    !             endif
    !           endif       ! flip

    !         enddo         ! j=1,Np

    ! endif           ! ifxp

    return
  end subroutine periodicity_chk
  !----------------------------------------------------------------------

  subroutine gen_bounded_k(Np,kx,ky,kz,kk,kmin,kmax,seed)

    implicit none

    integer Np
    real(kind=rp) :: kx(1),ky(1),kz(1)
    real(kind=rp) :: kk
    real(kind=rp) :: kmax,kmin
    real(kind=rp) :: theta,phi

    integer i,seed

    real(kind=rp) :: twopi
    logical inbounds

    twopi = 2.0*4.0*atan(1.0)

    do i=1,Np
       inbounds=.false.

       do while (.not.inbounds)
          theta=twopi*ran2(seed)
          phi=twopi*ran2(seed)

          ky(i)=kk*sin(theta)
          kx(i)=kk*cos(theta)*cos(phi)
          kz(i)=kk*cos(theta)*sin(phi)

          !          write(6,*) ky(i),kx(i),kz(i),kmin

          if ((abs(kx(i)).gt.kmin).and. &
               (abs(ky(i)).gt.kmin).and. &
               (abs(kz(i)).gt.kmin)) then
             inbounds=.true.
          endif
          !         Don't need check for kmax since kk is projected on a sphere
          !         Hence max amplitude == kk

       enddo
    enddo

    return
  end subroutine gen_bounded_k
  !----------------------------------------------------------------------

  !> NOTE: This modifies the value of Np!
  subroutine gen_dodeca_k(kx, ky, kz, K_tot, Np, seed)

    real(kind=rp), intent(inout) :: kx(:),ky(:),kz(:)
    real(kind=rp), intent(in) :: K_tot
    integer, intent(inout) :: Np
    integer, intent(in) :: seed

    real(kind=rp) :: rotx,roty,rotz

    rotx = ran2(seed)*2.0_rp*pi
    roty = ran2(seed)*2.0_rp*pi
    rotz = ran2(seed)*2.0_rp*pi
    call compute_sphere(Np, kx, ky, kz, K_tot, rotx, roty, rotz, .false.)

    return
  end subroutine gen_dodeca_k

  !----------------------------------------------------------------------
  real(kind=rp) function vlamin(vec,n)
    real(kind=rp) :: VEC(1)
    integer, intent(in) :: n
    integer :: i
    real(kind=rp) :: TMIN
    TMIN = 99.0E20

    do 100 I=1,N
       TMIN = MIN(TMIN,abs(VEC(I)))
100 CONTINUE
    VLAMIN = TMIN
    return
  end function vlamin

  real(kind=rp) function vlamax(vec,n)
    real(kind=rp) :: VEC(1)
    integer, intent(in) :: n
    integer :: i
    real(kind=rp) :: TMAX
    TMAX = 99.0E20

    do 100 I=1,N
       TMAX = MAX(TMAX,abs(VEC(I)))
100 CONTINUE
    VLAMAX = TMAX
    return
  end function vlamax
  !-----------------------------------------------------------------------

end module spec
