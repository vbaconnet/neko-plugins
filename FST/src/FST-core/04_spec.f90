!
!     Philipp Schlatter's routine to generate wavenumbers on a sphere.
!
!---------------------------------------------------------------------- 

module spec
  use math
  use sphere
  use spectrum, only: ek
  use global_params
  use num_types, only: rp
  use utils, only : neko_error
  
  implicit none
  
  contains
  
  subroutine spec_s(dlx, dly, dlz, periodic_x, periodic_y, periodic_z)
    real(kind=rp), intent(out) :: dlx, dly, dlz
    logical, intent(in) :: periodic_x, periodic_y, periodic_z
    
    character(len=LOG_SIZE) :: log_buf
    
    
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
    
    real(kind=rp) :: dmin,il,dkint
    
    real(kind=rp) :: co(2*Npmax,nshells,3)
    real(kind=rp) :: kk(0:nshells),q(nshells),dk(nshells)
    integer :: lu(nshells)
    real(kind=rp) :: tke_tot         ! 
    real(kind=rp) :: tke_tot1        ! 
    real(kind=rp) :: tke_shell(nshells)
    
    real(kind=rp) :: shell_energy
    real(kind=rp) :: tke_scaled
    
    real(kind=rp) :: vlsum        ! function
    
    real(kind=rp) :: kxmin,kxmax,kymin,kymax,kzmin,kzmax
    real(kind=rp) :: ktmp
    
    !---------------------------------------- 
    
    
    call print_param('integral length scale', fst_il)
    Np = Npmax
    
    tke_scaled = (3.0/2.*(fst_ti*glb_uinf)**2) ! = 3/2 * Tu**2 * Uinf**2
    
    !     Just initializing   
    kxmax = 1.0E-20
    kxmin = 1.0E+20
    
    kymax = 1.0E-20
    kymin = 1.0E+20
    
    kzmax = 1.0E-20
    kzmin = 1.0E+20
    
    !     spectrum discretization
    !      nshells=nshells
    
    !      write(6,*) 'FST - Largest wavenumber:', 2.0*pi/kstart
    !      write(6,*) 'FST - Smallest wavenumber:',2.0*pi/kend
    
    !  ------ integrate the energy spectrum (mimics continuous integral) ---
    Ndk   = 5000    ! just a large no of points on the spectrum
    dkint = (kend-kstart)/float(Ndk)
    
    tke_tot1 = (ek(kstart,fst_il,1._rp) + ek(kend,fst_il,1._rp))
    do i=1,Ndk-1
      tke_tot1 = tke_tot1 + ek(kstart + i*dkint, fst_il, 1._rp)
    end do
    tke_tot1 = tke_tot1*dkint
    call print_param('FST - integrated energy in spectrum ',tke_tot1)
    ! ------------------------------------------------------------------------
    
    ! ----- integrate the energy spectrum with nshells points ----------------
    dkint = (kend - kstart)/real(nshells-1)
    tke_tot = 0.
    do i=1,nshells
      tke_tot = tke_tot + ek(kstart + (i-1)*dkint,fst_il,1._rp)
    end do
    tke_tot = tke_tot*dkint
    write (log_buf, *) 'FST - discretized on ', nshells, ' shells :' , tke_tot
    call neko_log%message(log_buf)
    ! -------------------------------------------------------------------------
    
    !     Write wavenumbers to ffst_ile
    if (write_files) then
      open(file='sphere.dat', unit=10)
      
      write(10,*) 'energy shell parameters'
      write(10,'(a20,i18)') 'Nshells',nshells
      write(10,'(a20,f18.9)') 'kstart',kstart
      write(10,'(a20,f18.9)') 'kend',kend
      write(10,'(a20,i18)') 'Np',Np
      write(10,*) 'isotropic coordinates'
      write(10,'(2a5,3a18)') 'i','j','x','y','z'
    endif 
    
    kk(0) =0.
    !     compute the coordinates using two dodecaeder
    
    tke_tot1 = 0.
    !seed = -143
    
    call print_param("Truncated TKE",tke_scaled/tke_tot)
    
    do i=1,nshells
      
      k2 = ( kstart + (i-1)*(kend - kstart)/real(nshells-1, kind=rp) )**2
      kk(i) = sqrt(k2) ! kk = kstart, kstart+dk, kstart+2dk + ... + kend
      dk(i) = (kend - kstart)/real(nshells-1, kind=rp)
      
      q(i) = ek(kk(i),fst_il,tke_scaled/tke_tot)   ! 1/tke_tot so the total
      ! truncated energy = tke_scaled
      tke_shell(i) = q(i)*dk(i)
      tke_tot1 = tke_tot1 + tke_shell(i)
      
      ! Fill            co(1:Np,i,1), co(1:Np,i,2), co(1:Np,i,3)
      call gen_dodeca_k(co(1,i,1), co(1,i,2), co(1,i,3), &
      kk(i),Np,seed)
      
      ! Recompute wavenumbers in the periodic directions
      call periodicity_chk(co(1,i,1),co(1,i,2),co(1,i,3), &
      Np,kk(i),dlx,dly,dlz, periodic_x, periodic_y, periodic_z)
      
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
      
      !     Get smallest and largest fst modes in x,y,z        
      ! ktmp = vlamax(co(1,i,1),2*Np)
      ktmp = vlmax(co(1,i,1),2*Np)
      !        write(6,*) 'ktmp,kmax',ktmp,kxmax
      if (ktmp.gt.kxmax) kxmax = ktmp
      
      ktmp = vlamin(co(1,i,1),2*Np)
      !        write(6,*) 'ktmp,kmin',ktmp,kxmin
      if (ktmp.lt.kxmin) kxmin = ktmp
      
      ktmp = vlmax(co(1,i,2),2*Np)
      if (ktmp.gt.kymax) kymax = ktmp
      
      ktmp = vlamin(co(1,i,2),2*Np)
      if (ktmp.lt.kymin) kymin = ktmp
      
      ktmp = vlmax(co(1,i,3),2*Np)
      if (ktmp.gt.kzmax) kzmax = ktmp
      
      ktmp = vlamin(co(1,i,3),2*Np)
      if (ktmp.lt.kzmin) kzmin = ktmp
      !--------------------------------------------------
    end do            ! 1,Nshells
    !      write(6,'(A15,1x,E15.8E2)') 'q-shell total:',
    !     &      tke_tot1
    
    if (write_files) close(10) ! sphere.dat
    
    z1=0
    z2=0
    do i=1,nshells
      lu(i)=0
      shell_modes(i)=0
    end do
    
    l=0
    !      open(ffst_ile='chosen_coordinates.dat',unit=12)
    !      write(12,'(A21)') '# choosen coordinates'
    !      write(12,'(A10,2x,3(A18,2x))') '# Shell No', 'kx', 'ky', 'kz'
    do i=1, nshells
      do j=1, 2*Np
        !         If some modes need to be removed.
        !         Currently all modes except (0,0,0) are preserved.
        if (co(j,i,1).eq.0.0.and.co(j,i,2).eq.0.0.and.co(j,i,3).eq.0.0) then
          z2=z2+1
          continue
        else
          z1=z1+1
          
          ! write to chosen_coordinates.dat
          !              write(12,'(i10,2x,3(E18.9E2,2x))') i,
          !     &           co(j,i,1), co(j,i,2), co(j,i,3)
          
          lu(i)=lu(i)+1           ! no of modes in each shell
          shell_modes(i)=shell_modes(i)+1
          
          l=l+1
          do k=1,3
            k_num(l,k) = co(j,i,k)
            k_num_all(l,k) = co(j,i,k)
          enddo
          k_length = l
          shell(l) = i
          !              shell2(l) = i
          
        endif           ! if (.not.(0,0,0))
      end do          ! j=1,2*Np
    end do            ! i=1,nshells
    ! write(6,*) 'FST - (0,0,0) wavenumber removed'
    call neko_log%message('FST - (0,0,0) wavenumber removed')
    
    write(log_buf, *) 'Saved ',z1,' of ',z1+z2, ' fst modes.'
    call neko_log%message(log_buf)
    !      close(12)
    
    !     determine scaling according to the energy spectrum
    !      write(6,2014) 'Shell No', 'Amp','No modes','Tot Energy'
    ! 2014 format(A10,3x,A15,3x,A10,3x,A15)
    
    tke_tot1 = 0.
    do i=1,nshells
      !write(*,*) "YEYEY", tke_shell(i), shell_modes(i)
      shell_amp(i) = sqrt(2.*tke_shell(i)*2./ &
      (real(shell_modes(i), kind=rp)))
      
      shell_energy = real(shell_modes(i), kind=rp)* &
      ((shell_amp(i)**2.))/2.
      
      tke_tot1 = tke_tot1 + shell_energy
      
      !         write(6,2013) i, shell_amp(i),shell_modes(i),
      !     &      shell_energy
    end do
    
    ! 2013 format(i10,3x,E15.8E2,3x,i10,3x,E15.8E2)
    
    !      open(ffst_ile='force.dat',unit=11)
    !      write(11,*) 'energy spectrum parameters'
    !      call writedat(11)
    !      call hline(11)
    !      write(11,'(a20,i18)') 'nshells',nshells
    !      write(11,'(a20,i18)') 'Np',Np
    !      write(11,'(a20,i18)') 'Neig',z1
    !      write(11,'(a20,f18.9)') 'kstart',kstart
    !      write(11,'(a20,f18.9)') 'kend',kend
    !      write(11,'(a20,f18.9)') 'Lint',fst_il
    !      write(11,'(a20,f18.9)') 'tke_tot',tke_tot
    !!      write(11,'(a20,f18.9)') 'width',width
    !      write(11,'(a20,f18.9)') 'etastart',etastart
    !      write(11,'(a20,f18.9)') 'ymax',ymax
    !      write(11,'(a20,f18.9)') 'tke_tot',1./tke_tot*1.5
    !      write(11,'(a20,f18.9)') 'tke_tot1',1./tke_tot1*1.5
    
    !      write(11,'(5a18)') '#Shell','#perShell','k','dk','E(k)'
    
    !      start = 1.
    !      do i=1,nshells
    !         write(*,*) 'Scaling shell number ',i
    !         write(*,*) '  k : ',kk(i)
    !         write(*,*) '  q : ',q(i)
    !         write(*,*) '  # : ',lu(i)
    !         write(*,*) ' dk : ',kk(3)-kk(2)
    
    !         write(11,'(2i5,3E18.9)') i,lu(i),kk(i),kk(3)-kk(2), q(i)
    
    !         scf = 1./200.
    !         scf = q(i)*(kk(3)-kk(2))/lu(i)         ! what is this?
    
    !         shell1(i) = i                         ! not sure of the usage
    !         nshell(i)= lu(i)
    !         wn(i) = kk(i)
    !         dwn(i) = kk(3)-kk(2)
    !         spect(i) = q(i)
    !         spect(i) = q(i)*dk(i)/tke_tot             ! prabal.
    
    !      end do
    
    !      close(11)
    
    write (log_buf, *) 'FST - ',k_length,'wavenumbers generated'
    call neko_log%message(log_buf)
    ! 2012 format(A7,1x,i5,1x,A21)
    
    call print_param('FST - Largest wavelength in x',  2.0*pi/kxmin)
    call print_param('FST - Smallest wavelength in x', 2.0*pi/kxmax)
    call print_param('FST - Largest wavelength in y',  2.0*pi/kymin)
    call print_param('FST - Smallest wavelength in y', 2.0*pi/kymax)
    call print_param('FST - Largest wavelength in z',  2.0*pi/kzmin)
    call print_param('FST - Smallest wavelength in z', 2.0*pi/kzmax)
    
    ! write(6,2015) 'FST - Largest wavelength in x',  2.0*pi/kxmin! ,kxmin
    ! write(6,2016) 'FST - Smallest wavelength in x', 2.0*pi/kxmax! ,kxmax
    ! write(6,2015) 'FST - Largest wavelength in y',  2.0*pi/kymin! ,kymin
    ! write(6,2016) 'FST - Smallest wavelength in y', 2.0*pi/kymax! ,kymax
    ! write(6,2015) 'FST - Largest wavelength in z',  2.0*pi/kzmin! ,kzmin
    ! write(6,2016) 'FST - Smallest wavelength in z', 2.0*pi/kzmax! ,kzmax
    ! 2015 format(A30,6x,E13.5E2)  
    ! 2016 format(A31,5x,E13.5E2)  
    
    return
  end subroutine spec_s
  
  !---------------------------------------------------------------------- 
  
  !> Based on 3d wavenumbers, makes the wavenumbers in the direction
  !! kp periodic, based on the length Lp.
  !! kp is the array of wavenumbers in the periodic direction. kp is filled
  !! with wavenumbers that are multiple of 2pi/Lp, so kp = n*2pi/Lp
  subroutine make_periodic_1D(k1, k2, kp, np, K_total, Lp)
    real(kind=rp), intent(inout) :: k1(1), k2(1), kp(1)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: K_total
    real(kind=rp), intent(in) :: Lp
    
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
        flip = ran2(seed)            ! coin toss
        
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
          !print *, rtmp

          if (rtmp .gt. 1) then                 
            k1(j) = sign(1.0_rp, k1(j)) * sqrt(rtmp)
          else
            print *, rtmp
            rtmp = sqrt((K_total_sq - kp(j)**2.0_rp)/2.0_rp)
            k1(j) = sign(1.0_rp,k1(j))*rtmp
            k2(j) = sign(1.0_rp,k2(j))*rtmp
          endif
          
        endif       ! flip 
      enddo         ! j=1,Np
      
    end subroutine make_periodic_1D
    
  !> Based on 3d wavenumbers, makes the wavenumbers in the directions
  !! kp1 and kp2 periodic, based on the lengths Lp1 and Lp2.
  !! See make_periodic_1D for more details. 
  subroutine make_periodic_2D(k1, kp1, kp2, np, K_total, L1, L2)
    real(kind=rp), intent(inout) :: k1(1), kp1(1), kp2(1)
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
        
        print *, rtmp, n_j1_signed, n_j2_signed

        ! Always reduce the component that is highest
        if (kp1(j) .gt. kp2(j)) then
          n_j1_signed = n_j1_signed - sign(1.0_rp, kp1(j))
          kp1(j) = real(n_j1_signed, kind=rp) * twopi_over_L1     
        else
          n_j2_signed = n_j2_signed - sign(1.0_rp, kp2(j))
          kp2(j) = real(n_j2_signed, kind=rp) * twopi_over_L2  
        end if

        rtmp = K_total_sq - kp1(j)**2.0_rp - kp2(j)**2.0_rp
        print *, rtmp, n_j1_signed, n_j2_signed
        valid_config = (rtmp .gt. 0.0_rp)

      end do

      k1(j) = sign(1.0_rp, k1(j)) * sqrt( rtmp )
     
    enddo
    
  end subroutine make_periodic_2D
  

    subroutine periodicity_chk(kx, ky, kz, np, kk, dlx, dly, dlz, ifxp, ifyp, ifzp)
    real(kind=rp), intent(inout) :: kx(1),ky(1),kz(1)
    integer, intent(in) :: np
    real(kind=rp), intent(in) :: kk
    real(kind=rp), intent(in) :: dlx,dly,dlz
    logical, intent(in) :: ifxp,ifyp,ifzp
    
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
        call make_periodic_1D(ky, kz, kx, np, kk, dlx)
      elseif (ifyp) then
        call make_periodic_1D(kz, kx, ky, np, kk, dly)
      elseif (ifzp) then
        call make_periodic_1D(kx, ky, kz, np, kk, dlz)
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
          
  subroutine gen_dodeca_k(kx,ky,kz,K_tot,Np,seed)
    
    real(kind=rp), intent(inout) :: kx(1),ky(1),kz(1)
    real(kind=rp), intent(in) :: K_tot
    integer, intent(inout) :: Np
    integer, intent(in) :: seed
    
    real(kind=rp) :: pi
    real(kind=rp) :: rotx,roty,rotz
    
    pi = 4.0*atan(1.0)
    
    rotx = ran2(seed)*2.*pi
    roty = ran2(seed)*2.*pi
    rotz = ran2(seed)*2.*pi
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
      100  CONTINUE
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
      100  CONTINUE
      VLAMAX = TMAX
      return
  end function vlamax
  !-----------------------------------------------------------------------
          
end module spec
