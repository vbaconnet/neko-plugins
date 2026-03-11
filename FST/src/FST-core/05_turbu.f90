module turbu
    use neko
    use spec
    use global_params
    implicit none


contains

  !---------------------------------------------------------------------- 
    
  subroutine make_turbu(periodic_x, periodic_y, periodic_z, seed, &
                        write_file_path, dx, dy, dz, gdim, coef)
    logical, intent(in) :: periodic_x, periodic_y, periodic_z
integer, intent(inout) :: seed
    character(len=*), intent(in) :: write_file_path
    real(kind=rp), intent(in), optional :: dx, dy, dz
    integer, intent(in), optional :: gdim 
    type(coef_t), intent(in), optional :: coef

    integer :: k,i,j
    integer :: shellno
    integer :: gdim_ = 3
    !integer :: seed
    real(kind=rp) :: ue,ve,we
    real(kind=rp) :: uamp,vamp,wamp
    real(kind=rp) :: amp, bb1(fst_modes, 3), dlx, dly, dlz
    real(kind=rp) :: u_hat(fst_modes, 3), u_hat_p(fst_modes, 3)
    character(len=LOG_SIZE) :: log_buf

    if (present(gdim)) gdim_ = gdim

    call print_param("nshells", real(nshells, kind=rp))
    call print_param("Npmax", real(Npmax, kind=rp))
    call print_param("fst_ti", fst_ti)
    call print_param("fst_il", fst_il)
    call print_param("kstart", kstart)
    call print_param("kend", kend)

    if (.not. present(coef)) then
      if (.not. (present(dx) .and. present(dy) .and. present(dz))) then
        call neko_error("Please specify either coef or dx,dy,dz")
      end if
    end if

    if (present(dx)) then
      dlx = dx
    else 
    dlx = glmax(coef%dof%x, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv) - &
          glmin(coef%dof%x, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv)
    end if

    if (present(dy)) then
      dly = dy
    else 
    dly = glmax(coef%dof%y, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv) - &
          glmin(coef%dof%y, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv)
    end if

    if (present(dz)) then
      dlz = dz
    else 
    dlz = glmax(coef%dof%z, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv) - &
          glmin(coef%dof%z, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv)
    end if
 
    write (log_buf, *) "[FST] Length in x: ", dlx
    call neko_log%message(log_buf)
    write (log_buf, *) "[FST] Length in y: ", dly
    call neko_log%message(log_buf)
    write (log_buf, *) "[FST] Length in z: ", dlz
    call neko_log%message(log_buf)


    if ( pe_rank.eq.0 ) then
      
      ! Hardcoded seed for random # generation.This way the FST generated is always the same.  
      !seed = -143        

      if (write_files) open(unit=137,form='formatted', &
                            file=trim(write_file_path) // '/bb.txt')

      call spec_s(dlx, dly, dlz, periodic_x, periodic_y, periodic_z, seed, &
              write_file_path) ! get isotropically distributed wavenumbers in spheres

      do k=1,gdim_
        do i=1,fst_modes

                bb(i,k) = ran2(seed)*2.0*pi        ! random phase shift
                bb1(i,k) = 2.0*ran2(seed)-1.0      ! random amplitude

          if (write_files) write(137,*) bb(i,1), bb1(i,1)
          ! write(6,*) 'BB', bb(i,1)
        enddo
      enddo

      if (write_files) close(137)
      ! write(6,*) 'FST - Random amplitude generated'
      call neko_log%message("[FST] Random amplitude generated")
      
      !     make sure that continuity is enforced by ensuring u_vec.k_vec=(0 0 0)
      do i=1,k_length
         do j= 1,gdim_
            u_hat(i,j)=bb1(i,j)
         enddo

         do j=1,gdim_
            u_hat_p(i,j) = u_hat(i,j) &
                 - (u_hat(i,1)*k_num_all(i,1) &
                 + u_hat(i,2)*k_num_all(i,2) &
                 + u_hat(i,3)*k_num_all(i,3)) &
                 * k_num_all(i,j) &
                 / (k_num_all(i,1)**2 &
                 +  k_num_all(i,2)**2 &
                 +  k_num_all(i,3)**2)
         enddo

         do j=1,gdim_
            u_hat_pn(i,j) = u_hat_p(i,j) &
                 / sqrt(u_hat_p(i,1)**2 &
                 + u_hat_p(i,2)**2 &
                 + u_hat_p(i,3)**2)
         enddo
      enddo

      call neko_log%message('[FST] Amplitudes projection done')

      !           Check energy in individual modes
      ue=0.
      ve=0.
      we=0.
      !           Also write the modes
      if (write_files) open(file=trim(write_file_path) // '/fst_spectrum.csv', unit=13)
      if (write_files) write(13,'(7(A, ","),A)') 'ShellNo','kx','ky','kz', &
            'amp','u_hat_pn1','u_hat_pn2', 'u_hat_pn3'
      do i=1,k_length
        shellno = shell(i)
        amp = shell_amp(shellno)

        !write (*,*) "AMP: ", amp

        uamp = u_hat_pn(i,1)*amp
        vamp = u_hat_pn(i,2)*amp
        wamp = u_hat_pn(i,3)*amp
        
        if (write_files) write(13,'(7(g0, ","), g0)') shellno,k_num_all(i,1),k_num_all(i,2), &
        k_num_all(i,3),amp, u_hat_pn(i,1), u_hat_pn(i,2), u_hat_pn(i,3)

        ue =  ue + ((uamp)**2)/2.
        ve =  ve + ((vamp)**2)/2.
        we =  we + ((wamp)**2)/2.
      enddo
      
      if (write_files) close(13)
      
      write(log_buf,'(A18,10x,E12.5E2)') '[FST] Energy in u',ue 
      call neko_log%message(log_buf)      
      write(log_buf,'(A18,10x,E12.5E2)') '[FST] Energy in v',ve 
      call neko_log%message(log_buf)      
      write(log_buf,'(A18,10x,E12.5E2)') '[FST] Energy in w',we 
      call neko_log%message(log_buf)      
      write(log_buf,'(A20,8x,E12.5E2)') '[FST] Estimated tke', &
                                    (ue+ve+we)/2.
      call neko_log%message(log_buf)                                    
      write(log_buf,'(A24,9x,E12.5E2)') '[FST] Estimated Tu*U_inf', &
      sqrt((ue+ve+we)/3.)
      call neko_log%message(log_buf)                                  
    
    end if
                                 
    return
  end subroutine make_turbu
  !----------------------------------------------------------------------       
  
end module turbu
