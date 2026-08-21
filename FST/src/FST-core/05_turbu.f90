module turbu
  use num_types, only: rp
  use fst_utils, only : ran2
  use math, only: glmax, glmin, pi 
  use logger, only: LOG_SIZE, neko_log
  use coefs, only: coef_t
  use spec, only: spec_s
  use mpi_f08, only: MPI_IN_PLACE, MPI_MAX, MPI_MIN, MPI_INTEGER, MPI_Bcast, &
   MPI_Allreduce
  use comm, only: pe_rank, MPI_REAL_PRECISION, NEKO_COMM
  implicit none


contains

  !----------------------------------------------------------------------

  subroutine make_turbu(phase_shifts, random_vectors, Npeff, IL, Tu, U_inf, Npmax, Nshells, k_start, k_end, &
       k_x, k_y, k_z, n_modes, shell, shell_amp, periodic_x, periodic_y, &
       periodic_z, seed, write_file_path, write_files, coef, Lx, Ly, Lz)
    
    real(kind=rp), allocatable :: phase_shifts(:)
    real(kind=rp), allocatable :: random_vectors(:,:)
    integer, intent(out) :: Npeff
    real(kind=rp), intent(in) :: IL, Tu, U_inf
    integer, intent(in) :: Npmax
    integer, intent(in) :: Nshells
    real(kind=rp), intent(in) :: k_start, k_end
    real(kind=rp), allocatable, intent(inout) :: k_x(:), k_y(:), k_z(:)
    integer, intent(inout) :: n_modes ! n_modes = k_length
    integer, allocatable, intent(inout) :: shell(:)
    real(kind=rp), allocatable, intent(inout) :: shell_amp(:)
    logical, intent(in) :: periodic_x, periodic_y, periodic_z
    integer, intent(inout) :: seed
    character(len=*), intent(in) :: write_file_path
    logical, intent(in) :: write_files
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in), optional :: Lx, Ly, Lz

    integer :: k,i,j, ierr
    integer :: shellno
    real(kind=rp) :: ue,ve,we
    real(kind=rp) :: uamp,vamp,wamp, u_dot_k, norm_ki
    real(kind=rp) :: amp, dlx, dly, dlz
    real(kind=rp) :: bb(2*Npmax*Nshells, 3), bb1(2*Npmax*Nshells, 3), &
      u_hat(3), u_hat_p(3)
    character(len=LOG_SIZE) :: log_buf

    if (present(Lx)) then
      dlx = Lx
      write (log_buf, *) "[FST] Length in x: ", dlx
      call neko_log%message(log_buf)
    else
      dlx = glmax(coef%dof%x, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv) - &
          glmin(coef%dof%x, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv)
    end if

    if (present(Ly)) then
      dly = Ly
      write (log_buf, *) "[FST] Length in y: ", dly
      call neko_log%message(log_buf)
    else
      dly = glmax(coef%dof%y, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv) - &
          glmin(coef%dof%y, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv)
    end if

    if (present(Lz)) then
      dlz = Lz
      write (log_buf, *) "[FST] Length in z: ", dlz
      call neko_log%message(log_buf)
    else
      dlz = glmax(coef%dof%z, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv) - &
          glmin(coef%dof%z, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv)
    end if


    !
    ! Generate wavenumbers on rank 0
    !
    if ( pe_rank .eq. 0 ) then

       if (write_files) open(unit=137,form='formatted', &
            file=trim(write_file_path) // '/bb.txt')

       ! 
       ! Generate wavenumbers distributed on spheres
       ! note that Npmax may get modified
       ! 
       print *, "> kstart", k_start, k_end
       call spec_s(Npeff, IL, Tu, U_inf, Npmax, Nshells, k_start, k_end, &
         k_x, k_y, k_z, shell, shell_amp, dlx, dly, dlz, periodic_x, periodic_y, &
         periodic_z, seed, write_file_path, write_files) ! get isotropically distributed wavenumbers in spheres

       n_modes = 2*Npeff*Nshells

    end if

    !
    ! Broadcast total number of modes and allocate arrays
    !
    call MPI_Bcast(n_modes, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

    allocate(phase_shifts(n_modes))
    allocate(shell_amp(n_modes))
    allocate(random_vectors(n_modes,3))

    !
    ! Populate amplitudes, phase shifts, and project on div.-free psace
    !
    if (pe_rank .eq. 0) then

       do k=1,coef%msh%gdim

         ! this loop should be done with Npeff instead BUT we keep it this way
         ! so the ran2 function is called the exact same number of times as the
         ! original code
          do i=1,2*Npmax*Nshells 

             bb(i,k) = ran2(seed)*2.0*pi ! random phase shift
             phase_shifts(i) = bb(i,1)

             bb1(i,k) = 2.0*ran2(seed)-1.0 ! random amplitude

             if (write_files) write(137,*) bb(i,1), bb1(i,1)
          enddo
       enddo

       if (write_files) close(137)
       ! write(6,*) 'FST - Random amplitude generated'
       call neko_log%message("FST - Random amplitude generated")

       !
       ! Enforce continuity on the random unit vectors
       !              _   _
       !              u . k = 0
       !
       do i = 1, n_modes

          ! u_hat stores the random amplitudes between 0 and 1
          do j = 1, coef%msh%gdim
             u_hat(j) = bb1(i,j)
          enddo

          u_dot_k = u_hat(1)*k_x(i) + u_hat(2)*k_y(i) + u_hat(3)*k_z(i)
          norm_ki = k_x(i)**2 + k_y(i)**2 + k_z(i)**2
          
          ! Next, project onto divergence-free space
          u_hat_p(1) = u_hat(1) - k_x(i) * u_dot_k / norm_ki
          u_hat_p(2) = u_hat(2) - k_y(i) * u_dot_k / norm_ki
          u_hat_p(3) = u_hat(3) - k_z(i) * u_dot_k / norm_ki

          ! Finally, normalize so the vectors are unitary
          do j=1,coef%msh%gdim
             random_vectors(i,j) = u_hat_p(j) &
                  / sqrt(u_hat_p(1)**2 &
                  + u_hat_p(2)**2 &
                  + u_hat_p(3)**2)
          enddo

       enddo

       call neko_log%message('FST - Amplitudes projection done')

       !           Check energy in individual modes
       ue=0.
       ve=0.
       we=0.
       !           Also write the modes
       if (write_files) then
          open(file=trim(write_file_path) // '/fst_spectrum.csv', unit=13)
          write(13,'(9(A, ","),A)') 'ShellNo','kx','ky','kz', &
               'u_amp','v_amp','w_amp','u_hat_pn1','u_hat_pn2', 'u_hat_pn3'
       end if
       do i=1, n_modes
          shellno = shell(i)
          amp = shell_amp(shellno)

          uamp = random_vectors(i,1)*amp
          vamp = random_vectors(i,2)*amp
          wamp = random_vectors(i,3)*amp

          if (write_files) write(13,'(9(g0, ","), g0)') shellno, k_x(i), &
            k_y(i), k_z(i), uamp, vamp, wamp, random_vectors(i,1), &
            random_vectors(i,2), random_vectors(i,3)

          ue = ue + ((uamp)**2)/2.
          ve = ve + ((vamp)**2)/2.
          we = we + ((wamp)**2)/2.
       enddo

       if (write_files) close(13)

       write(log_buf,'(A18,10x,E12.5E2)') 'FST - Energy in u',ue
       call neko_log%message(log_buf)
       write(log_buf,'(A18,10x,E12.5E2)') 'FST - Energy in v',ve
       call neko_log%message(log_buf)
       write(log_buf,'(A18,10x,E12.5E2)') 'FST - Energy in w',we
       call neko_log%message(log_buf)
       write(log_buf,'(A20,8x,E12.5E2)') 'FST - Estimated tke', &
            (ue+ve+we)/2.
       call neko_log%message(log_buf)
       write(log_buf,'(A24,9x,E12.5E2)') 'FST - Estimated Tu*U_inf', &
            sqrt((ue+ve+we)/3.)
       call neko_log%message(log_buf)

    end if ! end if pe_rank .eq. 0

    !
    ! Broadcast variables so all ranks know what has been generated
    !
    call MPI_Bcast(k_x, n_modes, &
         MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(k_y, n_modes, &
         MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(k_z, n_modes, &
         MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)

    call MPI_Bcast(random_vectors, n_modes*3, &
         MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)

    call MPI_Bcast(phase_shifts, n_modes, &
         MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)

    call MPI_Bcast(shell, n_modes, &
         MPI_INTEGER , 0, NEKO_COMM, ierr)
    call MPI_Bcast(shell_amp, Nshells , &
         MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)

    return
  end subroutine make_turbu
  !----------------------------------------------------------------------

end module turbu
