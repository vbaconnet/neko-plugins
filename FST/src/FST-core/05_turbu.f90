module turbu
  use num_types, only: rp
  use fst_utils, only : ran2, print_param
  use math, only: pi, abscmp
  use utils, only: neko_error
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
       periodic_z, seed, write_file_path, write_files, gdim, Lx, Ly, Lz)
    
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
    real(kind=rp), intent(inout) :: shell_amp(Nshells)
    logical, intent(in) :: periodic_x, periodic_y, periodic_z
    integer, intent(inout) :: seed
    character(len=*), intent(in) :: write_file_path
    logical, intent(in) :: write_files
    integer, intent(in) :: gdim
    real(kind=rp), intent(in), optional :: Lx, Ly, Lz

    integer :: k,i,j, ierr
    integer :: shellno
    real(kind=rp) :: ue,ve,we
    real(kind=rp) :: uamp,vamp,wamp, u_dot_k, norm_ki
    real(kind=rp) :: amp
    real(kind=rp) :: u_hat(3), u_hat_p(3)
    character(len=LOG_SIZE) :: log_buf

    if (periodic_x .and. abscmp(Lx, 0.0_rp)) &
      call neko_error("Periodic in x requested but total length is zero!")
    if (periodic_y .and. abscmp(Ly, 0.0_rp)) &
      call neko_error("Periodic in y requested but total length is zero!")
    if (periodic_z .and. abscmp(Lz, 0.0_rp)) &
      call neko_error("Periodic in z requested but total length is zero!")

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
       ! NOTE: k_x,k_y,k_z and shell are allocated only on rank 0. This is 
       ! because the size Npeff if computed inside spec_s, which also populates
       ! these arrays! That is why after the end if we allocate the arrays on
       ! all other ranks.
       call spec_s(Npeff, IL, Tu, U_inf, Npmax, Nshells, k_start, k_end, &
         k_x, k_y, k_z, shell, shell_amp, Lx, Ly, Lz, periodic_x, periodic_y, &
         periodic_z, seed, write_file_path, write_files) ! get isotropically distributed wavenumbers in spheres

       ! This will be the total size of our arrays
       n_modes = 2*Npeff*Nshells

    end if

    !
    ! Broadcast total number of modes and allocate arrays
    !
    call MPI_Bcast(n_modes, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

    allocate(phase_shifts(n_modes))
    allocate(random_vectors(n_modes,3))

    !
    ! Populate amplitudes, phase shifts, and project on div.-free psace
    !
    if (pe_rank .eq. 0) then

      block

         ! Temporary arrays for random generation
         real(kind=rp) :: bb(2*Npmax*Nshells, 3), bb1(2*Npmax*Nshells, 3)

         call neko_log%section("Unit vectors & phase shifts")
         call neko_log%message("Generating random phase shifts in [0;2pi]")
         call neko_log%message("Generating random vector components in [-1;1]")
         do k=1, gdim

            ! this loop should be done with Npeff instead BUT we keep it this way
            ! so the ran2 function is called the exact same number of times as the
            ! original code
            do i=1,2*Npmax*Nshells 

               bb(i,k) = ran2(seed)*2.0_rp*pi ! random phase shift

               ! Load phase_shifts, but careful with the bounds
               ! if Npmax < Np, we risk overflow
               if (i .le. n_modes) phase_shifts(i) = bb(i,1)

               bb1(i,k) = 2.0*ran2(seed)-1.0 ! random amplitude

               if (write_files) write(137,*) bb(i,1), bb1(i,1)
            enddo
            
         enddo

       if (write_files) close(137)

       !
       ! Enforce continuity on the random unit vectors
       !              _   _
       !              u . k = 0
       !
       call neko_log%message("Projecting vector components on div.-free space")
       do i = 1, n_modes

          ! u_hat stores the random amplitudes between 0 and 1
          do j = 1,  gdim
             u_hat(j) = bb1(i,j)
          enddo

          u_dot_k = u_hat(1)*k_x(i) + u_hat(2)*k_y(i) + u_hat(3)*k_z(i)
          norm_ki = k_x(i)**2 + k_y(i)**2 + k_z(i)**2
          
          ! Next, project onto divergence-free space
          u_hat_p(1) = u_hat(1) - k_x(i) * u_dot_k / norm_ki
          u_hat_p(2) = u_hat(2) - k_y(i) * u_dot_k / norm_ki
          u_hat_p(3) = u_hat(3) - k_z(i) * u_dot_k / norm_ki

          ! Finally, normalize so the vectors are unitary
          do j=1, gdim
             random_vectors(i,j) = u_hat_p(j) &
                  / sqrt(u_hat_p(1)**2 &
                  + u_hat_p(2)**2 &
                  + u_hat_p(3)**2)
          enddo

       enddo


       !
       ! Write generated modes and amplitudes to file
       !
       if (write_files) then
          call neko_log%message("Writing generated vectors in " // &
                    trim(write_file_path) // '/fst_spectrum.csv')

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

       enddo

       if (write_files) close(13)

      end block

     call neko_log%end_section()

    end if ! end if pe_rank .eq. 0

    !
    ! Broadcast variables so all ranks know what has been generated
    !

    ! Allocate the missing arrays on all ranks != 0 since the allocation
    ! was done on rank 0 in spec_s
    if (pe_rank .ne. 0) then
      call neko_log%message("Allocating arrays on non-zero ranks")
      allocate(k_x(n_modes))
      allocate(k_y(n_modes))
      allocate(k_z(n_modes))
      allocate(shell(n_modes))
      call neko_log%message("Broadcasting generated FST to non-zero ranks")
    end if

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
