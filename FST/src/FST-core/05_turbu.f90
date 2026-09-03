module turbu
  use num_types, only: rp, xp
  use fst_utils, only : ran2, pi
  use utils, only : NEKO_FNAME_LEN
  use math, only: abscmp
  use utils, only: neko_error
  use logger, only: LOG_SIZE, neko_log
  use coefs, only: coef_t
  use spec, only: spec_s
  use mpi_f08, only: MPI_IN_PLACE, MPI_MAX, MPI_MIN, MPI_INTEGER, MPI_Bcast, &
   MPI_Allreduce
  use comm, only: pe_rank, MPI_EXTRA_PRECISION, NEKO_COMM
  implicit none


contains

  !----------------------------------------------------------------------

  subroutine make_turbu(phase_shifts, random_vectors, Npeff, IL, Tu, U_inf, Npmax, Nshells, k_start, k_end, &
       k_x, k_y, k_z, n_modes, shell, shell_amp, periodic_x, periodic_y, &
       periodic_z, seed, write_file_path, write_files, gdim, Lx, Ly, Lz)
    
    real(kind=xp), allocatable :: phase_shifts(:)
    real(kind=xp), allocatable :: random_vectors(:,:)
    integer, intent(out) :: Npeff
    real(kind=xp), intent(in) :: IL, Tu, U_inf
    integer, intent(in) :: Npmax
    integer, intent(in) :: Nshells
    real(kind=xp), intent(in) :: k_start, k_end
    real(kind=xp), allocatable, intent(inout) :: k_x(:), k_y(:), k_z(:)
    integer, intent(inout) :: n_modes ! n_modes = k_length
    integer, allocatable, intent(inout) :: shell(:)
    real(kind=xp), intent(inout) :: shell_amp(Nshells)
    logical, intent(in) :: periodic_x, periodic_y, periodic_z
    integer, intent(inout) :: seed
    character(len=*), intent(in) :: write_file_path
    logical, intent(in) :: write_files
    integer, intent(in) :: gdim
    real(kind=rp), intent(in), optional :: Lx, Ly, Lz

    integer :: k,i,j, ierr
    integer :: shellno
    real(kind=xp) :: ue,ve,we
    real(kind=xp) :: uamp,vamp,wamp, u_dot_k, norm_ki
    real(kind=xp) :: amp
    real(kind=xp) :: u_hat(3), u_hat_p(3)
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: Lx_, Ly_, Lz_

    if (present(Lx)) then
      if (periodic_x .and. abscmp(Lx, 0.0_rp)) &
        call neko_error("Periodic in x requested but total length is zero!")
      Lx_ = Lx
    else
      Lx_ = 1.0_rp
    end if
    if (present(Ly)) then
      if (periodic_y .and. abscmp(Ly, 0.0_rp)) &
        call neko_error("Periodic in y requested but total length is zero!")
      Ly_ = Ly
    else
      Ly_ = 1.0_rp
    end if
    if (present(Lz)) then
      if (periodic_z .and. abscmp(Lz, 0.0_rp)) &
        call neko_error("Periodic in z requested but total length is zero!")
      Lz_ = Lz
    else
      Lz_ = 1.0_rp
    end if

    !
    ! Generate wavenumbers on rank 0
    !
    if ( pe_rank .eq. 0 ) then

       ! 
       ! Generate wavenumbers distributed on spheres
       ! note that Npmax may get modified
       ! 
       ! NOTE: k_x,k_y,k_z and shell are allocated only on rank 0. This is 
       ! because the size Npeff if computed inside spec_s, which also populates
       ! these arrays! That is why after the end if we allocate the arrays on
       ! all other ranks.
       call spec_s(Npeff, IL, Tu, U_inf, Npmax, Nshells, k_start, k_end, &
         k_x, k_y, k_z, shell, shell_amp, Lx_, Ly_, Lz_, periodic_x, periodic_y, &
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
         real(kind=xp) :: bb(2*Npmax*Nshells, 3), bb1(2*Npmax*Nshells, 3)

         call neko_log%section("Unit vectors & phase shifts")
         call neko_log%message("Generating random phase shifts in [0;2pi]")
         call neko_log%message("Generating random vector components in [-1;1]")
         
         do k=1, gdim

            ! this loop should be done with Npeff instead BUT we keep it this way
            ! so the ran2 function is called the exact same number of times as the
            ! original code
            do i=1,2*Npmax*Nshells 

               bb(i,k) = ran2(seed)*2.0_xp*pi ! random phase shift

               ! Load phase_shifts, but careful with the bounds
               ! if Npmax < Np, we risk overflow
               if (i .le. n_modes) phase_shifts(i) = bb(i,1)

               bb1(i,k) = 2.0_xp*ran2(seed)-1.0_xp ! random amplitude
            enddo
            
         enddo

       ! Output random values to file
       call write_bb(write_file_path, bb(:,1), bb1(:,1))

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
       call write_fst_spectrum(write_file_path, shell, k_x, k_y, k_z, &
               shell_amp, random_vectors)
      
      end block ! end bb, bb1

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
         MPI_EXTRA_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(k_y, n_modes, &
         MPI_EXTRA_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(k_z, n_modes, &
         MPI_EXTRA_PRECISION, 0, NEKO_COMM, ierr)

    call MPI_Bcast(random_vectors, n_modes*3, &
         MPI_EXTRA_PRECISION, 0, NEKO_COMM, ierr)

    call MPI_Bcast(phase_shifts, n_modes, &
         MPI_EXTRA_PRECISION, 0, NEKO_COMM, ierr)

    call MPI_Bcast(shell, n_modes, &
         MPI_INTEGER , 0, NEKO_COMM, ierr)

    call MPI_Bcast(shell_amp, Nshells , &
         MPI_EXTRA_PRECISION, 0, NEKO_COMM, ierr)

    return
  end subroutine make_turbu
  !----------------------------------------------------------------------

  subroutine write_bb(path, rand_phase, rand_amp, name)
    character(len=*), intent(in) :: path
    real(kind=xp), intent(in) :: rand_phase(:)
    real(kind=xp), intent(in), optional :: rand_amp(:)
    character(len=*), intent(in), optional :: name

    integer :: unit, ierr, i, n

    character(len=NEKO_FNAME_LEN) :: name_

    name_ = "bb.txt"
    if (present(name)) name_ = trim(name)

    call neko_log%message("Writing " // trim(path) // "/" // trim(name_))
    open(file = trim(path)//'/'//trim(name_), newunit=unit, &
          status="replace", action="write", iostat=ierr)

    if (ierr .ne. 0) call neko_error("Error opening file " // trim(name_))

    n = size(rand_phase)

    if (present(rand_amp)) then
       do i = 1, n
          write(unit, *) rand_phase(i), rand_amp(i)
       end do
    else
       do i = 1, n
          write(unit, *) rand_phase(i), 0.0_xp
       end do
    end if

    close(unit)

  end subroutine write_bb

  subroutine write_fst_spectrum(path, shell, k_x, k_y, k_z, shell_amp, &
                  r_vec, name)
    character(len=*), intent(in) :: path
    integer, intent(in) :: shell(:)
    real(kind=xp), intent(in) :: k_x(:), k_y(:), k_z(:), shell_amp(:)
    real(kind=xp), intent(in) :: r_vec(:,:)
    character(len=*), intent(in), optional :: name

    integer :: unit, ierr, i, n, shellno
    real(kind=xp) :: amp

    character(len=NEKO_FNAME_LEN) :: name_

    name_ = "fst_spectrum.csv"
    if (present(name)) name_ = trim(name)

    call neko_log%message("Writing " // trim(path) // "/" // trim(name_))
    open(file = trim(path)//'/'//trim(name_), newunit=unit, &
          status="replace", action="write", iostat=ierr)

    if (ierr .ne. 0) call neko_error("Error opening file " // trim(name_))

    write(unit,'(7(A, ","),A)') 'ShellNo','kx','ky','kz', &
      'amp','u_hat_pn1','u_hat_pn2', 'u_hat_pn3'

    do i = 1, size(k_x)
     shellno = shell(i)
     amp = shell_amp(shellno)

      write(unit,'(7(g0, ","), g0)') shellno, k_x(i), k_y(i), k_z(i), amp, &
              r_vec(i,1), r_vec(i,2), r_vec(i,3)
    end do

    close(unit)

  end subroutine write_fst_spectrum

end module turbu
