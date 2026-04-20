program generate_FST
    use num_types, only: rp
    use turbu, only: make_turbu
    use global_params, only: kstart, kend, Npmax, Nshells
    use utils, only: NEKO_FNAME_LEN, neko_error
    use neko, only: neko_init, neko_finalize
    implicit none

    integer :: seed, argc, stat
    character(len=NEKO_FNAME_LEN) :: inputchar, output_dir
    real(kind=rp) :: dx, dy, dz

    call neko_init

    !
    ! Read arguments
    !
    argc = command_argument_count()

    if ((argc .lt. 2) .or. (argc .gt. 2)) then
        write(*,*) 'Usage: ./generate_FST seed output_dir'
        write(*,*) 'Example command: ./generate_FST -143 my/output/dir'
        stop
    end if
  
    call get_command_argument(1, inputchar)
    read(inputchar, *) seed    
    call get_command_argument(2, output_dir)
    !read(inputchar, *) output_dir
    
    !
    ! Build output directory
    !

    ! First stage: seed
    ! write(output_dir,'(A)') "generated_fst_files/seed_"
    ! if (seed < 0) write (output_dir, '(A,A)') trim(output_dir), 'm'
    ! write (output_dir, '(A,I0)') trim(output_dir), abs(seed)

    ! ! Second stage: parameters
    ! write(output_dir, '(A,A,g0,A,g0,A,I0,A,I0)') trim(output_dir), "_ks_", &
    !                 kstart, "_ke_", kend, "_Ns_", Nshells, "_Np_", Npmax

    write (*,*) "Output directory: ", trim(output_dir)

    ! Create output directory
    call execute_command_line('mkdir -p ' // trim(output_dir), exitstat=stat)
    if (stat /= 0) call neko_error("Failed to create directory" // &
                                    trim(output_dir))

    !
    ! Generate the FST
    !
    dx = 0.79994602501392365_rp
    dy = 0.40000000596046448_rp     
    dz = 0.15000000596046453_rp
    call make_turbu(.false., .false., .true., seed, trim(output_dir), dx=dx, &
                    dy=dy, dz=dz)

    call neko_finalize

end program generate_FST
