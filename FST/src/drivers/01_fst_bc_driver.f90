module fst_bc_driver

  use neko
  use FST, only: FST_t
  use fld_file_output, only: fld_file_output_t

  implicit none
  private

  ! ============ FST ===============================================
  !
  ! All of the variable in CAPS are global variables
  !
  !> This is our FST object, contains everything to generate and apply
  !! free stream turbulence
  type(FST_t) :: FST_OBJ
  !> This is a buffer variable to print out things properly
  character(len=LOG_SIZE) :: LOG_BUF
  ! ============================================================================

  logical :: ENABLED
  logical :: FST_GENERATED = .false.


  !
  ! For outputting the forcing as a field file
  !
  type(fld_file_output_t) :: FOUT
  logical :: OUTPUT_FORCING_FIELD
  integer :: OUTPUT_FREQ ! in timesteps
  integer :: COUNTER = 0

  !> This is basically a copy of the forcing zone mask, but with an added
  !! item at STUPID_MASK(0) = size(STUPID_MASK), this is to be able
  !! to use `masked_copy` (see later for clarifications).
  integer, allocatable :: STUPID_MASK(:)

  !> Path to the fst files.
  character(len=:), allocatable :: PATH

  ! ============================================================================

  public :: fst_bc_driver_initialize, fst_bc_driver_finalize, &
       fst_bc_driver_apply

contains

  !> Initialize user variables or external objects
  subroutine fst_bc_driver_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    logical :: px, py, pz
    real(kind=rp) :: x, ymin, ymax, zmin, zmax, delta_y, delta_z, Ly, Lz
    real(kind=rp) :: ystart, yend, zstart, zend
    integer :: i, ierr, n
    real(kind=rp) :: alpha
    alpha = 0.1

    n = coef%dof%size()

    call json_get_or_default(params, "case.FST.enabled", ENABLED, .true.)

    if (.not. ENABLED) then
       call neko_warning("FST is disabled")
       return
    end if

    call json_get_or_default(params, "case.FST.periodic_x", px, .false.)
    call json_get_or_default(params, "case.FST.periodic_y", py, .false.)
    call json_get_or_default(params, "case.FST.periodic_z", pz, .false.)

    if (px) call neko_error("Periodicity in x is not supported yet!")

    ! Read parameters for the FST fringe in space
    if (.not. py .or. .not. pz) call json_get(params, "case.FST.alpha", alpha)

    !
    ! Compute bounds of the domain to set the fringes properly
    ! Note that these values will be overridden if any of these directions
    ! are periodic.
    !
    ymin = glmin(coef%dof%y, n); ymax = glmax(coef%dof%y, n)
    zmin = glmin(coef%dof%z, n); zmax = glmax(coef%dof%z, n)

    call json_get_or_default(params, "case.FST.ystart", ystart, ymin)
    call json_get_or_default(params, "case.FST.yend", yend, ymax)
    delta_y = alpha * (ymax - ymin)

    call json_get_or_default(params, "case.FST.zstart", zstart, zmin)
    call json_get_or_default(params, "case.FST.zend", zend, zmax)
    delta_z = alpha * (zmax - zmin)

    !
    ! Read all FST parameters
    !
    block
      real(kind=rp) :: Uinf, Tu, L, k_start, k_end, t_ramp, t_start
      integer :: Nshells, Npmax, seed

      call json_get_or_default(params, "case.FST.seed", seed, -143)
      if (seed .ge. 0) call neko_error("Seed must be negative!")

      call json_get(params, "case.FST.Uinf", Uinf)
      call json_get(params, "case.FST.Tu", Tu)
      call json_get(params, "case.FST.L", L)
      call json_get(params, "case.FST.k_start", k_start)
      call json_get(params, "case.FST.k_end", k_end)
      call json_get(params, "case.FST.n_shells", Nshells)
      call json_get(params, "case.FST.n_pts_per_shell", Npmax)

      ! Read parameters for the FST fringe in time
      call json_get_or_default(params, "case.FST.t_start", t_start, 0.0_rp)
      call json_get(params, "case.FST.t_ramp", t_ramp)

      call FST_OBJ%init_bc(Uinf, Tu, L, k_start, k_end, Nshells, Npmax, &
         px, py, pz, zmin, zmax, zstart, zend, delta_z, delta_z, &
         ymin, ymax, ystart, yend, delta_y, delta_y, &
         t_start, t_ramp, &
         seed)

    end block

    call json_get_or_default(params, 'case.FST.files_output_path', PATH, &
         "./FST_output_files")
    call system("mkdir -p " // trim(PATH))

  end subroutine fst_bc_driver_initialize

  subroutine fst_bc_driver_apply(u, v, w, bc, coef, t, tstep, angle, on_cpu)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    class(bc_t), intent(in) :: bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(in) :: angle
    logical, intent(in), optional :: on_cpu

    integer :: i, idx

    if (.not. ENABLED) return

    !
    ! At the first timestep we generate the FST based
    ! on the boundry mask!
    !
    if (.not. FST_GENERATED) then
       call FST_obj%generate_bc(coef%dof%x, coef%dof%y, coef%dof%z, &
         bc%msk, bc%msk(0), u, v, w, PATH, coef%msh%gdim)
       FST_GENERATED = .true.
    end if

    ! Then, apply the free stream turbulence that will add on
    ! top of the existing baseflow.
    ! NOTE: if on_cpu is true, memory is not copied to the GPU (you need
    ! to do it yourself)
    call FST_obj%apply_BC(bc%msk, bc%msk(0), &
         u%dof%x, u%dof%y, u%dof%z, t, u%x, v%x, w%x, angle, on_cpu)

    ! If we compute on cpu, copy memory. This is slower!
    if (NEKO_BCKND_DEVICE .eq. 1 .and. on_cpu) then
       if (pe_rank .eq. 0) call neko_warning("You are computing FST on CPU")
       if (pe_rank .eq. 0) call neko_warning("You are copying memory to GPU")
       if (pe_rank .eq. 0) call neko_warning("This is very slow! Check on_host")
       call device_memcpy(u%x, u%x_d, u%size(), HOST_TO_DEVICE, .false.)
       call device_memcpy(v%x, v%x_d, v%size(), HOST_TO_DEVICE, .false.)
       call device_memcpy(w%x, w%x_d, w%size(), HOST_TO_DEVICE, .false.)
    end if

  end subroutine fst_bc_driver_apply

  ! Finalize user variables or external objects
  subroutine fst_bc_driver_finalize()

    if (.not. ENABLED) return

    call FST_OBJ%free()

  end subroutine fst_bc_driver_finalize

end module fst_bc_driver
