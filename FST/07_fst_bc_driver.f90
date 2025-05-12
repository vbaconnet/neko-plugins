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

  real(kind=rp) :: DT
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

    ! Our variables
    type(box_point_zone_t), pointer :: bpz, spz
    type(field_t), pointer :: fu,fv,fw
    character(len=:), allocatable :: read_str, fname
    logical :: px, py, pz
    real(kind=rp) :: x, ymin, ymax, zmin, zmax, amp
    integer :: i, idx, ierr, n
    real(kind=rp) :: alpha, t_ramp, t_start, delta_y, delta_z, delta_x

    call json_get_or_default(params, "case.FST.enabled", ENABLED, .true.)

    if (.not. ENABLED) then
       call neko_warning("FST is disabled")
       return
    end if

    call json_get(params, "case.timestep", DT)
    call json_get_or_default(params, "case.FST.periodic_x", px, .false.)
    call json_get_or_default(params, "case.FST.periodic_y", py, .false.)
    call json_get_or_default(params, "case.FST.periodic_z", pz, .false.)

    ! Compute the bounds of inlet plane for the fringe. Depending
    ! on which one is periodic we will set a fringe or not.
    ymin = 99.0_rp
    ymax = -99_rp
    zmin = -99.0_rp
    zmax = 99.0_rp

    ! Search for min and max
    do i = 1, u%msh%mpts
       ymin = min(ymin, u%msh%points(i)%x(2))
       ymax = max(ymax, u%msh%points(i)%x(2))
       zmin = min(zmin, u%msh%points(i)%x(3))
       zmax = max(zmax, u%msh%points(i)%x(3))
    end do

    call MPI_Allreduce(MPI_IN_PLACE, ymin, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, ymax, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, zmin, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, zmax, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)

    ! Read parameters for the FST fringe
    call json_get(params, "case.FST.alpha", alpha)
    call json_get(params, "case.FST.t_ramp", t_ramp)
    call json_get_or_default(params, "case.FST.t_start", t_start, 0.0_rp)

    ! Set the fringe ramp_in length depending on wether or not we are periodic.
    ! If periodic, then we set a negative to artificially "apply everywhere"
    if (py) then
       delta_y = -99.0_rp * (ymax - ymin)
    else
       delta_y = alpha * (ymax - ymin)
    end if

    if (pz) then
       delta_z = -99.0_rp * (zmax - zmin)
    else
       delta_z = alpha * (zmax - zmin)
    end if

    ! Initialize the fst parameters
    call FST_OBJ%init_bc(zmin, zmax, delta_z, delta_z, ymin, ymax, delta_y, &
         delta_y, t_start, t_ramp, px, py, pz)

  end subroutine fst_bc_driver_initialize

  subroutine fst_bc_driver_apply(u, v, w, bc, coef, t, tstep, angle)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    class(bc_t), target, intent(in) :: bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(in) :: angle

    integer :: i, idx
    class(bc_t), pointer :: u_bc

    if (.not. ENABLED) return

    u_bc => bc

    !
    ! At the first timestep we generate the FST based
    ! on the boundry mask!
    !
    if (tstep .eq. 1) then
       call FST_obj%generate_bc(coef, u_bc%msk, u_bc%msk(0), u=u, v=v, w=w)
    end if

    ! Then, apply the free stream turbulence that will add on
    ! top of the existing baseflow.
    ! NOTE: it does not copy to the GPU
    call FST_obj%apply_BC(u_bc%msk, u_bc%msk(0), &
         u%dof%x, u%dof%y, u%dof%z, t, u%x, v%x, w%x, angle)

  end subroutine fst_bc_driver_apply

  ! Finalize user variables or external objects
  subroutine fst_bc_driver_finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    if (.not. ENABLED) return

    call FST_OBJ%free()

  end subroutine fst_bc_driver_finalize

end module fst_bc_driver
