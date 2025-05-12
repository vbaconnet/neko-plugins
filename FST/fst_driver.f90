module fst_driver

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
  !> This is a "point zone" that contains all of the indices of the points in our
  !! forcing zone, where we will force the FST.
  class(point_zone_t), pointer :: FSTBOX
  ! ============================================================================

  ! ============ Sponge =======================================
  !> This is a "point zone" that contains all of the indices of the points in our
  !! sponge zone.
  class(point_zone_t), pointer :: SPONGEBOX
  !> Array containing the forcing fringe for the spong, since it is not
  !! time dependent we store it.
  real(kind=rp), allocatable :: SPONGEBOX_FRINGE(:)
  !> Baseflow in the sponge
  real(kind=rp), allocatable :: SPONGEBOX_BF(:,:)
  logical :: USE_SPONGE
  !======================================================================

  logical :: ENABLED
  logical :: USE_FORCING
  logical :: USE_BASEFLOW

  real(kind=rp) :: UINF
  real(kind=rp) :: VINF
  real(kind=rp) :: WINF

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

  public :: fst_driver_initialize, fst_driver_finalize, fst_driver_apply_BC, &
       fst_driver_apply_forcing

  contains

  !> Initialize user variables or external objects
  subroutine fst_driver_initialize(t, u, v, w, p, coef, params, &
       forcing, baseflow, sponge)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    logical, intent(in) :: forcing
    logical, intent(in) :: baseflow
    logical, intent(in), optional :: sponge

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

    USE_FORCING = forcing
    USE_BASEFLOW = baseflow

    USE_SPONGE = .false.
    if (present(sponge)) USE_SPONGE = sponge

    call json_get(params, "case.timestep", DT)

    !
    ! Read the point zone information from the case file
    !
    if (USE_FORCING) then
       call json_get(params, "case.FST.point_zone", read_str)
       FSTBOX => neko_point_zone_registry%get_point_zone(trim(read_str))

       select type(FSTBOX)
       type is (box_point_zone_t)
          bpz => FSTBOX
       class default
          call neko_error("Wrong point zone type")
       end select
    end if

    !
    ! Read the sponge point zone information
    !
    if (USE_SPONGE) then
       call json_get(params, "case.FST.sponge_point_zone", read_str)
       SPONGEBOX => neko_point_zone_registry%get_point_zone(trim(read_str))

       allocate(SPONGEBOX_FRINGE(SPONGEBOX%size))
       allocate(SPONGEBOX_BF(SPONGEBOX%size, 3))

       select type(SPONGEBOX)
       type is (box_point_zone_t)
          spz => SPONGEBOX
       class default
          call neko_error("Wrong point zone type")
       end select

    end if

    !
    ! Check if we want to output the forcing as an entire field
    !
    call json_get_or_default(params, "case.FST.output_forcing", OUTPUT_FORCING_FIELD, .false.)

    if (OUTPUT_FORCING_FIELD .and. (USE_FORCING .or. USE_SPONGE)) then

       call json_get(params, "case.FST.output_frequency", OUTPUT_FREQ)

       n = 0
       if (USE_FORCING) n = n + FSTBOX%size
       if (USE_BASEFLOW) n = n + SPONGEBOX%size

       ! Artifact to be able to use masked_copy later in forcing,
       ! we store the total mask of the FST region + sponge region
       allocate(STUPID_MASK(0:n))
       STUPID_MASK(0) = n

       if (USE_BASEFLOW .and. USE_FORCING) then
          do i = 1, n
             if ( i .le. FSTBOX%size) then
                STUPID_MASK(i) = FSTBOX%mask(i)
             else
                STUPID_MASK(i) = SPONGEBOX%mask(i - FSTBOX%size)
             end if
          end do
       else if (USE_FORCING) then
          do i = 1, n
             STUPID_MASK(i) = FSTBOX%mask(i)
          end do
       else if (USE_BASEFLOW) then
          do i = 1, n
             STUPID_MASK(i) = SPONGEBOX%mask(i)
          end do
       end if

       call json_get(params, "case.output_directory", read_str)
       call FOUT%init(sp, trim(read_str) // "/forcing", 3)

       call neko_field_registry%add_field(u%dof, 'fu')
       call neko_field_registry%add_field(u%dof, 'fv')
       call neko_field_registry%add_field(u%dof, 'fw')
       fu => neko_field_registry%get_field('fu')
       fv => neko_field_registry%get_field('fv')
       fw => neko_field_registry%get_field('fw')

       call rzero(fu%x, fu%dof%size())
       call rzero(fv%x, fv%dof%size())
       call rzero(fw%x, fw%dof%size())
    else if (OUTPUT_FORCING_FIELD .and. .not. USE_SPONGE) then
       call neko_warning("There is no forcing/sponge, therefore I cannot output&
the forcing :)")
    end if

    call json_get_or_default(params, "case.FST.periodic_x", px, .false.)
    call json_get_or_default(params, "case.FST.periodic_y", py, .false.)
    call json_get_or_default(params, "case.FST.periodic_z", pz, .false.)

    !
    ! Initialize FST stuff
    !
    if (.not. USE_FORCING) then

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
       call FST_OBJ%init_bc(&
            zmin, zmax, delta_z, delta_z, &
            ymin, ymax, delta_y, delta_y, &
            t_start, t_ramp, &
            px, py, pz)

       if (.not. USE_BASEFLOW) then
         call json_get(params, "case.fluid.initial_condition.value(1)", UINF)
         call json_get(params, "case.fluid.initial_condition.value(2)", VINF)
         call json_get(params, "case.fluid.initial_condition.value(3)", WINF)
       end if

    else if (USE_FORCING) then

       !
       ! Read parameters for the FST fringe
       !
       call json_get(params, "case.FST.alpha", alpha)
       delta_x = alpha * (bpz%xmax - bpz%xmin)
       delta_y = alpha * (bpz%ymax - bpz%ymin)
       call json_get(params, "case.FST.fringe_amplitude", amp)
       call json_get(params, "case.FST.t_ramp", t_ramp)
       call json_get_or_default(params, "case.FST.t_start", t_start, 0.0_rp)

       !
       ! Initialize the fst parameters
       !
       call FST_OBJ%init_forcing(&
            bpz%xmin, bpz%xmax, delta_x, delta_x, & ! xstart, xend
            bpz%ymin, bpz%ymax, delta_y, delta_y, & ! xmin, xmax
            amp, &
            t_start, t_ramp, &
            px, py, pz) ! amplitude, delta rise, delta fall

       call FST_obj%generate_forcing(coef, FSTBOX, u, v, w)

    else
       call neko_error("Error initializing the FST")
    end if

    if (USE_SPONGE) then

       ! Build the sponge forcing and store the baseflow inside as well
       ! for the baseflow we use the initial condition right now
       do i = 1, SPONGEBOX%size
          idx = SPONGEBOX%mask(i)
          x = u%dof%x(idx,1,1,1)

          SPONGEBOX_FRINGE(i) = outflow_ramp(x, spz%xmin, spz%xmax, 3000.0_rp)
          SPONGEBOX_BF(i,1) = u%x(idx, 1, 1, 1)
          SPONGEBOX_BF(i,2) = v%x(idx, 1, 1, 1)
          SPONGEBOX_BF(i,3) = w%x(idx, 1, 1, 1)
       end do
    end if

  end subroutine fst_driver_initialize

  subroutine fst_driver_apply_BC(u, v, w, bc, coef, t, tstep, which_solver, angle)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    class(bc_t), target, intent(inout) :: bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    character(len=*), intent(in) :: which_solver
    real(kind=rp), intent(in) :: angle

    integer :: i, idx
    class(bc_t), pointer :: u_bc

    if (.not. ENABLED) return

    u_bc => bc

    if (USE_FORCING) then
       call neko_error("FST was initialized for forcing but is being used as BC")
    end if

    ! Check that we are being called by `fluid`
    if (trim(which_solver) .eq. "fluid") then

         if (.not. allocated(u%x) .or. &
              .not. allocated(v%x) .or. &
              .not. allocated(w%x)) call neko_error("BC not allocated!")

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

    end if

  end subroutine fst_driver_apply_BC

  ! Define the forcing here
  subroutine fst_driver_apply_forcing(f,t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t

    ! Access velocity fields without having them as arguments
    type(field_t), pointer :: u, v, w, fu, fv, fw
    integer :: i, idx

    if (.not. ENABLED) return

    if (.not. USE_SPONGE .and. .not. USE_FORCING) then
       call neko_error("There is no forcing defined for this FST")
    end if

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')

    if (OUTPUT_FORCING_FIELD) then
       fu => neko_field_registry%get_field('fu')
       fv => neko_field_registry%get_field('fv')
       fw => neko_field_registry%get_field('fw')
    end if

    ! Sponge region at the outlet forcing to baseflow
    if (USE_SPONGE) then
       do i = 1, SPONGEBOX%size
          idx = SPONGEBOX%mask(i)
          f%u(idx, 1,1,1) = SPONGEBOX_FRINGE(i)*( SPONGEBOX_BF(i,1) - u%x(idx,1,1,1) )
          f%v(idx, 1,1,1) = SPONGEBOX_FRINGE(i)*( SPONGEBOX_BF(i,2) - v%x(idx,1,1,1) )
          f%w(idx, 1,1,1) = SPONGEBOX_FRINGE(i)*( SPONGEBOX_BF(i,3) - w%x(idx,1,1,1) )
       end do
    end if

    ! Forcing
    if (USE_FORCING) then
       call FST_obj%apply_forcing(FSTBOX, &
            u%dof%x, u%dof%y, u%dof%z, t, &
            u%x, v%x, w%x, &
            f%u, f%v, f%w)
    end if

    ! Output
    if (OUTPUT_FORCING_FIELD .and. mod(COUNTER, OUTPUT_FREQ) .eq. 0) then

       call masked_copy(fu%x, f%u, STUPID_MASK, fu%dof%size(), STUPID_MASK(0))
       call masked_copy(fv%x, f%v, STUPID_MASK, fv%dof%size(), STUPID_MASK(0))
       call masked_copy(fw%x, f%w, STUPID_MASK, fw%dof%size(), STUPID_MASK(0))

       FOUT%fields%items(1)%ptr => fu
       FOUT%fields%items(2)%ptr => fv
       FOUT%fields%items(3)%ptr => fw
       call FOUT%file_%write(FOUT%fields, t - DT)

    end if
    COUNTER = COUNTER + 1

  end subroutine fst_driver_apply_forcing

  ! Finalize user variables or external objects
  subroutine fst_driver_finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    if (.not. ENABLED) return

    call FST_OBJ%free()
    if (allocated(SPONGEBOX_FRINGE)) deallocate(SPONGEBOX_FRINGE)
    if (allocated(SPONGEBOX_BF)) deallocate(SPONGEBOX_BF)

  end subroutine fst_driver_finalize

  ! Linear ramp for the sponge region
  !
  ! ----------------------------> x
  !            ____________
  !           /
  !         /  |          |
  !       /    |          |
  !   __/      |          |
  !     |      |          |
  !   xstart   xramp    xend
  function outflow_ramp(x, xstart, xend, amp) result(frng)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: xstart, xend, amp

    real(kind=rp) :: xramp, frng

    xramp = (xend - xstart)/3.0_rp

    if (x .ge. xstart .and. x .le. xramp) then
       frng = amp * ( x - xstart )/( xramp - xstart )
    else if (x .gt. xramp .and. x .lt. xend) then
       frng = amp
    else
       frng = 0.0_rp
    end if

  end function outflow_ramp
  
end module fst_driver
