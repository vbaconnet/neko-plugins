! Copyright (c) 2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!
!> Implements the `sponge_t` type.

module sponge
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use utils, only: extract_fld_file_index, filename_chsuffix
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use global_interpolation, only: global_interpolation_t
  use comm, only: pe_rank
  use json_utils, only : json_get, json_get_or_default
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use fluid_user_source_term, only: fluid_user_source_term_t
  use math, only: sub3, col2, copy, cmult
  use device_math, only: device_sub3, device_col2, device_copy, device_cmult
  use utils, only: neko_error
  use logger, only: neko_log, LOG_SIZE, NEKO_LOG_DEBUG, NEKO_LOG_INFO, NEKO_LOG_VERBOSE
  use field_math, only: field_copy
  use vector, only: vector_t
  use device, only: device_map, device_memcpy, HOST_TO_DEVICE
  use neko_config, only: NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated
  use box_point_zone, only: box_point_zone_t
  use dofmap, only: dofmap_t
  use file, only: file_t
  use field_list, only: field_list_t
  use space, only: space_t, GLL
  use interpolation, only: interpolator_t
  use fld_file_data, only: fld_file_data_t
  implicit none
  private

  type, public :: sponge_t

     !> X velocity component.
     type(field_t), pointer :: u => null()
     !> Y velocity component.
     type(field_t), pointer :: v => null()
     !> Z velocity component.
     type(field_t), pointer :: w => null()

     !> Base flow components
     type(field_t), pointer :: u_bf => null()
     type(field_t), pointer :: v_bf => null()
     type(field_t), pointer :: w_bf => null()
     logical :: baseflow_set = .false.

     !> Fringe array and parameters.
     type(field_t), pointer :: fringe => null()
     real(kind=rp) :: alpha_in
     real(kind=rp) :: amplitudes(3)
     logical :: fringe_set = .false.

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => sponge_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          sponge_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => sponge_free
     !> Assign the base flow field directly from fields
     procedure, pass(this) :: assign_baseflow_field => sponge_assign_baseflow_field
     !> Assign the base flow field from a fld file
     procedure, pass(this) :: assign_baseflow_file => sponge_assign_baseflow_file
     !> Compute the fringe field.
     procedure, pass(this) :: compute_fringe => sponge_compute_fringe
     !> Compute the sponge field.
     procedure, pass(this) :: compute => sponge_compute
  end type sponge_t

contains

  !> Constructor from json.
  subroutine sponge_init_from_json(this, json)
    class(sponge_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    real(kind=rp) :: alpha_in
    real(kind=rp), allocatable :: amplitudes(:)

    call json_get(json, "alpha", alpha_in)
    call json_get(json, "amplitudes", amplitudes)

    call sponge_init_from_attributes(this, alpha_in, amplitudes)
  end subroutine sponge_init_from_json

  !> Actual constructor.
  subroutine sponge_init_from_attributes(this, alpha_in, amplitudes)
    class(sponge_t), intent(inout) :: this
    real(kind=rp), intent(in) :: alpha_in, amplitudes(:)

    call this%free()

    this%alpha_in = alpha_in

    this%amplitudes(1) = amplitudes(1)
    this%amplitudes(2) = amplitudes(2)
    this%amplitudes(3) = amplitudes(3)

    if (pe_rank .eq. 0) write (*,*) "(SPONGE) Amplitudes", this%amplitudes

    call neko_log%message("Initializing sponge", lvl = NEKO_LOG_INFO)

    call neko_log%message("(SPONGE) Pointing at fields u,v,w", lvl = NEKO_LOG_INFO)
    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")

    call neko_log%message("(SPONGE) Initializing bf and fringe fields", lvl = NEKO_LOG_INFO)
    call neko_field_registry%add_field(this%u%dof, "sponge_u_bf")
    call neko_field_registry%add_field(this%u%dof, "sponge_v_bf")
    call neko_field_registry%add_field(this%u%dof, "sponge_w_bf")
    call neko_field_registry%add_field(this%u%dof, "sponge_fringe") 
    this%u_bf => neko_field_registry%get_field("sponge_u_bf")
    this%v_bf => neko_field_registry%get_field("sponge_v_bf")
    this%w_bf => neko_field_registry%get_field("sponge_w_bf")
    this%fringe => neko_field_registry%get_field("sponge_fringe")
    
    call neko_log%message("Done initializing sponge", lvl = NEKO_LOG_INFO)

  end subroutine sponge_init_from_attributes

  !> Copy the baseflow from a set of fields
  !! param @fname e.g. "field0.f00010"
  subroutine sponge_assign_baseflow_file(this, raw_fname, interpolate)
    class(sponge_t), intent(inout) :: this
    character(len=*), intent(in) :: raw_fname
    logical, intent(in) :: interpolate

    character(len=1024) :: fname
    type(global_interpolation_t) :: global_interp
    type(file_t) :: f
    type(fld_file_data_t) :: fld
    integer :: sample_idx
    type(space_t) :: prev_Xh
    type(interpolator_t) :: space_interp

    call neko_log%message("(SPONGE) Assigning baseflow from file", lvl = NEKO_LOG_INFO)
    ! Extract idx from ".f00010"
    call neko_log%message("(SPONGE) reading file " // trim(raw_fname), lvl = NEKO_LOG_INFO)
    sample_idx = extract_fld_file_index(trim(raw_fname), -1)
    if (sample_idx .eq. -1) &
         call neko_error("(SPONGE) Invalid file name. The&
      & file format must be e.g. 'mean0.f00001'")

    call filename_chsuffix(raw_fname, fname, 'fld')

    ! Read field file
    call f%init(trim(fname))
    call fld%init
    call f%set_counter(sample_idx)
    call f%read(fld)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_log%message("(SPONGE) fld data copying to device" , lvl = NEKO_LOG_INFO)
       call device_memcpy(fld%x%x, fld%x%x_d, fld%x%n, HOST_TO_DEVICE, .false.)
       call device_memcpy(fld%y%x, fld%y%x_d, fld%y%n, HOST_TO_DEVICE, .false.)
       call device_memcpy(fld%z%x, fld%z%x_d, fld%z%n, HOST_TO_DEVICE, .true.)
       call device_memcpy(fld%u%x, fld%u%x_d, fld%u%n, HOST_TO_DEVICE, .false.)
       call device_memcpy(fld%v%x, fld%v%x_d, fld%v%n, HOST_TO_DEVICE, .false.)
       call device_memcpy(fld%p%x, fld%p%x_d, fld%p%n, HOST_TO_DEVICE, .true.)
       call device_memcpy(fld%w%x, fld%w%x_d, fld%w%n, HOST_TO_DEVICE, .true.)
    end if

    if (interpolate) then

       ! Generates an interpolator object and performs the point search
       call neko_log%message("(SPONGE) Generating global interpolator " // trim(fname), lvl = NEKO_LOG_INFO)
       call fld%generate_interpolator(global_interp, this%u%dof, this%u%msh, &
            0.000001_rp)

       ! Evaluate velocities and pressure
       call neko_log%message("(SPONGE) Interpolating " // trim(fname), lvl = NEKO_LOG_INFO)
       call global_interp%evaluate(this%u_bf%x, fld%u%x , .false. )
       call global_interp%evaluate(this%v_bf%x, fld%v%x , .false. )
       call global_interp%evaluate(this%w_bf%x, fld%w%x , .false. )

       call global_interp%free

    else

       ! Build a space_t object from the data in the fld file
       call prev_Xh%init(GLL, fld%lx, fld%ly, fld%lz)
       call space_interp%init(this%u%Xh, prev_Xh)

       ! Do the space-to-space interpolation
       call space_interp%map(this%u_bf%x, fld%u%x, fld%nelv, this%u%Xh)
       call space_interp%map(this%v_bf%x, fld%v%x, fld%nelv, this%v%Xh)
       call space_interp%map(this%w_bf%x, fld%w%x, fld%nelv, this%w%Xh)

       call space_interp%free
 
       !if (NEKO_BCKND_DEVICE .eq. 1) then
       !   call neko_log%message("(SPONGE) memcopying data on gpu", lvl = NEKO_LOG_INFO)
       !   call device_memcpy(this%u_bf%x, this%u_bf%x_d, this%u_bf%size(), HOST_TO_DEVICE, .false.)
       !   call device_memcpy(this%v_bf%x, this%v_bf%x_d, this%u_bf%size(), HOST_TO_DEVICE, .false.)
       !   call device_memcpy(this%w_bf%x, this%w_bf%x_d, this%u_bf%size(), HOST_TO_DEVICE, .true.)
       !end if

    end if
    
    call neko_log%message("(SPONGE) done assigning baseflow file", lvl = NEKO_LOG_INFO)

    this%baseflow_set = .true.

  end subroutine sponge_assign_baseflow_file

  !> Copy the baseflow from a set of fields
  subroutine sponge_assign_baseflow_field(this, u,v,w)
    class(sponge_t), intent(inout) :: this
    type(field_t), intent(in) :: u,v,w

    type(field_list_t) :: fld_list
    type(file_t) :: f
    !type(field_t), pointer :: us,vs,ws

    call neko_log%message("(SPONGE) Assigning baseflow", lvl = NEKO_LOG_INFO)
    call field_copy(this%u_bf, u)
    call field_copy(this%v_bf, v)
    call field_copy(this%w_bf, w)
    call neko_log%message("(SPONGE) done assigning baseflow", lvl = NEKO_LOG_INFO)

    this%baseflow_set = .true.

    !call neko_log%message("Output BF to SPNG_BF.fld")
    !f = file_t("SPNG_BF.fld")
    !call fld_list%init(3)
    !
    !us => u 
    !vs => v
    !ws => w
    !call fld_list%assign(1, us)
    !call fld_list%assign(2, vs)
    !call fld_list%assign(3, ws)
    !call f%write(fld_list)
    !call fld_list%free()
    !call neko_log%message("... Done")
    !nullify(us)
    !nullify(vs)
    !nullify(ws)

  end subroutine sponge_assign_baseflow_field

  !> Build the fringe based on the type of point zone
  ! Linear ramp for the sponge region. xramp is controlled
  ! by alpha_in which is the percentage of the sponge region
  ! that should be linear
  !
  ! ----------------------------> x
  !            ____________
  !           /
  !         /  |          |
  !       /    |          |
  !   __/      |          |
  !     |      |          |
  !   xstart   xramp    xend
  ! Currently only supports box point zones
  subroutine sponge_compute_fringe(this, zone_name)
    class(sponge_t), intent(inout) :: this
    character(len=*), intent(in) :: zone_name

    !type(box_point_zone_t) :: b
    class(point_zone_t), pointer :: zone
    integer :: i, n, imask
    real(kind=rp) :: x, dx

    call neko_log%message("(SPONGE) Computing fringe", lvl = NEKO_LOG_INFO)

    zone => neko_point_zone_registry%get_point_zone(trim(zone_name))
    call neko_log%message("(SPONGE) Fringing zone " // trim(zone%name), lvl = NEKO_LOG_INFO)

       select type(b => zone)
       type is (box_point_zone_t)

          n = b%size
          dx = this%alpha_in*(b%xmax - b%xmin)

          call neko_log%message("(SPONGE) doing ramp calculation", lvl = NEKO_LOG_INFO)
          do i = 1, b%size
             imask = b%mask(i)
             x = this%u%dof%x(imask, 1, 1, 1)

             this%fringe%x(imask,1,1,1) = S((x - b%xmin)/dx)
          
          end do
       call neko_log%message("(SPONGE) done fringing zone " // trim(zone%name), lvl = NEKO_LOG_INFO)
       class default
          call neko_error("Wrong point zone type")
       end select
    call neko_log%message("(SPONGE) done fringing zone " // trim(zone%name), lvl = NEKO_LOG_INFO)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_log%message("(SPONGE) fringe, copying to device" , lvl = NEKO_LOG_INFO)
       call device_memcpy(this%fringe%x, this%fringe%x_d, &
            this%fringe%size(), HOST_TO_DEVICE, .false.)
    end if

    this%fringe_set = .true.

  end subroutine sponge_compute_fringe

  !> Destructor.
  subroutine sponge_free(this)
    class(sponge_t), intent(inout) :: this

    nullify(this%u_bf)
    nullify(this%v_bf)
    nullify(this%w_bf)
    nullify(this%fringe)
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)

    this%fringe_set = .false.
    this%baseflow_set = .false.

  end subroutine sponge_free

  !> Compute the sponge field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine sponge_compute(this, f, t)
    class(sponge_t), intent(inout) :: this
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t

    if (.not. this%baseflow_set) call neko_error("SPONGE: No baseflow set")
    if (.not. this%fringe_set) call neko_error("SPONGE: No fringe set")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_log%message("(SPONGE) Computing sponge device" , lvl = NEKO_LOG_INFO)
       call device_sub3(f%u_d, this%u_bf%x_d, this%u%x_d, this%u%size())
       call device_sub3(f%v_d, this%v_bf%x_d, this%v%x_d, this%v%size())
       call device_sub3(f%w_d, this%w_bf%x_d, this%w%x_d, this%w%size())

       call device_col2(f%u_d, this%fringe%x_d, this%fringe%size())
       call device_cmult(f%u_d, this%amplitudes(1), this%fringe%size())

       call device_col2(f%v_d, this%fringe%x_d, this%fringe%size())
       call device_cmult(f%v_d, this%amplitudes(2), this%fringe%size())
       
       call device_col2(f%w_d, this%fringe%x_d, this%fringe%size())
       call device_cmult(f%w_d, this%amplitudes(3), this%fringe%size())
    else
       call neko_log%message("(SPONGE) Computing sponge host" , lvl = NEKO_LOG_INFO)
       call sub3(f%u, this%u_bf%x, this%u%x, this%u%size())
       call sub3(f%v, this%v_bf%x, this%v%x, this%v%size())
       call sub3(f%w, this%w_bf%x, this%w%x, this%w%size())

       call col2(f%u, this%fringe%x, this%fringe%size())
       call cmult(f%u, this%amplitudes(1), this%fringe%size())
       call col2(f%v, this%fringe%x, this%fringe%size())
       call cmult(f%v, this%amplitudes(2), this%fringe%size())
       call col2(f%w, this%fringe%x, this%fringe%size())
       call cmult(f%w, this%amplitudes(3), this%fringe%size())
    end if

  end subroutine sponge_compute

  ! Smooth step function, 0 if x <= 0, 1 if x >= 1, 1/exp(1/(x-1) + 1/x) between 0 and 1
  function S(x) result(y)
    real(kind=rp), intent(in) :: x
    real(kind=rp)             :: y

    if ( x.le.0._rp ) then
       y = 0._rp
    else if ( x.ge.1._rp ) then
       y = 1._rp
    else
       y = 1._rp / (1._rp + exp( 1._rp/(x-1._rp) + 1._rp/x))
    end if

  end function S
end module sponge
