module sponge_driver
  use num_types, only: rp
  use sponge, only: sponge_t
  use field, only: field_t
  use coefs, only: coef_t
  use json_module, only: json_file
  use logger, only: neko_log
  use json_utils, only: json_get, json_get_or_default
  use fluid_user_source_term, only: fluid_user_source_term_t

  implicit none
  type(sponge_t) :: spng

contains

  subroutine sponge_driver_initialize(t,u,v,w,coef,param)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: param

    real(kind=rp) :: alpha_in
    real(kind=rp), allocatable :: amplitudes(:)
    character(len=:), allocatable :: fname, pz_name

    call neko_log%section("Initializing sponge")

    call json_get(param, "case.sponge.alpha_in", alpha_in)
    call json_get(param, "case.sponge.amplitudes", amplitudes)
    call spng%init_from_attributes(alpha_in, amplitudes) 

    call json_get(param, "case.sponge.baseflow_file_name", fname)
    call spng%assign_baseflow_file(trim(fname), .false.)

    call json_get(param, "case.sponge.point_zone_name", pz_name)
    call spng%compute_fringe(trim(pz_name))
    call neko_log%end_section()

  end subroutine sponge_driver_initialize
 
  subroutine sponge_driver_compute(f, t, tmin)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t, tmin

    if (t .lt. tmin) return

    call spng%compute(f, t)
  end subroutine sponge_driver_compute

  subroutine sponge_driver_free()
     call spng%free()
  end subroutine sponge_driver_free

end module sponge_driver
