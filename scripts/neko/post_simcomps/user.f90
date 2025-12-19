module user
  use neko
  implicit none

  type(fld_file_data_t) :: fld
  type(file_t) :: f

contains
 
  ! Register user defined functions here (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u

    u%user_finalize_modules => finalize 
    u%user_init_modules => initialize
    u%user_check => usercheck
  end subroutine user_setup

  ! Sync a vector to device
  subroutine sync_vec(vec, direction, sync)
    type(vector_t), intent(inout) :: vec
    integer, intent(in), optional :: direction
    logical, intent(in), optional :: sync

    integer :: dir_ = host_to_device
    logical :: sync_ = .false.

    if (present(direction)) dir_ = direction
    if (present(sync)) sync_ = sync

    !call vec%copyto(dir_, sync_)
    if (NEKO_BCKND_DEVICE .eq. 1) call device_memcpy(vec%x, vec%x_d, vec%size(), dir_, sync_)
  end subroutine sync_vec
  
  ! Sync a fld to device 
  subroutine sync_fld(fld, sync_mesh, direction, sync)
    type(fld_file_data_t), intent(inout) :: fld
    logical, intent(in) :: sync_mesh
    integer, intent(in), optional :: direction
    logical, intent(in), optional :: sync

    call sync_vec(fld%u, direction, sync)
    call sync_vec(fld%v, direction, sync)
    call sync_vec(fld%w, direction, sync)
    call sync_vec(fld%p, direction, .true.)
    
    if (sync_mesh) then
      call sync_vec(fld%x, direction, sync)
      call sync_vec(fld%y, direction, sync)
      call sync_vec(fld%z, direction, .true.)
    end if
    
  end subroutine sync_fld

  ! Initialize user variables or external objects
  subroutine initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    integer :: i, i_c

    character(len=:), allocatable :: read_str

    call json_get(params, "case.post.file_series", read_str)

    call f%init(trim(read_str))
    call f%read(fld) ! read fld file (on CPU)
    call sync_fld(fld, .true.) ! sync all fields to device (GPU)
 
    ! Copy contents of the fld file into  
    call device_copy(u%x_d, fld%u%x_d, fld%u%size())
    call device_copy(v%x_d, fld%v%x_d, fld%v%size())
    call device_copy(w%x_d, fld%w%x_d, fld%w%size())

  end subroutine initialize
  
  ! This is called at the end of every time step
  subroutine usercheck(t, tstep, u, v, w, p, coef, param)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: param
   
    call f%read(fld)
    call sync_fld(fld, .false.)    
    
    call device_copy(u%x_d, fld%u%x_d, fld%u%size())
    call device_copy(v%x_d, fld%v%x_d, fld%v%size())
    call device_copy(w%x_d, fld%w%x_d, fld%w%size())

  end subroutine usercheck


 ! user-defined boundary condition
  subroutine user_bc(field_bc_list, bc_bc_list, coef, t, tstep, which_solver)
   type(field_list_t), intent(inout) :: field_bc_list
   type(bc_list_t), intent(inout) :: bc_bc_list
  type(coef_t), intent(inout) :: coef
  real(kind=rp), intent(in) :: t
  integer, intent(in) :: tstep
  character(len=*), intent(in) :: which_solver


 end subroutine user_bc

  ! Finalize user variables or external objects
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    call fld%free()

  end subroutine finalize
end module user
