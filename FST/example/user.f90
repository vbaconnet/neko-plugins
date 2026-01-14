module user
  use neko
  use fst_bc_driver
  implicit none

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%initialize => initialize ! Initialize parameters
    u%finalize => finalize   ! Finalize
    u%dirichlet_conditions => user_bc    ! Update boundary conditions
  end subroutine user_setup

  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    type(field_t)  , pointer :: u
    type(field_t)  , pointer :: v
    type(field_t)  , pointer :: w
    type(field_t)  , pointer :: p
    type(coef_t)   , pointer :: coef

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")
    p => neko_registry%get_field("p")

    coef => neko_user_access%case%fluid%c_Xh
   
    call fst_bc_driver_initialize(time%t,u,v,w,p,coef,neko_user_access%case%params) 

  end subroutine initialize

  !> Set boundary conditions
  subroutine user_bc(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    
    type(field_t), pointer :: ui,vi,wi,pi
    type(coef_t)   , pointer :: coef

    ui => null()
    vi => null()
    wi => null()
    pi => null()
    coef => null()

    !
    ! Velocity BC
    !
    if (trim(fields%items(1)%ptr%name) .eq. "u") then

       associate(u => fields%items(1)%ptr, &
            v => fields%items(2)%ptr, &
            w => fields%items(3)%ptr)

         !
         ! At the first time step, apply velocity BC by copying
         ! what has been applied in the initial condition
         !
         if (time%tstep .eq. 1) then    
          
            ! Get solution fields (u,v,w) which contain the initial condition
            ! if we are at the first timestep 
            ui => neko_registry%get_field("u")
            vi => neko_registry%get_field("v")
            wi => neko_registry%get_field("w")
            ! Pressure can be handled directly on GPU
            if (neko_bcknd_device .eq. 1) then
               call device_masked_copy_0(u%x_d, ui%x_d, bc%msk_d, u%dof%size(), bc%msk(0))
               call device_masked_copy_0(v%x_d, vi%x_d, bc%msk_d, v%dof%size(), bc%msk(0))
               call device_masked_copy_0(w%x_d, wi%x_d, bc%msk_d, w%dof%size(), bc%msk(0))
            else
               call masked_copy_0(u%x, ui%x, bc%msk, u%dof%size(), bc%msk(0))
               call masked_copy_0(v%x, vi%x, bc%msk, v%dof%size(), bc%msk(0))
               call masked_copy_0(w%x, wi%x, bc%msk, w%dof%size(), bc%msk(0))
            end if
         
         end if


         !
         ! Apply FST
         !
         coef => neko_user_access%case%fluid%c_Xh
         call fst_bc_driver_apply(u, v, w, bc, coef, time%t, time%tstep, 0.0_rp, .false.)
       end associate

    !
    ! Pressure BC
    !
    else if (trim(fields%items(1)%ptr%name) .eq. "p") then

       associate(p => fields%items(1)%ptr)
         
         if (time%tstep .eq. 1) then    
            pi => neko_registry%get_field("p")
         
            ! Pressure can be handled directly on GPU
            if (neko_bcknd_device .eq. 1) then
               call device_masked_copy_0(p%x_d, pi%x_d, bc%msk_d, p%dof%size(), bc%msk(0))
            else
               call masked_copy_0(p%x, pi%x, bc%msk, p%dof%size(), bc%msk(0))
            end if

         end if

       end associate
    end if

  end subroutine user_bc

  ! Free relevant objects 
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    call fst_bc_driver_finalize()

  end subroutine finalize

  end module user
