module user
  use neko
  use fst_bc_driver
  implicit none

  character(len=LOG_SIZE) :: log_buf
contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules     => initialize ! Initialize parameters
    u%user_finalize_modules => finalize   ! Finalize
    u%user_dirichlet_update => user_bc    ! Update boundary conditions
  end subroutine user_setup

  subroutine initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    call fst_bc_driver_initialize(t,u,v,w,p,coef,params) 

  end subroutine initialize

  !> Set boundary conditions
  subroutine user_bc(field_bc_list, bc, coef, t, tstep)
    type(field_list_t), intent(inout) :: field_bc_list
    type(field_dirichlet_t), intent(in) :: bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
   
    type(field_t), pointer :: ui,vi,wi,pi

    !
    ! Velocity BC
    !
    if (trim(field_bc_list%items(1)%ptr%name) .eq. "u") then

       associate(u => field_bc_list%items(1)%ptr, &
            v => field_bc_list%items(2)%ptr, &
            w => field_bc_list%items(3)%ptr)

         !
         ! At the first time step, apply velocity BC by copying
         ! what has been applied in the initial condition
         !
         if (tstep .eq. 1) then    
          
            ! Get solution fields (u,v,w) which contain the initial condition
            ! if we are at the first timestep 
            ui => neko_field_registry%get_field("u")
            vi => neko_field_registry%get_field("v")
            wi => neko_field_registry%get_field("w")
            ! Pressure can be handled directly on GPU
            if (neko_bcknd_device .eq. 1) then
               call device_masked_copy(u%x_d, ui%x_d, bc%msk_d, u%dof%size(), bc%msk(0))
               call device_masked_copy(v%x_d, vi%x_d, bc%msk_d, v%dof%size(), bc%msk(0))
               call device_masked_copy(w%x_d, wi%x_d, bc%msk_d, w%dof%size(), bc%msk(0))
            else
               call masked_copy(u%x, ui%x, bc%msk, u%dof%size(), bc%msk(0))
               call masked_copy(v%x, vi%x, bc%msk, v%dof%size(), bc%msk(0))
               call masked_copy(w%x, wi%x, bc%msk, w%dof%size(), bc%msk(0))
            end if

            ! We need to do this on the CPU since the FST is applied 
            ! on the CPU on top of the baseflow given by u%x, v%x, w%x.
!!$            call masked_copy(u%x, ui%x, bc%msk, u%dof%size(), bc%msk(0))
!!$            call masked_copy(v%x, vi%x, bc%msk, v%dof%size(), bc%msk(0))
!!$            call masked_copy(w%x, wi%x, bc%msk, w%dof%size(), bc%msk(0))
         end if

         !
         ! Apply FST
         !
         call fst_bc_driver_apply(u,v,w,bc,coef,t,tstep,angle=0.0_rp, memcpy=.false.)

         ! Synchronize to GPU
!!$         if (neko_bcknd_device .eq. 1) then
!!$            call device_memcpy(u%x, u%x_d, u%dof%size(), &
!!$                    host_to_device, sync = .false.)
!!$            call device_memcpy(v%x, v%x_d, v%dof%size(), &
!!$                    host_to_device, sync = .false.)
!!$            call device_memcpy(w%x, w%x_d, w%dof%size(), &
!!$                    host_to_device, sync = .false.)
!!$
!!$         end if
       end associate

    !
    ! Pressure BC
    !
    else if (trim(field_bc_list%items(1)%ptr%name) .eq. "p") then

       associate(p => field_bc_list%items(1)%ptr)
         
         if (tstep .eq. 1) then    
            pi => neko_field_registry%get_field("p")
         
            ! Pressure can be handled directly on GPU
            if (neko_bcknd_device .eq. 1) then
               call device_masked_copy(p%x_d, pi%x_d, bc%msk_d, p%dof%size(), bc%msk(0))
            else
               call masked_copy(p%x, pi%x, bc%msk, p%dof%size(), bc%msk(0))
            end if

         end if

       end associate
    end if

  end subroutine user_bc

  ! Free relevant objects 
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    call fst_bc_driver_finalize(t, params)

  end subroutine finalize

  end module user
