! Copyright (c) 2021-2022, The Neko Authors
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
! "AS IS" AND ANY ERPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
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
!> Operators accelerator backends
module opr_fst_device
  use num_types, only : rp, xp, c_rp, c_xp
  use device, only : device_get_ptr
  use utils, only : neko_error
  use comm
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: opr_fst_device_fst

#ifdef HAVE_CUDA
  interface
     subroutine cuda_fst(t, Uinf, u_d,v_d,w_d, mask_d,n_mask, &
            ubf_d, vbf_d, wbf_d, k_x_d, n_total_modes, phi_0_d, shell_d, &
            shell_amp_d, randvec_d, cosa, sina, fringe_time, fs_d) &
            bind(c, name = 'cuda_fst')
       use, intrinsic :: iso_c_binding
       import c_rp, c_xp
       type(c_ptr), value :: u_d,v_d,w_d,mask_d,ubf_d,vbf_d,wbf_d,&
            k_x_d, phi_0_d, shell_d, shell_amp_d, randvec_d, fs_d
       real(kind=c_rp) :: t
       real(kind=c_xp) :: Uinf, cosa, sina, fringe_time
       integer(c_int) :: n_mask, n_total_modes
     end subroutine cuda_fst
  end interface
#elif HAVE_HIP
  interface
     subroutine hip_fst(t, Uinf, u_d,v_d,w_d, mask_d,n_mask, &
            ubf_d, vbf_d, wbf_d, k_x_d, n_total_modes, phi_0_d, shell_d, &
            shell_amp_d, randvec_d, cosa, sina, fringe_time, fs_d) &
            bind(c, name = 'hip_fst')
       use, intrinsic :: iso_c_binding
       import c_rp, c_xp
       type(c_ptr), value :: u_d,v_d,w_d,mask_d,ubf_d,vbf_d,wbf_d,&
            k_x_d, phi_0_d, shell_d, shell_amp_d, randvec_d, fs_d
       real(kind=c_rp) :: t
       real(kind=c_xp) :: Uinf, cosa, sina, fringe_time
       integer(c_int) :: n_mask, n_total_modes
     end subroutine hip_fst
  end interface
#endif

contains

  subroutine opr_fst_device_fst(t, Uinf,u_d,v_d,w_d, mask_d,n_mask, &
       ubf_d, vbf_d, wbf_d, k_x_d, n_total_modes, phi_0_d, shell_d, &
       shell_amp_d, randvec_d, cosa, sina, fringe_time, fs_d)
    real(kind=c_rp) :: t
    real(kind=c_xp) :: Uinf, cosa, sina, fringe_time
    integer :: n_mask, n_total_modes
    type(c_ptr) :: fs_d, ubf_d, vbf_d, wbf_d, phi_0_d, randvec_d, shell_d, &
         shell_amp_d, k_x_d, u_d, v_d, w_d, mask_d
#ifdef HAVE_CUDA
    call cuda_fst(t, Uinf, u_d,v_d,w_d, mask_d,n_mask, &
            ubf_d, vbf_d, wbf_d, k_x_d, n_total_modes, phi_0_d, shell_d, &
            shell_amp_d, randvec_d, cosa, sina, fringe_time, fs_d)
#elif HAVE_HIP
    call hip_fst(t, Uinf, u_d,v_d,w_d, mask_d,n_mask, &
            ubf_d, vbf_d, wbf_d, k_x_d, n_total_modes, phi_0_d, shell_d, &
            shell_amp_d, randvec_d, cosa, sina, fringe_time, fs_d)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine opr_fst_device_fst

end module opr_fst_device
