! Copyright (c) 2020-2024, The Neko Authors
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
!> Operators
module fst_operator
  use neko_config, only : NEKO_BCKND_SX, NEKO_BCKND_DEVICE, NEKO_BCKND_XSMM, &
       NEKO_DEVICE_MPI
  use num_types, only : rp, xp
  use opr_fst_cpu, only: opr_fst_cpu_fst
  use opr_fst_device, only: opr_fst_device_fst
  use device, only: device_get_ptr
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  public :: fst_bc_compute

contains

  subroutine fst_bc_compute(t, Uinf, u, v, w, mask, n_mask, &
       u_baseflow, v_baseflow, w_baseflow, &
       wavenumbers_x, n_total_modes, phi_0, shell, shell_amplitudes, &
       random_vectors, angleXY, fringe_time, fringe_space, on_host)
    real(kind=rp), intent(inout), dimension(:,:,:,:) :: u, v, w
    real(kind=rp), intent(in) :: u_baseflow(:), v_baseflow(:), w_baseflow(:)
    integer, intent(in) :: n_mask, n_total_modes
    integer, intent(in) :: mask(0:n_mask), shell(:)
    real(kind=rp), intent(in), dimension(:) :: wavenumbers_x, shell_amplitudes, &
         fringe_space
    real(kind=rp), intent(in), dimension(:,:) :: phi_0, random_vectors
    real(kind=rp), intent(in) :: t, Uinf, fringe_time, angleXY
    logical, intent(in) :: on_host

    type(c_ptr) :: fs_d, ubf_d, vbf_d, wbf_d, phi_0_d, randvec_d, shell_d, &
         shell_amp_d, k_x_d, u_d, v_d, w_d, mask_d

    integer :: i
    real(kind=rp) :: sina, cosa
    cosa = cos(angleXY)
    sina = sin(angleXY)

    if (NEKO_BCKND_DEVICE .eq. 1 .and. .not. on_host) then
       u_d = device_get_ptr(u)
       v_d = device_get_ptr(v)
       w_d = device_get_ptr(w)
       fs_d = device_get_ptr(fringe_space)
       ubf_d = device_get_ptr(u_baseflow)
       vbf_d = device_get_ptr(v_baseflow)
       wbf_d = device_get_ptr(w_baseflow)
       phi_0_d = device_get_ptr(phi_0)
       randvec_d = device_get_ptr(random_vectors)
       shell_d = device_get_ptr(shell)
       shell_amp_d = device_get_ptr(shell_amplitudes)
       k_x_d = device_get_ptr(wavenumbers_x)
       mask_d = device_get_ptr(mask)

       call opr_fst_device_fst(t, Uinf, u_d,v_d,w_d, mask_d,n_mask, &
            ubf_d, vbf_d, wbf_d, k_x_d, n_total_modes, phi_0_d, shell_d, &
            shell_amp_d, randvec_d, cosa, sina, fringe_time, fs_d)
    else
       call opr_fst_cpu_fst(t, Uinf, u,v,w, mask,n_mask, &
            u_baseflow, v_baseflow, w_baseflow, &
            wavenumbers_x, n_total_modes, phi_0, shell, shell_amplitudes, &
            random_vectors, cosa, sina, fringe_time, fringe_space)
    end if

  end subroutine fst_bc_compute
end module fst_operator
