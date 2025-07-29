! Copyright (c) 2021-2024, The Neko Authors
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
!> Operators CPU backend
module opr_fst_cpu
  use num_types, only : rp, dp, xp
  implicit none
  private

  public :: opr_fst_cpu_fst

contains

  !> Compute FST on CPU
  subroutine opr_fst_cpu_fst(t, Uinf, u,v,w, mask,n_mask, &
            u_baseflow, v_baseflow, w_baseflow, &
            wavenumbers_x, n_total_modes, phi_0, shell, shell_amplitudes, &
            random_vectors, cosa, sina, fringe_time, fringe_space)

    real(kind=rp), intent(inout), dimension(:,:,:,:) :: u, v, w
    real(kind=rp), intent(in) :: u_baseflow(:), v_baseflow(:), w_baseflow(:)
    integer, intent(in) :: n_mask, n_total_modes
    integer, intent(in) :: mask(0:n_mask), shell(:)
    real(kind=xp), intent(in), dimension(:) :: wavenumbers_x, shell_amplitudes, &
         fringe_space
    real(kind=xp), intent(in), dimension(:,:) :: phi_0, random_vectors
    real(kind=xp), intent(in) :: t, Uinf, fringe_time, sina, cosa

    real(kind=xp) :: phi_t, f, rand_vec(3), pert, urand, vrand, wrand, phi

    integer :: idx, m, shellno

    phi_t = t*Uinf

    ! Loop on all points in the point zone
    do idx = 1, mask(0)

       !> This vector will contain the sum of all Fourier modes
       rand_vec = 0.0_rp

       ! Sum all sin modes for each gll point
       do m = 1, n_total_modes

          ! Random phase shifts
          !phase_shft = bb(m,1)

          ! kx(x - U*t) + ky*y + kz*z + phi
          ! = kx*x + ky*y + kz*z - kx*U*t + phi
          !phi = k_num_all(m,1) * (x(i,1,1,1) - glb_uinf*t) + &
          !     k_num_all(m,2) *  y(i,1,1,1) + &
          !     k_num_all(m,3) *  z(i,1,1,1) + &
          !     phase_shft

          phi = phi_0(m,idx) - wavenumbers_x(m)*phi_t

          shellno = shell(m)

          pert = shell_amplitudes(shellno)*sin(phi)

          rand_vec(1) = rand_vec(1) + random_vectors(m,1)*pert
          rand_vec(2) = rand_vec(2) + random_vectors(m,2)*pert
          rand_vec(3) = rand_vec(3) + random_vectors(m,3)*pert

       enddo

       ! Project the rand_vec if there is a rotation
       urand = rand_vec(1) !rand_vec(1)*cosa - rand_vec(2)*sina
       vrand = rand_vec(2) !rand_vec(1)*sina + rand_vec(2)*cosa
       wrand = rand_vec(3) !rand_vec(3)

       f = fringe_time*fringe_space(idx)
       u(mask(idx),1,1,1) = u_baseflow(idx) + f*urand
       v(mask(idx),1,1,1) = v_baseflow(idx) + f*vrand
       w(mask(idx),1,1,1) = w_baseflow(idx) + f*wrand

    end do

  end subroutine opr_fst_cpu_fst

end module opr_fst_cpu
