/*
 Copyright (c) 2025, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas.h>
#include <device_config.h>
#include <device/cuda/check.h>
#include <math/bcknd/device/cuda/mathops_kernel.h>

/**
 * Device kernel for FST
 */
template< typename T, typename U>
__global__ void fst_kernel(
                           const U t,
                           const U Uinf,
                           T * __restrict__ u,
                           T * __restrict__ v,
                           T * __restrict__ w,
                           const int * __restrict__ mask,
                           const int n_mask,
                           const T * __restrict__ ubf,
                           const T * __restrict__ vbf,
                           const T * __restrict__ wbf,
                           const U * __restrict__ k_x,
                           const int n_total_modes,
                           const U * __restrict__ phi_0,
                           const int * __restrict__ shell,
                           const U * __restrict__ shell_amp,
                           const U * __restrict__ randvec,
                           const U cosa,
                           const U sina,
                           const U fringe_time,
                           const U * __restrict__ fs
                           ) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  U ra, rb, rc;
  U phi = 0.0;
  U phi_t = Uinf*t;
  U pert = 0.0;
  U f = 0.0;
  //int shellno;

  // phi_0 is size (k_length, n_mask) in fortran, so that means size
  // [n_mask][k_length] in C and therefore.

  for (int i = idx; i < n_mask; i += str) {

    ra = 0.0;
    rb = 0.0;
    rc = 0.0;
    for (int m = 0; m < n_total_modes; m += 1) {

      phi = phi_0[m + i*n_total_modes] - k_x[m] * phi_t;
      //shellno = shell[m];

      pert = shell_amp[shell[m]-1]*sin(phi);
      //ra += randvec[m + 0*n_total_modes]*pert;
      //rb += randvec[m + 1*n_total_modes]*pert;
      //rc += randvec[m + 2*n_total_modes]*pert;
      ra += randvec[m + 0*n_total_modes]*pert;
      rb += randvec[m + 1*n_total_modes]*pert;
      rc += randvec[m + 2*n_total_modes]*pert;
    }

    f = fringe_time*fs[i];
    u[mask[i+1]-1] = ubf[i] + f*ra;
    v[mask[i+1]-1] = vbf[i] + f*rb;
    w[mask[i+1]-1] = wbf[i] + f*rc;

  }
}

extern "C" {

void cuda_fst(real *t, real *Uinf,
              void *u_d, void *v_d, void *w_d, int *mask_d, int *n_mask,
              void *ubf_d, void *vbf_d, void *wbf_d, void *k_x_d,
              int *n_total_modes, void *phi_0_d, int *shell_d,
              void *shell_amp_d, void *randvec_d, real *cosa, real *sina,
              real *fringe_time, void *fs_d) {

  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(*n_mask, 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

  fst_kernel<real, real>
    <<<nblcks, nthrds, 0, stream>>>(
                                    *t, *Uinf,
                                    (real *) u_d, (real *) v_d, (real *) w_d,
                                    (int *) mask_d, *n_mask,
                                    (real *) ubf_d, (real *) vbf_d, (real *) wbf_d,
                                    (real *) k_x_d, *n_total_modes,
                                    (real *) phi_0_d, (int *) shell_d, (real *) shell_amp_d,
                                    (real *) randvec_d,
                                    *cosa, *sina, *fringe_time, (real *) fs_d
                                    );
  CUDA_CHECK(cudaGetLastError());

}

}
