#include "hamiltonian_rlm_gpu.h"

__global__ void get_matrix_elements(const GpuComplex* coeff_jm,
                                    const GpuComplex* coeff_m_theta_1,
                                    const GpuComplex* coeff_m_theta_2,
                                    const GpuComplex* Vmn,
                                    GpuComplex* wfs,
                                    int dim_n,
                                    int off_head) {
    int nphi = gridDim.x * blockDim.x;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    GpuComplex sum;
    int wfs_idx = j * nphi + k;
    int idx_m, n, m;
    wfs[wfs_idx].x = wfs[wfs_idx].y = 0;
    if (k >= j) {
        for (m = 0; m < dim_n; m++) {
            sum.x = sum.y = 0;
            idx_m = m * dim_n;
            if (k - j <= off_head) {
                n = off_head + k - j;
                sum = GPU_CX_ADD(sum, GPU_CX_MUL(Vmn[idx_m + n], GPU_CX_CONJ(coeff_m_theta_1[n])));
            }
            if ((nphi + j - k) <= off_head) {
                n = off_head + k - j - nphi;
                sum = GPU_CX_ADD(sum, GPU_CX_MUL(Vmn[idx_m + n], GPU_CX_CONJ(coeff_m_theta_1[n])));
            }
            GpuComplex term = GPU_CX_MUL(GPU_CX_MUL(sum, coeff_jm[j * dim_n + m]), coeff_m_theta_2[m]);
            wfs[wfs_idx] = GPU_CX_ADD(wfs[wfs_idx], term);
        }
        wfs[k * nphi + j] = GPU_CX_CONJ(wfs[wfs_idx]);
    }
}

void set_hamil_matrix(const GpuComplex* dev_coeff_m_theta_1,
                      const GpuComplex* dev_coeff_m_theta_2,
                      const GpuComplex* dev_coeff_jm,
                      const GpuComplex* dev_Vmn,
                      GpuComplex* dev_wfs,
                      int n_phi,
                      int dim_n,
                      int off_head) {
    dim3 grids(n_phi / 8, n_phi / 8);
    dim3 blocks(8, 8);
    get_matrix_elements<<<grids, blocks, 0, NULL>>>(
        dev_coeff_jm, dev_coeff_m_theta_1, dev_coeff_m_theta_2, dev_Vmn, dev_wfs, dim_n, off_head);
    gpuDeviceSynchronize();
}
