#include <math.h>
#include "../diag_wrappers/gpu_runtime.h"

__global__ void get_matrix_elements(const GpuComplex* coeff_jm,
                                    const GpuComplex* coeff_m_theta_1,
                                    const GpuComplex* coeff_m_theta_2,
                                    const GpuComplex* Vmn,
                                    GpuComplex* wfs,
                                    int dim_n,
                                    int off_head);

void set_hamil_matrix(const GpuComplex* dev_coeff_m_theta_1,
                      const GpuComplex* dev_coeff_m_theta_2,
                      const GpuComplex* dev_coeff_jm,
                      const GpuComplex* dev_Vmn,
                      GpuComplex* dev_wfs,
                      int n_phi,
                      int dim_n,
                      int off_head);
