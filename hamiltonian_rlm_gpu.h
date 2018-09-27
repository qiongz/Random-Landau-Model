#include<math.h>
#include<cuda.h>
#include<cuda_runtime.h>
#include"cublas_v2.h"

__global__ void get_matrix_elements(const cuFloatComplex* coeff_jm, const cuFloatComplex* coeff_m_theta_1,const cuFloatComplex* coeff_m_theta_2, const cuFloatComplex* Vmn, cuFloatComplex* wfs, int dim_n, int off_head);

void set_hamil_matrix(const cuFloatComplex* dev_coeff_m_theta_1,const cuFloatComplex* dev_coeff_m_theta_2, const cuFloatComplex* dev_coeff_jm, const cuFloatComplex *dev_Vmn, cuFloatComplex* dev_wfs, int n_phi,int dim_n,int off_head);

