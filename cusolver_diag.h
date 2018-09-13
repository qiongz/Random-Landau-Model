#include<cuda.h>
#include<cuda_runtime.h>
#include<cusolverDn.h>
#include"cublas_v2.h"

void cusolver_diag_hamil(cuFloatComplex* dev_wfs, float*dev_energy,int n_phi);
