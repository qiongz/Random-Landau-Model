#include<iostream>
#include<cuda.h>
#include<cuda_runtime.h>
#include"cublas_v2.h"
#include"magma_v2.h"
#include"magma_operators.h"
using namespace std;

void magma_diag_hamil(magmaFloatComplex* dev_wfs, float*dev_energy,int n_phi);
