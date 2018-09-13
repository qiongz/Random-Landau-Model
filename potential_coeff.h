#include<math.h>
#include<cuda.h>
#include<cuda_runtime.h>
#include"cublas_v2.h"
#define PI2 (M_PI*2.0)

void prepare_potential_coeff(float* impurity_x, float * impurity_y, float *impurity_intensity, cuFloatComplex * Vmn, cuFloatComplex* coeff_mn, int dim_m, int dim_n, int off_head,int impurity_num);
