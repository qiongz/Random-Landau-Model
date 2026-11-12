#include <math.h>
#include "../diag_wrappers/gpu_runtime.h"

#define PI2 (M_PI * 2.0)

void prepare_potential_coeff(float* impurity_x,
                             float* impurity_y,
                             float* impurity_intensity,
                             GpuComplex* Vmn,
                             GpuComplex* coeff_mn,
                             int dim_m,
                             int dim_n,
                             int off_head,
                             int impurity_num);
