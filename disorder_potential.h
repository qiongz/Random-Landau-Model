#include"mkl.h"
#define METHOD_INT VSL_RNG_METHOD_UNIFORMBITS_STD
#define METHOD_FLOAT VSL_RNG_METHOD_UNIFORM_STD

void generate_disorder_potential(float * impurity_x, float *impurity_y, float* impurity_intensity, int impurity_num, float L1, float L2, unsigned long seed, int sample_offset);
