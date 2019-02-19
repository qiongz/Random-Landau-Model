#include"mkl.h"
#define METHOD_INT VSL_RNG_METHOD_UNIFORMBITS_STD
#define METHOD_DTYPE VSL_RNG_METHOD_UNIFORM_STD
#if defined DP
#define dtype double
#define VS_RNG_DTYPE vdRngUniform
#else
#define dtype float
#define VS_RNG_DTYPE vsRngUniform
#endif

void generate_disorder_potential(dtype * impurity_x, dtype *impurity_y, dtype* impurity_intensity, int impurity_num, dtype L1, dtype L2, unsigned long seed, int sample_offset);
