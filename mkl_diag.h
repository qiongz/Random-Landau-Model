#include<complex>
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include"mkl.h"
#if defined DP
#define dtype double
#define HEEVD LAPACKE_zheevd
#else
#define dtype float
#define HEEVD LAPACKE_cheevd
#endif
using namespace std;
void mkl_cheevd(complex<dtype> *hamiltonian, dtype * energy, int l);
