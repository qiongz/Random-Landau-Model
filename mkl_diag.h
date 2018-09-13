#include<complex>
#define MKL_Complex8 std::complex<float>
#include"mkl.h"
using namespace std;
void mkl_cheevd(complex<float> *hamiltonian, float * energy, int l);
