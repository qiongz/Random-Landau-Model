#include<complex>
#include<cmath>
#if defined DP
#define dtype double
#else
#define dtype float
#endif
#define PI2 (M_PI*2.0)
using namespace std;
void prepare_coeff(complex<dtype> *coeff_mn, complex<dtype> *coeff_m_theta, complex<dtype> *coeff_jm, 
		int n_phi, int off_head, int dim_m, int dim_n, int n_mesh, dtype L1, dtype L2);
void get_potential_coeff(dtype* impurity_x, dtype * impurity_y, dtype *impurity_intensity, complex<dtype> * Vmn, 
		complex<dtype>* coeff_mn, int dim_m, int dim_n, int off_head,int impurity_num);
