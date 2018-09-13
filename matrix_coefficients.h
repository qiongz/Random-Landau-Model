#include<complex>
#include<cmath>
#define PI2 (M_PI*2.0)
using namespace std;
void prepare_coeff(complex<float> *coeff_mn, complex<float> *coeff_m_theta, complex<float> *coeff_jm, 
		int n_phi, int off_head, int dim_m, int dim_n, int n_mesh, float L1, float L2);
void get_potential_coeff(float* impurity_x, float * impurity_y, float *impurity_intensity, complex<float> * Vmn, 
		complex<float>* coeff_mn, int dim_m, int dim_n, int off_head,int impurity_num);
