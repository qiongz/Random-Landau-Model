#include"../diag_wrappers/mkl_diag.h"
#include<complex>
#include<pthread.h>
#if defined DP
#define dtype double
#else
#define dtype float
#endif
using namespace std;

/* Structure to store address of each variable in peer_solve_projected() */
struct struct_solve {
    int theta_1;
    int n_phi;
    int n_mesh;
    int dim_m;
    int dim_n;
    int off_head;
    int *theta_2;
    int theta_len;
    dtype *energy;
    complex<dtype>*wave_function;
    complex<dtype> *v_mn;
    complex<dtype> *coeff_jm;
    complex<dtype> *coeff_m_theta;
};

// multi-threaded version
void *peer_solve_projected( void* peer_solve_params);
void solve_projected(int theta_1, int n_mesh,int n_phi, int dim_m, int dim_n,int off_head, 
complex<dtype> *wave_function,complex<dtype> *coeff_m_theta, complex<dtype> *coeff_jm,complex<dtype> *v_mn, dtype *energy);
