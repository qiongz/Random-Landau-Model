#include"mkl_diag.h"
#include<complex>
#include<pthread.h>
#if defined DP
#define dtype double
#else
#define dtype float
#endif
using namespace std;

/* Structure to store address of each variable in peer_solve_projected() */
typedef struct {
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
} peer_solve_paramsT;

void *peer_solve_projected( void* peer_solve_params);
