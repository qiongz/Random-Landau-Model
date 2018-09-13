#include<complex>
#include<vector>
#include<pthread.h>
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#if defined DP
#define dtype double
#define DOT zdotc
#else
#define dtype float
#define DOT cdotc
#endif
using namespace std;

/* Structure to store address of each variable in peer_cal_Chern() */
typedef struct { 
    int n_phi;
    int n_mesh;
    int theta_1;
    int *theta_2;
    int theta_len;
    complex<dtype> *wave_function;
    dtype *chern_numbers_theta; 
} peer_Chern_paramsT;


void cal_Chern( complex<dtype>* wfs_full,dtype* chern_numbers_theta, int n_phi, int n_mesh,int theta_1);
void cal_Chern_wfs_IO(dtype* chern_numbers_theta, int n_phi, int n_mesh, int theta_1);
void *peer_cal_Chern( void *peer_Chern_params);
