#include<complex>
#include<vector>
#include<pthread.h>
#define MKL_Complex8 std::complex<float>
#include "mkl.h"
using namespace std;

/* Structure to store address of each variable in peer_cal_Chern() */
typedef struct { 
    int n_phi;
    int n_mesh;
    int theta_1;
    int *theta_2;
    int theta_len;
    complex<float> *wave_function;
    float *chern_numbers_theta; 
} peer_Chern_paramsT;


void cal_Chern( complex<float>* wfs_full,float* chern_numbers_theta, int n_phi, int n_mesh,int theta_1);
void cal_Chern_wfs_IO(float* chern_numbers_theta, int n_phi, int n_mesh, int theta_1);
void *peer_cal_Chern( void *peer_Chern_params);
