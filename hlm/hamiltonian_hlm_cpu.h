#include"mkl_diag.h"
#include<complex>
#include<cstring>
#include<pthread.h>
#include"mkl.h"
#define PI2 (2.0*M_PI)
#if defined DP
#define dtype double
#define mkl_gemm3m cblas_zgemm3m
#else
#define dtype float
#define mkl_gemm3m cblas_cgemm3m
#endif
using namespace std;


/* Structure to store address of each variable in peer_solve_projected() */
struct struct_solve { 
    int lx;
    int ly;
    int Q;
    int band;
    int n_mesh;
    int theta_1;
    int *theta_2;
    int theta_len;
    complex<dtype> *wfs_full;
    dtype *potential;
    dtype *energy_theta;
    complex<dtype> *wfs;
    complex<dtype> *wfs_clean;
    complex<dtype> *truc_hamiltonian;
    dtype *truc_energy;
};

void get_wfs_ktheta(int theta_1,int theta_2,complex<dtype> *wfs);
// used to calculate the prestored wave_functions
void get_wfs_ktheta(complex<dtype> *wfs_ktheta);
void get_wfs_rspace(complex<dtype>* wfs,dtype* potential);
// get the twisted clean wave function by diagonalizating the real-space hamiltonian
void get_wfs_rspace(int theta_1, int theta_2,complex<dtype>* wfs);
void *peer_solve_projected(void* peer_solve_params);
