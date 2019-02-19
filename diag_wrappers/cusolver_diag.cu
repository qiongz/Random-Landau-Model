#include"cusolver_diag.h"

void cusolver_diag_hamil(cuFloatComplex* dev_wfs, float*dev_energy,int n_phi){
    cusolverDnHandle_t cusolverH ;
    cusolverDnCreate (& cusolverH );
    cusolverDnSetStream(cusolverH,NULL);
    cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR;
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    cuFloatComplex*dev_work;
    int *info;
    int lwork;
    // querying workspace
    cusolverDnCheevd_bufferSize (cusolverH,jobz, uplo, n_phi, dev_wfs, n_phi, dev_energy, & lwork );
    cudaMalloc (( void **)& dev_work, sizeof (*dev_work)* lwork );
    cudaMalloc((void**)&info, sizeof(*info));
    // diagonalization the hamiltonian
    cusolverDnCheevd(cusolverH,jobz,uplo,n_phi,dev_wfs,n_phi,dev_energy,dev_work,lwork,info);
    cudaStreamSynchronize(NULL);
    cusolverDnDestroy(cusolverH);
    cudaFree(dev_work);
    cudaFree(info);
}

