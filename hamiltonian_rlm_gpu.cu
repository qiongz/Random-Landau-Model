#include"hamiltonian_gpu.h"
__device__ __host__ cuFloatComplex  operator*(cuFloatComplex a,cuFloatComplex b) {
    return cuCmulf(a,b);
}
__device__ __host__ cuFloatComplex  operator+(cuFloatComplex a,cuFloatComplex b) {
    return cuCaddf(a,b);
}

__global__ void get_matrix_elements(const cuFloatComplex* coeff_jm, const cuFloatComplex* coeff_m_theta_1,const cuFloatComplex* coeff_m_theta_2, const cuFloatComplex* Vmn, cuFloatComplex* wfs, int dim_n, int off_head) {
    int nphi=gridDim.x*blockDim.x;
    int j=blockIdx.x*blockDim.x+threadIdx.x;
    int k=blockIdx.y*blockDim.y+threadIdx.y;
    cuFloatComplex sum;
    int wfs_idx=j*nphi+k;
    int idx_m,n,m;
    wfs[wfs_idx].x=wfs[wfs_idx].y=0;
    if(k>=j) {
        for( m = 0; m < dim_n; m++) {
            sum.x=sum.y=0;
            idx_m=m*dim_n;
            if(k-j<=off_head ) {
                n=off_head+k-j;
                sum=sum+Vmn[idx_m + n]*cuConjf(coeff_m_theta_1[n]);
            }
            if((nphi+j-k)<=off_head) {
                n=off_head+k-j-nphi;
                sum=sum+Vmn[idx_m + n]*cuConjf(coeff_m_theta_1[n]);
            }
            // old version of Kronecker delta
            /*
                for(int n = 0; n < dim_n; n++)
                    if(abs(j-k+n-off_head)%nphi==0) {
                        sum=sum+Vmn[idx_m + n]*coeff_m_theta_1[n];
                    }
            */
            wfs[wfs_idx]= wfs[wfs_idx]+ sum* coeff_jm[j*dim_n+m]*coeff_m_theta_2[m];
        }
        wfs[k*nphi+j]=cuConjf(wfs[wfs_idx]);
    }

}

void set_hamil_matrix(const cuFloatComplex* dev_coeff_m_theta_1,const cuFloatComplex* dev_coeff_m_theta_2, const cuFloatComplex* dev_coeff_jm, const cuFloatComplex *dev_Vmn, cuFloatComplex* dev_wfs, int n_phi,int dim_n,int off_head){
    // n_phi should be multiplier of 8 
    dim3 grids(n_phi/8,n_phi/8);
    dim3 blocks(8,8);
    // calculating the hamiltonian matrix elements
    get_matrix_elements<<<grids,blocks,0,NULL>>>(dev_coeff_jm,dev_coeff_m_theta_1,dev_coeff_m_theta_2,dev_Vmn,dev_wfs,dim_n,off_head);
    cudaStreamSynchronize(NULL);
}

