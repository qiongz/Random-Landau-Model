#include"potential_coeff.h"

__global__ void get_potential_matrix_elements(float *impurity_x,float *impurity_y,float *impurity_intensity, cuFloatComplex* Vmn, cuFloatComplex* coeff_mn, int dim_m, int dim_n, int off_head, int impurity_num){
    int m=blockIdx.x*blockDim.x+threadIdx.x;
    int n=blockIdx.y*blockDim.y+threadIdx.y;
    if(m<dim_m && n<dim_n){
    float kx=(m-off_head)*PI2;
    float ky=(n-off_head)*PI2;
    float phase;
    cuFloatComplex sum;
    sum.x=sum.y=0;
    for(int i=0;i<impurity_num;i++){
       phase=kx*impurity_x[i]+ky*impurity_y[i];
       sum.x+=cos(phase)*impurity_intensity[i];
       sum.y+=sin(phase)*impurity_intensity[i];
    }
    Vmn[m*dim_n+n]=cuCmulf(sum,coeff_mn[m*dim_n+n]);
    }
}

void prepare_potential_coeff(float* impurity_x, float * impurity_y, float *impurity_intensity, cuFloatComplex * Vmn, cuFloatComplex* coeff_mn, int dim_m, int dim_n, int off_head,int impurity_num) {
    dim3 grids(dim_m/4+1,dim_n/4+1);
    dim3 blocks(4,4);
    get_potential_matrix_elements<<<grids,blocks,0,NULL>>>(impurity_x, impurity_y,impurity_intensity,Vmn,coeff_mn,dim_m,dim_n,off_head, impurity_num);
    cudaStreamSynchronize(NULL);
}
