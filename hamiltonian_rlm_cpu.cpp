#include"hamiltonian_rlm_cpu.h"
void *peer_solve_projected( void* peer_solve_params){
    int m, n, j, k, l ;
    int wfs_index, en_index, kx, ky;
    int idx_m,idx_kx,idx_jk,delta_jk;
    int n_phi,dim_n,dim_m,off_head;
    complex<dtype> coeff_theta_jm,sum;
    struct struct_solve *params=(struct struct_solve *)peer_solve_params;
    n_phi=params->n_phi;
    dim_n=params->dim_n;
    dim_m=params->dim_m;
    off_head=params->off_head;

    for(l = 0; l < params->theta_len; l++) {
        kx = params->theta_1;
        ky = params->theta_2[l];
        wfs_index = ((kx%2) * (params->n_mesh + 1) + ky) * n_phi*n_phi;
        en_index = ky * n_phi;
        idx_kx = kx*dim_n;
        for(j = 0; j < n_phi; j++)
            for(k = j; k < n_phi; k++) {
                params->wave_function[wfs_index + j * n_phi + k] = 0;
		idx_jk = j*n_phi*dim_n+k*dim_n;
		delta_jk=j-k-off_head;
                for(m = 0; m < dim_m; m++) {
		    idx_m = m*dim_n;
                    coeff_theta_jm=conj(params->coeff_m_theta[ky*dim_m+m])*params->coeff_jm[j*dim_m+m];
		    sum=0;
		    // kronecker delta (j,k-(n-off_head)), j-(k-(n-off_head))= N*n_phi, N=...-2,-1,0,1,2,...
		    // j-k \in [-n_phi,n_phi], (n-off_head) \in [-off_head,off_head], and off_head<n_phi
		    // possible (n-off_head) values: n-off_head = { -n_phi+(k-j), (k-j)} for some k and j's 
		    if(k-j<=off_head){
		      n=off_head+k-j;
                      sum += params->v_mn[idx_m + n] *params->coeff_m_theta[idx_kx + n];
		    }
		    if((n_phi+j-k)<=off_head){
		      n=off_head+k-j-n_phi;
                      sum += params->v_mn[idx_m + n] *params->coeff_m_theta[idx_kx + n];
		    }
		    // old version
                    //for(n = 0; n < dim_n; n++) 
		      //if(abs(delta_jk+n)%n_phi==0)
                      //sum += params->v_mn[idx_m + n] *params->coeff_m_theta[idx_kx + n];
                    params->wave_function[wfs_index + j * n_phi + k] += sum*coeff_theta_jm;

                }
                params->wave_function[wfs_index + k * n_phi + j] = conj(params->wave_function[wfs_index + j * n_phi + k]);
            }
        mkl_heevd(params->wave_function + wfs_index, params->energy + en_index, n_phi);
    }
    pthread_exit((void*)params);
}

void solve_projected(int theta_1, int n_mesh,int n_phi, int dim_m, int dim_n,int off_head, 
complex<dtype> *wave_function,complex<dtype> *coeff_m_theta, complex<dtype> *coeff_jm,complex<dtype> *v_mn, dtype *energy){
    int m, n, j, k, l ;
    int wfs_index, en_index, kx, ky;
    int idx_m,idx_kx,idx_jk,delta_jk;
    complex<dtype> coeff_theta_jm,sum;

    for(ky = 0; ky<= n_mesh; ky++) {
        kx = theta_1;
        wfs_index = ((kx%2) * (n_mesh + 1) + ky) * n_phi*n_phi;
        en_index = ky * n_phi;
        idx_kx = kx*dim_n;
        for(j = 0; j < n_phi; j++)
            for(k = j; k < n_phi; k++) {
                wave_function[wfs_index + j * n_phi + k] = 0;
		idx_jk = j*n_phi*dim_n+k*dim_n;
		delta_jk=j-k-off_head;
                for(m = 0; m < dim_m; m++) {
		    idx_m = m*dim_n;
                    coeff_theta_jm=conj(coeff_m_theta[ky*dim_m+m])*coeff_jm[j*dim_m+m];
		    sum=0;
		    // kronecker delta (j,k-(n-off_head)), j-(k-(n-off_head))= N*n_phi, N=...-2,-1,0,1,2,...
		    // j-k \in [-n_phi,n_phi], (n-off_head) \in [-off_head,off_head], and off_head<n_phi
		    // possible (n-off_head) values: n-off_head = { -n_phi+(k-j), (k-j)} for some k and j's 
		    if(k-j<=off_head){
		      n=off_head+k-j;
                      sum += v_mn[idx_m + n] *coeff_m_theta[idx_kx + n];
		    }
		    if((n_phi+j-k)<=off_head){
		      n=off_head+k-j-n_phi;
                      sum += v_mn[idx_m + n] *coeff_m_theta[idx_kx + n];
		    }
		    // old version
                    //for(n = 0; n < dim_n; n++) 
		      //if(abs(delta_jk+n)%n_phi==0)
                      //sum += params->v_mn[idx_m + n] *params->coeff_m_theta[idx_kx + n];
                    wave_function[wfs_index + j * n_phi + k] += sum*coeff_theta_jm;

                }
                wave_function[wfs_index + k * n_phi + j] = conj(wave_function[wfs_index + j * n_phi + k]);
            }
        mkl_heevd(wave_function + wfs_index, energy + en_index, n_phi);
    }

}
