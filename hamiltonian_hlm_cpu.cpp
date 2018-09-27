#include"hamiltonian_hlm_cpu.h"
void get_wfs_ktheta(int lx, int ly, int Q, int band, int n_mesh, int theta_1,int theta_2,complex<dtype> *wfs) {
    // wavefunction_ktheta.shape= (n_mesh+1, n_mesh+1, lx/Q, ly, ly, lx/Q, Q), or (n_mesh+1, n_mesh+1, n_phi, n_sites)
    int kx,ky,p,q,m,j,r;
    int idx_1,idx;
    int n_sites=lx*ly;
    complex<dtype> wfs_value;
    complex<dtype> *hamiltonian = new complex<dtype> [Q * Q];
    dtype *energy = new dtype[Q];
    dtype norm;
    dtype inv_n_mesh=1.0/n_mesh;
    dtype inv_Q=1.0/Q;
    dtype pi2_inv_lx=PI2/lx;
    dtype pi2_inv_ly=PI2/ly;
    dtype kx_shifted,ky_shifted;
    dtype phase;
    // kx, ky in dimension n_phi
    for(kx=0; kx< lx/Q; kx++) {
        idx_1 = kx* ly* n_sites;
        kx_shifted=kx+theta_1*inv_n_mesh;
        for(ky=0; ky<ly; ky++) {
            ky_shifted=ky+theta_2*inv_n_mesh;
            // idx is the column index
            idx = ky * n_sites+idx_1;
            memset(hamiltonian,0,Q*Q*sizeof(complex<dtype>));
            for(p = 0; p < Q; p++) {
                hamiltonian[p*Q+p]= -2.0 * cos(ky_shifted * pi2_inv_ly + PI2 * p *inv_Q);
                if(p+1<Q)
                    hamiltonian[p*Q+p+1]=hamiltonian[(p+1)*Q+p]=-1;
            }
            hamiltonian[Q * (Q - 1)] = -complex<dtype>(cos(kx_shifted * Q * pi2_inv_lx), sin(kx_shifted * Q *pi2_inv_lx));
            hamiltonian[Q - 1] = conj(hamiltonian[Q * (Q-1)]);
            // diagonalize the k-space QxQ hamiltonian
            mkl_heevd(hamiltonian, energy, Q);
            // calculate the wavefunction in the n-th Landau level
            // loop over real-space indices r=j *lx +m * Q +q
            for(j = 0; j < ly; j++)
                for(m = 0; m < lx; m+=Q)
                    for(q = 0; q < Q; q++) {
                        phase=kx_shifted*pi2_inv_lx*m+ky_shifted*pi2_inv_ly*j;
                        wfs[idx+ j*lx +m+q] = complex<dtype>(cos(phase),-sin(phase)) * hamiltonian[band*Q+q];
                    }
            norm=0;
            // wavefunction normalization
            for(r = 0; r < n_sites; r++) {
                wfs_value = wfs[idx+r];
                norm += wfs_value.real() * wfs_value.real() + wfs_value.imag() * wfs_value.imag();
            }
            norm=1.0/sqrt(norm);
            for(r = 0; r < n_sites; r++)
                wfs[idx+r] *=norm;
        }
    }
    delete [] energy;
    delete [] hamiltonian;
}

// used to calculate the prestored wave_functions
void get_wfs_ktheta(int lx,int ly, int Q, int band, int n_mesh, complex<dtype> *wfs_ktheta) {
    // wavefunction_ktheta.shape= (n_mesh+1, n_mesh+1, lx/Q, ly, ly, lx/Q, Q), or (n_mesh+1, n_mesh+1, n_phi, n_sites)
    int theta_1,theta_2,kx,ky,p,q,i,m,j,r;
    int idx_1,idx_2,idx_3,idx;
    int n_sites=lx*ly;
    int n_phi=n_sites/Q;
    complex<dtype> wfs_value;
    complex<dtype> *hamiltonian = new complex<dtype> [Q * Q];
    dtype *energy = new dtype[Q];
    dtype norm;
    for(theta_1 = 0 ; theta_1 <= n_mesh; theta_1++) {
        idx_1= theta_1* (n_mesh+1) *n_sites* n_phi;
        for(theta_2 = 0; theta_2 <= n_mesh; theta_2++) {
            idx_2 = theta_2*n_sites* n_phi;
            // kx, ky in dimension n_phi
            for(kx=0; kx< lx/Q; kx++) {
               idx_3 = kx* ly* n_sites;
               for(ky=0; ky<ly; ky++) {
                   idx = ky * n_sites+idx_1+idx_2+idx_3;
                   memset(hamiltonian,0,Q*Q*sizeof(complex<dtype>));
                   for(p = 0; p < Q; p++) {
                      hamiltonian[p*Q+p]= -2.0 * cos((theta_2 * 1.0 / n_mesh + ky) * 2.0*M_PI / ly + 2.0*M_PI * p / Q);
                      if(p+1<Q)
                      hamiltonian[p*Q+p+1]=hamiltonian[(p+1)*Q+p]=-1;
                    }
                    hamiltonian[Q * (Q - 1)] = -complex<dtype>(cos((kx + theta_1 * 1.0 / n_mesh) * PI2 * Q / lx), sin((kx + theta_1 * 1.0 / n_mesh) * PI2 * Q / lx));
                    hamiltonian[Q - 1] = conj(hamiltonian[Q * (Q-1)]);
                    // diagonalize the k-space QxQ hamiltonian
                    mkl_heevd(hamiltonian, energy, Q);
                    // calculate the wavefunction in the lowest Landau level
                    // loop over real-space indices r=j *lx +m * Q +q
                    for(j = 0; j < ly; j++)
                        for(m = 0; m < lx / Q; m++)
                            for(q = 0; q < Q; q++) {
                                wfs_ktheta[idx+ j*lx + m* Q +q] = complex<dtype>(cos((kx + theta_1 * 1.0 /n_mesh) * PI2 * Q / lx * m + (ky + theta_2 * 1.0 / n_mesh) * PI2 / ly * j), -sin((kx + theta_1 * 1.0 / n_mesh) * PI2 * Q / lx * m + (ky + theta_2 * 1.0 / n_mesh) * PI2 /  ly* j)) * hamiltonian[q];
                            }
                    // wavefunction normalization
                    norm = 0;
                    for(r = 0; r < n_sites; r++) {
                        wfs_value = wfs_ktheta[idx+r];
                        norm += wfs_value.real() * wfs_value.real() + wfs_value.imag() * wfs_value.imag();
                    }
                    norm=sqrt(norm);
                    for(r = 0; r < n_sites; r++)
                        wfs_ktheta[idx+r] /=norm;
                }
            }
        }
    }
    delete [] energy;
    delete [] hamiltonian;
}

// get the twisted clean wave function by diagonalizating the real-space hamiltonian
void get_wfs_rspace(int lx, int ly, int Q, int band, int n_mesh, int theta_1, int theta_2,complex<dtype>* wfs){
    int n_sites=lx*ly;
    int n_phi=n_sites/Q;
    complex<dtype> *hamiltonian = new complex<dtype> [n_sites * n_sites];
    dtype *energy = new dtype[n_sites];
    memset(hamiltonian,0,sizeof(complex<dtype>)*n_sites*n_sites);
    for(int i=0;i<lx;i++)
      for(int j=0;j<ly;j++){
	 int ir=(i+1>=lx?i+1-lx:i+1);
	 int il=(i-1<0?i-1+lx:i-1);
	 hamiltonian[(i*ly+j)*n_sites+ir*ly+j]=-1;
	 hamiltonian[(ir*ly+j)*n_sites+i*ly+j]=-1;
	 hamiltonian[(i*ly+j)*n_sites+il*ly+j]=-1;
	 hamiltonian[(il*ly+j)*n_sites+i*ly+j]=-1;
	 int jr=(j+1>=ly?j+1-ly:j+1);
	 int jl=(j-1<0?j-1+ly:j-1);
         hamiltonian[(i*ly+j)*n_sites+i*ly+jr]=-complex<dtype>(cos(PI2*i/Q),sin(PI2*i/Q));
         hamiltonian[(i*ly+jr)*n_sites+i*ly+j]=-complex<dtype>(cos(PI2*i/Q),-sin(PI2*i/Q));
         hamiltonian[(i*ly+j)*n_sites+i*ly+jl]=-complex<dtype>(cos(PI2*i/Q),-sin(PI2*i/Q));
         hamiltonian[(i*ly+jl)*n_sites+i*ly+j]=-complex<dtype>(cos(PI2*i/Q),sin(PI2*i/Q));
       }
    for(int i=0;i<lx;i++){
	hamiltonian[(i*ly+ly-1)*n_sites+i*ly]=-complex<dtype>(cos(PI2*i/Q+theta_2*PI2/n_mesh),sin(PI2*i/Q+theta_2*PI2/n_mesh));
	hamiltonian[(i*ly)*n_sites+i*ly+ly-1]=-complex<dtype>(cos(PI2*i/Q+theta_2*PI2/n_mesh),-sin(PI2*i/Q+theta_2*PI2/n_mesh)); 
    }
    for(int j=0;j<ly;j++){
	hamiltonian[((lx-1)*ly+j)*n_sites+j]=-complex<dtype>(cos(theta_1*PI2/n_mesh),sin(theta_1*PI2/n_mesh));
	hamiltonian[j*n_sites+(lx-1)*ly+j]=-complex<dtype>(cos(theta_1*PI2/n_mesh),-sin(theta_1*PI2/n_mesh)); 
    }
    mkl_heevd(hamiltonian, energy, n_sites);
    for(int i=0;i<n_phi;i++)
       for(int j=0;j<n_sites;j++)
	     wfs[i*n_sites+j]=hamiltonian[(band%Q)*n_sites*n_phi+i*n_sites+j];
    delete [] energy,hamiltonian;
}

void *peer_solve_projected(void* peer_solve_params) {
    int theta_1,theta_2,l,i,j;
    int n_sites,n_phi,n_mesh;
    unsigned long dim_wfs;
    int wfs_idx,en_idx;
    peer_solve_paramsT *params = (peer_solve_paramsT*) peer_solve_params;
    complex<dtype> alpha{1.0f,0.0f};
    complex<dtype> beta{0.0f,0.0f};
    theta_1 = params->theta_1;
    n_sites=params->lx*params->ly; 
    n_phi=n_sites/params->Q;
    dim_wfs=n_sites*n_phi;
    n_mesh=params->n_mesh;
    for(l=0; l < params->theta_len; l++) {
        theta_2 = params->theta_2[l];
        wfs_idx = ((theta_1%2) * (n_mesh + 1) + theta_2) * dim_wfs;
        en_idx = theta_2* n_phi;
        // get the clean wave function for the lowest Landau level
        get_wfs_ktheta(params->lx, params->ly, params->Q, params->band, n_mesh,theta_1, theta_2,params->wfs_clean);
	//check if the k-space transformated wfs is correct 
	// by comparing it with r-space version
        //get_wfs_rspace(params->lx, params->ly, params->Q, params->band, n_mesh,theta_1, theta_2,params->wfs_clean);
        // project the potential to the lowest Landau level
	// since there's no diagonal matrix multiplication routine in MKL (like cdgem in cublas)
	for(i=0;i<n_phi;i++)
           for(j=0;j<n_sites;j++)
              params->wfs[i*n_sites+j]=params->potential[j]*params->wfs_clean[i*n_sites+j]; 
        mkl_gemm3m(CblasColMajor,CblasConjTrans,CblasNoTrans,n_phi,n_phi,n_sites,&alpha,params->wfs_clean,n_sites,params->wfs,n_sites,&beta,params->truc_hamiltonian,n_phi);	  
	// diagonalize the truncated hamiltonian
        mkl_heevd(params->truc_hamiltonian, params->energy_theta + en_idx, n_phi);
        // transform back to the full hamiltonian space
        mkl_gemm3m(CblasColMajor,CblasNoTrans,CblasNoTrans,n_sites,n_phi,n_phi,&alpha,params->wfs_clean,n_sites,params->truc_hamiltonian,n_phi,&beta,params->wfs_full+wfs_idx,n_sites);	
    }
    // clear tempory memory space
    free(params);
    pthread_exit((void*) 0);
}
