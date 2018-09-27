#include<iostream>
#include<fstream>
#include<iomanip>
#include<complex>
#include<cmath>
#include<cstring>
#include"chern.h"
#include"wfs_file.h"
#include"hamiltonian_hlm_cpu.h"
#include"mkl_diag.h"
#include"init_hlm.h"
#include<pthread.h>
#if defined DP
#define dtype double
#define VS_RNG_DTYPE vdRngUniform
#else
#define dtype float
#define VS_RNG_DTYPE vsRngUniform
#endif
#define METHOD_DTYPE VSL_RNG_METHOD_UNIFORM_STD

int main(int argc, char *argv[]) {
    /************************************  PARAMETERS INITIALIZATION  *******************************************/
    int lx,ly,Q,band,n_phi,n_sites,n_mesh,num_threads,n_sample;
    unsigned long seed,dim_wfs;
    int i,j,k,id,n;
    init(argc, argv, lx,ly,Q,band,n_mesh,n_sample,seed,num_threads);
    n_sites=lx*ly;
    n_phi=n_sites/Q;
    dim_wfs=n_phi*n_sites;

    /*************************************  MEMORY ALLOCATION    ******************************************/
    complex<dtype> *wfs_full = new complex<dtype>[(n_mesh+1)*2* dim_wfs];
    dtype *potential = new dtype[n_sites];
    complex<dtype>* wfs=new complex<dtype> [ num_threads * dim_wfs];
    complex<dtype>* wfs_clean = new complex<dtype> [num_threads * dim_wfs];
    complex<dtype>* truc_hamiltonian= new complex<dtype> [ num_threads* n_phi *n_phi];

    // theta_1, theta_2 averaged eigenvalues
    dtype *energy_theta = new dtype [ (n_mesh+1)*n_phi];
    dtype *energy_levels = new dtype [n_phi];
    dtype *chern_numbers_theta = new dtype [ n_mesh*n_phi];
    dtype *chern_numbers = new dtype [n_phi];

    int **thds_theta, *theta_len;
    thds_theta = new int*[num_threads];
    for(i = 0; i < num_threads; i++) {
        thds_theta[i] = new int[n_mesh / num_threads + 1];
    }
    theta_len = new int[num_threads];
    for(i=0; i<num_threads; i++)
        theta_len[i]=0;

    // ctheta_1, ctheta_2 for peer_cal_chern
    for(i = 0; i < n_mesh; i++) {
        thds_theta[i% num_threads][theta_len[i% num_threads]++] = i;
    }

    ofstream fchern;
    fchern.open("chern_hlm_cpu.dat");

    double kwfs_time,diag_time,chern_time,total_time;

    // initialize random number generator
    VSLStreamStatePtr rndStream;
    vslNewStream(&rndStream, VSL_BRNG_MT19937, seed);

    void *thread_results;
    pthread_t *peer_thds = new pthread_t[num_threads];

    /*********************************** MAIN PROGRAM ************************************************************/
    for(n = 0; n < n_sample; n++) {
	memset(chern_numbers,0,sizeof(dtype)*n_phi);
	memset(energy_levels,0,sizeof(dtype)*n_phi);
	// initialize the disorder potential
        VS_RNG_DTYPE(METHOD_DTYPE, rndStream, n_sites, potential,-1,1);
        // solve the (theta_1=0)-line of wave functions
        // solve all the disroder potential hamiltonians
        for(id = 0; id < num_threads; id++) {
            peer_solve_paramsT * params;
            params = (peer_solve_paramsT *) malloc(sizeof(peer_solve_paramsT));
	    params->lx=lx;
	    params->ly=ly;
	    params->Q=Q;
	    params->band=band;
	    params->n_mesh=n_mesh;
	    params->theta_1 = 0;
            params->theta_2 = thds_theta[id];
            params->theta_len = theta_len[id];
            params->wfs_full = wfs_full;
            params->energy_theta = energy_theta;
            params->potential = potential;
	    params->wfs=wfs+id*dim_wfs;
	    params->wfs_clean=wfs_clean+id*dim_wfs;
	    params->truc_hamiltonian=truc_hamiltonian+id*n_phi*n_phi;
            pthread_create(&(peer_thds[id]), NULL, peer_solve_projected,  (void* )params);
        }
        // join all the threads
        for(id = 0; id < num_threads; id++) 
            pthread_join(peer_thds[id], &thread_results);
        
        // average the energy over the first line of theta_2
        // energy will be divided by (n_mesh+1)*(n_mesh+1) finally
        for(i = 0; i < n_phi; i++) {
            for(j = 0; j < n_mesh; j++)
                energy_levels[i] += energy_theta[j* n_phi + i];
            energy_levels[i]+= energy_theta[i];
        }
        // periodic boundary conditions for wave functions
        for(k=0; k<dim_wfs; k++)
            wfs_full[n_mesh*dim_wfs+k]=wfs_full[k];

        // solve theta_1-line of wfs from 1 to n_mesh, chern number calculations are also performed
        for(int theta_1=1; theta_1<=n_mesh; theta_1++) {
            for(id = 0; id < num_threads; id++) {
                peer_solve_paramsT * params;
                params = (peer_solve_paramsT *) malloc(sizeof(peer_solve_paramsT));
	        params->lx=lx;
	        params->ly=ly;
	        params->Q=Q;
	        params->band=band;
	        params->n_mesh=n_mesh;
	        params->theta_1 = theta_1;
                params->theta_2 = thds_theta[id];
                params->theta_len = theta_len[id];
                params->wfs_full = wfs_full;
                params->energy_theta = energy_theta;
                params->potential = potential;
	        params->wfs=wfs+id*dim_wfs;
   	        params->wfs_clean=wfs_clean+id*dim_wfs;
	        params->truc_hamiltonian=truc_hamiltonian+id*n_phi*n_phi;
                pthread_create(&(peer_thds[id]), NULL, peer_solve_projected,  (void* )params);
            }
            // join all the threads
            for(id = 0; id < num_threads; id++) 
                pthread_join(peer_thds[id], &thread_results);
            
            // average the energy over the first line of theta_2
            for(i = 0; i < n_phi; i++) {
                for(j = 0; j < n_mesh; j++)
                    energy_levels[i] += energy_theta[j* n_phi + i];
                energy_levels[i]+= energy_theta[i];
            }
            // periodic boundary conditions
            for(k=0; k<dim_wfs; k++)
                wfs_full[((theta_1%2)*(n_mesh+1)+n_mesh)*dim_wfs+k]=wfs_full[((theta_1%2)*(n_mesh+1))*dim_wfs+k];

            // calculate Chern numbers
            for(id = 0; id < num_threads; id++) {
                peer_Chern_paramsT *params;
                params = (peer_Chern_paramsT *) malloc(sizeof(peer_Chern_paramsT));
	        params-> n_phi = n_phi;
	        params-> dim_vec = n_sites;
	        params-> n_mesh = n_mesh;
		params-> theta_1 = theta_1;
                params-> theta_2 = thds_theta[id];
                params-> theta_len = theta_len[id];
                params-> wave_function = wfs_full;
                params-> chern_numbers_theta = chern_numbers_theta;
                pthread_create(&(peer_thds[id]), NULL, peer_cal_Chern,  (void*)params);
            }
            // join the threads
            for(id = 0; id < num_threads; id++)
                pthread_join(peer_thds[id], NULL);

            // collect the chern number contributions from the theta-1= (global_theta_1-1) line
            for(k=0; k<n_phi; k++) {
                for(i=0; i<n_mesh; i++)
                    chern_numbers[k]+=chern_numbers_theta[k*n_mesh+i];
            }
        }
        // average the energy levels by dividing (n_mesh+1)*(n_mesh+1) at the final stage
        for(k=0; k<n_phi; k++) 
           energy_levels[k]/=(n_mesh+1)*(n_mesh+1);
        
        // print averaged energy levels
        for(i = 0; i < n_phi; i++){
            fchern << setw(8) << n << " " << setw(15) << setprecision(8) << energy_levels[i] << " " << chern_numbers[i] << endl;
	}
    }

    delete[] wfs;
    delete[] wfs_clean;
    delete[] truc_hamiltonian;
    delete[] wfs_full;
    delete[] energy_levels;
    delete[] chern_numbers;
    delete[] potential;
    for(id=0; id<num_threads; id++)
        delete[] thds_theta[id];
    delete[] thds_theta;
    delete[] theta_len;
    vslDeleteStream(&rndStream);
    fchern.close();

    return 0;
}
