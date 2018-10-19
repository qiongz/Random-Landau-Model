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
    fchern.open("energy_Chern_No.txt");

    double kwfs_time,diag_time,chern_time;
    long total_time;

    // initialize random number generator
    VSLStreamStatePtr rndStream;
    vslNewStream(&rndStream, VSL_BRNG_MT19937, seed);

    void *thread_results;
    pthread_t *peer_thds = new pthread_t[num_threads];
    struct struct_solve params[num_threads]; 
    struct struct_Chern Chern_params[num_threads];

    /*********************************** MAIN PROGRAM ************************************************************/
    for(n = 0; n < n_sample; n++) {
	memset(chern_numbers,0,sizeof(dtype)*n_phi);
	memset(energy_levels,0,sizeof(dtype)*n_phi);
	// initialize the disorder potential
        VS_RNG_DTYPE(METHOD_DTYPE, rndStream, n_sites, potential,-1,1);

        auto _t1=std::chrono::high_resolution_clock::now();
        // solve the (theta_1=0)-line of wave functions
        // solve all the disroder potential hamiltonians
        for(id = 0; id < num_threads; id++) {
	    params[id].lx=lx;
	    params[id].ly=ly;
	    params[id].Q=Q;
	    params[id].band=band;
	    params[id].n_mesh=n_mesh;
	    params[id].theta_1 = 0;
            params[id].theta_2 = thds_theta[id];
            params[id].theta_len = theta_len[id];
            params[id].wfs_full = wfs_full;
            params[id].energy_theta = energy_theta;
            params[id].potential = potential;
	    params[id].wfs=wfs+id*dim_wfs;
	    params[id].wfs_clean=wfs_clean+id*dim_wfs;
	    params[id].truc_hamiltonian=truc_hamiltonian+id*n_phi*n_phi;
            pthread_create(&(peer_thds[id]), NULL, peer_solve_projected,  (void* )&params[id]);
        }
        // join all the threads
        for(id = 0; id < num_threads; id++) 
            pthread_join(peer_thds[id], &thread_results);
        
        auto _t2=std::chrono::high_resolution_clock::now();
        diag_time=std::chrono::duration_cast<chrono::microseconds>(_t2-_t1).count()/1.0E6;
	if(n==0)
	cerr<<std::left<<setw(40)<<"diagonalization time/k-point: "<<diag_time/(n_mesh+1)<<" s\n";
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
	    params[id].lx=lx;
	    params[id].ly=ly;
	    params[id].Q=Q;
	    params[id].band=band;
	    params[id].n_mesh=n_mesh;
	    params[id].theta_1 = 0;
            params[id].theta_2 = thds_theta[id];
            params[id].theta_len = theta_len[id];
            params[id].wfs_full = wfs_full;
            params[id].energy_theta = energy_theta;
            params[id].potential = potential;
	    params[id].wfs=wfs+id*dim_wfs;
	    params[id].wfs_clean=wfs_clean+id*dim_wfs;
	    params[id].truc_hamiltonian=truc_hamiltonian+id*n_phi*n_phi;
            pthread_create(&(peer_thds[id]), NULL, peer_solve_projected,  (void* )&params[id]);
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

            _t1=std::chrono::high_resolution_clock::now();
            //cal_Chern(wfs_full,chern_numbers_theta,n_phi,n_sites,n_mesh,theta_1);
            // calculate Chern numbers
            for(id = 0; id < num_threads; id++) {
	        Chern_params[id].n_phi = n_phi;
	        Chern_params[id].dim_vec = n_sites;
	        Chern_params[id].n_mesh = n_mesh;
		Chern_params[id].theta_1 = theta_1;
                Chern_params[id].theta_2 = thds_theta[id];
                Chern_params[id].theta_len = theta_len[id];
                Chern_params[id].wave_function = wfs_full;
                Chern_params[id].chern_numbers_theta = chern_numbers_theta;
                pthread_create(&(peer_thds[id]), NULL, peer_cal_Chern,  (void*)&Chern_params[id]);
            }
            // join the threads
            for(id = 0; id < num_threads; id++)
                pthread_join(peer_thds[id], NULL);
            _t2=std::chrono::high_resolution_clock::now();
            chern_time=std::chrono::duration_cast<chrono::microseconds>(_t2-_t1).count()/1.0E6;
           
	    if(theta_1==1 && n==0){
	    total_time=round((chern_time*n_mesh+diag_time*(n_mesh+1))*n_sample);
            cerr<<setw(40)<<"Chern No. calculation time/k-point: "<<chern_time/(n_mesh+1)<<" s\n";
	    cerr<<setw(40)<<"estimated total time:"<<std::right<<setw(3)<<setfill('0')<<total_time/3600<<setw(1)<<":";
	    cerr<<setw(2)<<std::right<<setfill('0')<<(total_time%3600)/60<<setw(1)<<":";
	    cerr<<setw(2)<<std::right<<setfill('0')<<total_time%60<<endl<<setfill(' ');
	    cerr<<std::left<<setw(40)<<"estimated # of samples/hour:"<<n_sample/(total_time/3600.0)<<endl;
	    cerr<<setw(40)<<"# module"<<"%(time)"<<endl;
	    cerr<<setw(40)<<"diagonalization"<<diag_time*(n_mesh+1)*n_sample/total_time*100.0<<endl;
	    cerr<<setw(40)<<"Chern No. calculation"<<chern_time*n_mesh*n_sample/total_time*100.0<<endl;
	    }
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
