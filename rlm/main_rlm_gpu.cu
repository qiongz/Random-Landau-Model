#include<iostream>
#include<fstream>
#include<iomanip>
#include<complex>
#include<cuda.h>
#include<cuda_runtime.h>
#include"init_rlm.h"
#include"matrix_coefficients.h"
#include"disorder_potential.h"
#include"potential_coeff.h"
#include"hamiltonian_rlm_gpu.h"
#include"wfs_file.h"
#include"chern.h"
// macros for using different
// diagonalization routines
#if defined magma
#include"magma_diag.h"
#elif defined mkl
#include "mkl_diag.h"
#elif defined cusolver
#include"cusolver_diag.h"
#endif

using namespace std;
int main(int argc, char *argv[]) {
    /************************************  PARAMETERS INITIALIZATION  *******************************************/
    long n_phi,n_mesh,num_threads,n_sample,impurity_num,off_head, dim_m,dim_n;
    unsigned long seed,dim_wfs;
    float quanta_concentration,impurity_concentration,L1,L2;
    init(argc,argv,n_phi,quanta_concentration,impurity_concentration,n_mesh,n_sample,seed,num_threads);
    L1=sqrt(n_phi*quanta_concentration);
    L2=sqrt(n_phi*quanta_concentration);
    impurity_num = (impurity_concentration * n_phi / quanta_concentration);
    off_head =sqrt(impurity_num)/2+1;
    dim_m =dim_n= off_head * 2 + 1;
    dim_wfs=n_phi*n_phi;
    

    /*************************************  HOST MEMORY ALLOCATION    ******************************************/
    // Coefficients for calculating hamiltonian--> disorder independent coefficients
    complex<float> *coeff_mn, *coeff_m_theta, *coeff_jm;
    // disorder dependent coefficients
    // impurity positions
    float *impurity_x, *impurity_y;
    // impurity intensities
    float *impurity_intensity;
    // wave functions for all theta_1, theta_2
    complex<float> *wfs_full;
    // theta_1, theta_2 averaged eigenvalues
    float *energy_theta;
    // page-locked host memory, for data transfer
    // if use wfs write and read, reduces wfs in memory to 1*dim_wfs, otherwise 2*(n_mesh+1)*dim_wfs
    #ifdef wfsIO
    cudaHostAlloc((void**)&wfs_full,dim_wfs*sizeof(*wfs_full),cudaHostAllocMapped);
    #else
    unsigned long wfs_size=(n_mesh+1)*2*dim_wfs*sizeof(complex<float>);
    cudaHostAlloc((void**)&wfs_full,wfs_size,cudaHostAllocMapped);
    #endif
    cudaHostAlloc((void**)&coeff_m_theta,(n_mesh+1)*dim_m*sizeof(*coeff_m_theta),cudaHostAllocDefault);
    cudaHostAlloc((void**)&coeff_jm,n_phi*dim_m*sizeof(*coeff_jm),cudaHostAllocDefault);
    cudaHostAlloc((void**)&energy_theta,(n_mesh+1)*n_phi*sizeof(*energy_theta),cudaHostAllocDefault);

    coeff_mn=new complex<float>[dim_m*dim_n];
    impurity_x = new float[impurity_num];
    impurity_y = new float[impurity_num];
    impurity_intensity = new float[impurity_num];

    float *energy_levels = new float[n_phi];
    float *chern_numbers = new float[n_phi];
    float *chern_numbers_theta = new float[n_mesh *n_phi];

    /*************************************  DEVICE MEMORY ALLOCATION    ******************************************/
    complex<float> *dev_wfs;
    complex<float> *dev_coeff_m_theta_1;
    complex<float> *dev_coeff_m_theta_2;
    complex<float> *dev_coeff_jm;
    complex<float> *dev_coeff_mn;
    complex<float> *dev_Vmn;
    float *dev_energy;

    float * dev_impurity_x;
    float * dev_impurity_y;
    float * dev_impurity_intensity;


    // allocate device memory
    cudaMalloc ((void**)&dev_wfs,n_phi*n_phi*sizeof(*dev_wfs));
    cudaMalloc ((void**)&dev_coeff_m_theta_1,dim_m*sizeof(*dev_coeff_m_theta_1));
    cudaMalloc ((void**)&dev_coeff_m_theta_2,dim_m*sizeof(*dev_coeff_m_theta_2));
    cudaMalloc ((void**)&dev_coeff_jm,n_phi*dim_m*sizeof(*dev_coeff_jm));
    cudaMalloc ((void**)&dev_coeff_mn,dim_n*dim_m*sizeof(*dev_coeff_mn));
    cudaMalloc ((void**)&dev_Vmn,dim_m*dim_n*sizeof(*dev_Vmn));
    cudaMalloc ((void**)&dev_energy,n_phi*sizeof(*dev_energy));
    cudaMalloc ((void**)&dev_impurity_x,impurity_num*sizeof(*dev_impurity_x));
    cudaMalloc ((void**)&dev_impurity_y,impurity_num*sizeof(*dev_impurity_y));
    cudaMalloc ((void**)&dev_impurity_intensity,impurity_num*sizeof(*dev_impurity_intensity));


    /************************************* PTHREAD PARAMETERS  **************************************************/
    int **thds_ctheta;
    int *ctheta_len;
    // theta threads allocate
    thds_ctheta = new int*[num_threads];
    for(int i = 0; i < num_threads; i++) {
        thds_ctheta[i] = new int[n_mesh / num_threads + 1];
    }
    ctheta_len = new int[num_threads];

    // divide works into num_threads of pieces
    memset(ctheta_len, 0, num_threads*sizeof(int));
    for(int i = 0; i < n_mesh; i++)
        thds_ctheta[i% num_threads][ctheta_len[i%num_threads]++]=i;
    pthread_t *peer_thds = new pthread_t[num_threads];

    cudaEvent_t start,stop;
    float gpu_time,pot_time,coeff_time,chern_time,total_time;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    /***********************************  COEFFICIENTS INITIALIZATION  *******************************************/
    cudaEventRecord(start,0);
    // initialize and store the disorder independent coefficients
    prepare_coeff(coeff_mn, coeff_m_theta, coeff_jm,n_phi,off_head,dim_m,dim_n,n_mesh, L1, L2);
    // copy coefficients from host to device
    cudaMemcpy(dev_coeff_jm,coeff_jm,sizeof(*coeff_jm)*n_phi*dim_m,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_coeff_mn,coeff_mn,sizeof(*coeff_mn)*dim_n*dim_m,cudaMemcpyHostToDevice);
    // initialize data_file output
    ofstream fchern;
    fchern.open("energy_Chern_No.dat");

    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&coeff_time,start,stop);
    coeff_time/=1000.0;
    cerr<<std::left<<setw(40)<<"coefficients initialization time: "<<coeff_time<<" s \n";

    /*********************************** MAIN PROGRAM ************************************************************/
    for(int n = 0; n < n_sample; n++) {
        memset(energy_levels, 0, n_phi*sizeof(float));
        memset(chern_numbers, 0, n_phi*sizeof(float));

        cudaEventRecord(start,0);
        // initialize the disorder potential
        generate_disorder_potential(impurity_x,impurity_y,impurity_intensity,impurity_num,L1,L2,seed,n);
        cudaMemcpy(dev_impurity_x,impurity_x,sizeof(*impurity_x)*impurity_num,cudaMemcpyHostToDevice);
        cudaMemcpy(dev_impurity_y,impurity_y,sizeof(*impurity_y)*impurity_num,cudaMemcpyHostToDevice);
        cudaMemcpy(dev_impurity_intensity,impurity_intensity,
		sizeof(*impurity_intensity)*impurity_num,cudaMemcpyHostToDevice);
        // calculate the coefficient part from potential and transfer to the device
	prepare_potential_coeff(dev_impurity_x, dev_impurity_y, dev_impurity_intensity,
			(cuFloatComplex*)dev_Vmn,(cuFloatComplex*)dev_coeff_mn,dim_m,dim_n,off_head,impurity_num);

        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&pot_time,start,stop);
	pot_time/=1000.0;
	if(n==0)
          cerr<<setw(40)<<"potential initialization time: "<<pot_time<<" s\n";

        // solve the (theta_1=0)-line of wave functions
        int theta_1=0;
        for(int theta_2=0; theta_2<=n_mesh; theta_2++) {
	    int wfs_index;
            #if defined wfsIO
	    wfs_index = 0;
            #else 
            wfs_index = ((theta_1%2) * (n_mesh + 1) + theta_2) * dim_wfs;
	    #endif
            int en_index = theta_2 * n_phi;
            cudaMemcpy(dev_coeff_m_theta_1,coeff_m_theta+theta_1*dim_n,dim_n*sizeof(*coeff_m_theta),cudaMemcpyHostToDevice);
            cudaMemcpy(dev_coeff_m_theta_2,coeff_m_theta+theta_2*dim_n,dim_n*sizeof(*coeff_m_theta),cudaMemcpyHostToDevice);
            set_hamil_matrix((cuFloatComplex*)dev_coeff_m_theta_1,(cuFloatComplex*)dev_coeff_m_theta_2,
                             (cuFloatComplex*) dev_coeff_jm,(cuFloatComplex*)dev_Vmn, (cuFloatComplex*)dev_wfs,
                             n_phi,dim_n,off_head);
            #if defined magma
            magma_diag_hamil((magmaFloatComplex*)dev_wfs,energy_theta+en_index,n_phi);
            cudaMemcpy(wfs_full+wfs_index, dev_wfs, sizeof(*(dev_wfs))*dim_wfs, cudaMemcpyDeviceToHost);
            #elif defined mkl
            cudaMemcpy(wfs_full+wfs_index, dev_wfs, sizeof(*(dev_wfs))*dim_wfs, cudaMemcpyDeviceToHost);
            mkl_heevd(wfs_full+wfs_index,energy_theta+en_index,n_phi);
            #elif defined cusolver
            cusolver_diag_hamil((cuFloatComplex*)dev_wfs,dev_energy,n_phi);
            cudaMemcpy(wfs_full+wfs_index, dev_wfs, sizeof(*(dev_wfs))*dim_wfs, cudaMemcpyDeviceToHost);
            cudaMemcpy(energy_theta+en_index, dev_energy, sizeof(*(dev_energy))*n_phi, cudaMemcpyDeviceToHost);
            #endif
            #if defined wfsIO
	    write_wfs(theta_1%2,theta_2,wfs_full,dim_wfs);
            #endif
        }

        // average the (theta_1=0)-line of energy over theta_2
        // energy will be divided by (n_mesh+1)*(n_mesh+1) finally
        for(int i=0; i<n_phi; i++) {
            for(int j=0; j<=n_mesh; j++)
                energy_levels[i] += energy_theta[j* n_phi +i];
        }
        // solve theta_1-line of wfs from 1 to n_mesh, chern number calculations are also performed
        for(int theta_1=1; theta_1<=n_mesh; theta_1++) {
            if(theta_1==1)
              cudaEventRecord(start,0);
            // solve for another line
            for(int theta_2=0; theta_2<=n_mesh; theta_2++) {
		int wfs_index; 
                #if defined wfsIO
                wfs_index = 0;
                #else 
                wfs_index = ((theta_1%2) * (n_mesh + 1) + theta_2) * dim_wfs;
                #endif
                int en_index = theta_2 * n_phi;
                cudaMemcpy(dev_coeff_m_theta_1,coeff_m_theta+theta_1*dim_n,dim_n*sizeof(*coeff_m_theta),cudaMemcpyHostToDevice);
                cudaMemcpy(dev_coeff_m_theta_2,coeff_m_theta+theta_2*dim_n,dim_n*sizeof(*coeff_m_theta),cudaMemcpyHostToDevice);
                set_hamil_matrix((cuFloatComplex*)dev_coeff_m_theta_1,(cuFloatComplex*)dev_coeff_m_theta_2,
                                 (cuFloatComplex*) dev_coeff_jm,(cuFloatComplex*)dev_Vmn, (cuFloatComplex*)dev_wfs,
                                 n_phi,dim_n,off_head);
                #if defined magma
                magma_diag_hamil((magmaFloatComplex*)dev_wfs,energy_theta+en_index,n_phi);
                cudaMemcpy(wfs_full+wfs_index, dev_wfs, sizeof(*(dev_wfs))*dim_wfs, cudaMemcpyDeviceToHost);
                #elif defined mkl
                cudaMemcpy(wfs_full+wfs_index, dev_wfs, sizeof(*(dev_wfs))*dim_wfs, cudaMemcpyDeviceToHost);
                mkl_heevd(wfs_full+wfs_index,energy_theta+en_index,n_phi);
                #elif defined cusolver
                cusolver_diag_hamil((cuFloatComplex*)dev_wfs,dev_energy,n_phi);
                cudaMemcpy(wfs_full+wfs_index, dev_wfs, sizeof(*(dev_wfs))*dim_wfs, cudaMemcpyDeviceToHost);
                cudaMemcpy(energy_theta+en_index, dev_energy, sizeof(*(dev_energy))*n_phi, cudaMemcpyDeviceToHost);
                #endif
                #if defined wfsIO
	        write_wfs(theta_1%2,theta_2,wfs_full,dim_wfs);
                #endif
            }
            // average the energy
            for(int i=0; i<n_phi; i++) {
                for(int j=0; j<n_mesh+1; j++)
                    energy_levels[i] += energy_theta[j* n_phi +i];
            }
	    if(theta_1==1 && n==0){
              cudaEventRecord(stop,0);
              cudaEventSynchronize(stop);
              cudaEventElapsedTime(&gpu_time,start,stop);
	      gpu_time/=1000.0;
              cerr<<setw(40)<<"diagonalization time/k-point: "<<gpu_time/(n_mesh+1)<<" s\n";
              cudaEventRecord(start,0);
	    }
            #if defined wfsIO
            cal_Chern_wfs_IO(chern_numbers_theta,n_phi,n_phi, n_mesh, theta_1);
            #else
            // calculate Chern numbers for two lines
	    /*
            for(int id = 0; id < num_threads; id++) {
                peer_Chern_paramsT *params;
                params = (peer_Chern_paramsT *) malloc(sizeof(peer_Chern_paramsT));
                params-> n_phi = n_phi;
                params-> dim_vec = n_phi;
                params-> n_mesh = n_mesh;
                params-> theta_1 = theta_1;
                params-> theta_2 = thds_ctheta[id];
                params-> theta_len = ctheta_len[id];
                params-> wave_function = wfs_full;
                params-> chern_numbers_theta = chern_numbers_theta;
                pthread_create(&(peer_thds[id]), NULL, peer_cal_Chern,  (void*)params);
            }
            // join all the threads
            for(int id = 0; id < num_threads; id++)
                pthread_join(peer_thds[id], NULL);
		*/
            cal_Chern(wfs_full,chern_numbers_theta,n_phi,n_phi, n_mesh, theta_1);
            #endif
	    if(theta_1==1 && n==0){
              cudaEventRecord(stop,0);
              cudaEventSynchronize(stop);
              cudaEventElapsedTime(&chern_time,start,stop);
	      chern_time/=1000.0;
	      total_time=coeff_time+(pot_time+(chern_time*n_mesh+gpu_time*(n_mesh+1)))*n_sample;
              cerr<<setw(40)<<"Chern No. calculation time/k-point: "<<chern_time/(n_mesh+1)<<" s\n";
	      cerr<<setw(40)<<"estimated total time:"<<std::right<<setw(3)<<setfill('0')<<(int(total_time))/3600<<setw(1)<<":";
	      cerr<<setw(2)<<setfill('0')<<((int(total_time))%3600)/60<<setw(1)<<":";
	      cerr<<setw(2)<<setfill('0')<<(int(total_time))%60<<endl<<setfill(' ');
	      cerr<<std::left<<setw(40)<<"estimated # of samples/hour:"<<n_sample/(total_time/3600.0)<<endl;
	      cerr<<setw(40)<<"# module "<<" %(time)"<<endl;
	      cerr<<setw(40)<<"coefficients"<<coeff_time/total_time*100.0<<"\n";
	      cerr<<setw(40)<<"potential"<<pot_time*n_sample/total_time*100.0<<"\n";
	      cerr<<setw(40)<<"diagonalization"<<gpu_time*(n_mesh+1)*n_sample/total_time*100.0<<endl;
	      cerr<<setw(40)<<"Chern No. calculation"<<chern_time*n_mesh*n_sample/total_time*100.0<<endl;
	    }
            //serial version
            //cal_Chern(wfs_full,chern_numbers_theta,n_phi,n_phi, n_mesh, theta_1);
            // collect the chern number contributions from the theta_1 line
            for(int k=0; k<n_phi; k++) {
                for(int i=0; i<n_mesh; i++)
                    chern_numbers[k]+=chern_numbers_theta[k*n_mesh+i];
            }
        }
        //average the energy levels by dividing (n_mesh+1)*(n_mesh+1) at the final stage
        for(int k=0; k<n_phi; k++) {
            energy_levels[k]/=(n_mesh+1)*(n_mesh+1);
            chern_numbers[k]+=1.0/n_phi;
        }
        // print averaged energy levels
        for(int i = 0; i < n_phi; i++)
            fchern << setw(8) << n << " " << setw(15) << setprecision(8) << energy_levels[i] << " " << chern_numbers[i] << endl;
    }

    /***********************************  CLEAR DEVICE MEMORY  ***************************************/
    cudaFree(dev_wfs);
    cudaFree(dev_coeff_jm);
    cudaFree(dev_coeff_m_theta_1);
    cudaFree(dev_coeff_m_theta_2);
    cudaFree(dev_coeff_mn);
    cudaFree(dev_Vmn);
    cudaFree(dev_energy);
    cudaFree(dev_impurity_x);
    cudaFree(dev_impurity_y);
    cudaFree(dev_impurity_intensity);

    /***********************************  CLEAR HOST MEMORY  ****************************************/
    delete[] chern_numbers;
    delete[] energy_levels;
    delete[] impurity_intensity;
    delete[] impurity_y;
    delete[] impurity_x;
    delete[] coeff_mn;

    cudaFreeHost(wfs_full);
    cudaFreeHost(coeff_jm);
    cudaFreeHost(coeff_m_theta);
    cudaFreeHost(energy_theta);

    fchern.close();

    return 0;
}
