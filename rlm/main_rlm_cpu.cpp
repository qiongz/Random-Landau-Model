#include<iostream>
#include<fstream>
#include<iomanip>
#include<complex>
#include<cmath>
#include<cstring>
#include"../diag_wrappers/chern.h"
#include"../diag_wrappers/wfs_file.h"
#include"../diag_wrappers/mkl_diag.h"
#include"hamiltonian_rlm_cpu.h"
#include"disorder_potential.h"
#include"matrix_coefficients.h"
#include"init_rlm.h"
#include<pthread.h>
#if defined DP
#define dtype double
#else
#define dtype float
#endif

using namespace std;

int main(int argc, char *argv[]) {
    /************************************  PARAMETERS INITIALIZATION  *******************************************/
    long n_phi,n_mesh,num_threads,n_sample,impurity_num,off_head, dim_m,dim_n;
    unsigned long seed,dim_wfs;
    dtype quanta_concentration,impurity_concentration,L1,L2;
    init(argc,argv,n_phi,quanta_concentration,impurity_concentration,n_mesh,n_sample,seed,num_threads);
    L1=sqrt(n_phi*quanta_concentration);
    L2=sqrt(n_phi*quanta_concentration);
    impurity_num = (impurity_concentration * n_phi / quanta_concentration);
    off_head =sqrt(impurity_num)/2+1;
    dim_m =dim_n= off_head * 2 + 1;
    dim_wfs=n_phi*n_phi;

    /*************************************  MEMORY ALLOCATION    ******************************************/
    // Coefficients for calculating hamiltonian--> disorder independent coefficients
    complex<dtype> *coeff_mn, *coeff_m_theta, *coeff_jm;

    // disorder dependent coefficients
    // impurity positions
    dtype *impurity_x, *impurity_y;
    // impurity intensities
    dtype *impurity_intensity;
    unsigned int *i_buffer;
    complex<dtype> *v_mn;
    // wave functions for all theta_1, theta_2
    complex<dtype> *wave_functions;
    // theta_1, theta_2 averaged eigenvalues
    dtype *energy_theta;
    dtype *energy_levels;
    dtype *chern_numbers;
    dtype *chern_numbers_theta;
    // paramsters for peer_solve
    int **thds_theta, **thds_ctheta;
    int *ctheta_len,*theta_len;
    unsigned long wfs_size=(n_mesh+1)*2*n_phi*n_phi;

    // allocate memory space
    coeff_mn = new complex<dtype> [dim_m * dim_n];
    coeff_m_theta = new complex<dtype> [(n_mesh + 1)  * dim_m];
    coeff_jm = new complex<dtype> [n_phi * dim_m];
    impurity_x = new dtype[impurity_num];
    impurity_y = new dtype[impurity_num];
    impurity_intensity = new dtype[impurity_num];
    i_buffer = new unsigned int[impurity_num];
    v_mn = new complex<dtype>[dim_m * dim_n];
    wave_functions = new complex<dtype>[wfs_size];
    energy_theta = new dtype[(n_mesh + 1) * n_phi];
    energy_levels = new dtype[n_phi];
    chern_numbers = new dtype[n_phi];
    chern_numbers_theta = new dtype[n_mesh *n_phi];

    // theta threads allocate
    thds_theta = new int*[num_threads];
    thds_ctheta = new int*[num_threads];
    for(int i = 0; i < num_threads; i++) {
        thds_theta[i] = new int[(n_mesh + 1) / num_threads + 1];
        thds_ctheta[i] = new int[n_mesh / num_threads + 1];
    }
    theta_len = new int[num_threads];
    ctheta_len = new int[num_threads];

    // divide works into num_threads of pieces
    memset(theta_len, 0, num_threads*sizeof(int));
    memset(ctheta_len, 0, num_threads*sizeof(int));
    for(int i = 0; i < (n_mesh + 1); i++)
        thds_theta[i% num_threads][theta_len[i%num_threads]++]=i;
    for(int i = 0; i < n_mesh; i++)
        thds_ctheta[i% num_threads][ctheta_len[i%num_threads]++]=i;
    pthread_t *peer_thds = new pthread_t[num_threads];
    struct struct_solve params[num_threads]; 
    struct struct_Chern Chern_params[num_threads];


    // initialize output data file
    ofstream fchern;
    fchern.open("energy_Chern_No.txt");

    double coeff_time,pot_time,diag_time,chern_time,total_time;
    auto _t1=std::chrono::high_resolution_clock::now();
    // initialize and store the disorder independent coefficients
    prepare_coeff(coeff_mn, coeff_m_theta, coeff_jm,n_phi,off_head,dim_m,dim_n,n_mesh, L1, L2);
    auto _t2=std::chrono::high_resolution_clock::now();
    coeff_time=std::chrono::duration_cast<chrono::microseconds>(_t2-_t1).count()/1.0E6;
    cerr<<std::left<<setw(40)<<"coefficients initialization time: "<<coeff_time<<" s\n";

    /*********************************** MAIN PROGRAM ************************************************************/
    for(int n = 0; n < n_sample; n++) {
        memset(energy_levels, 0, n_phi*sizeof(dtype));
        memset(chern_numbers, 0, n_phi*sizeof(dtype));
        _t1=std::chrono::high_resolution_clock::now();
	// initialize the disorder potential
        generate_disorder_potential(impurity_x,impurity_y,impurity_intensity,impurity_num,L1,L2,seed,n);
        // calculate the coefficient part from potential and transfer to the device
        get_potential_coeff(impurity_x, impurity_y, impurity_intensity, v_mn,coeff_mn, dim_m,dim_n,off_head,impurity_num);
        _t2=std::chrono::high_resolution_clock::now();
        pot_time=std::chrono::duration_cast<chrono::microseconds>(_t2-_t1).count()/1.0E6;
	if(n==0)
        cerr<<setw(40)<<"potential initialization time: "<<pot_time<<" s\n";
        _t1=std::chrono::high_resolution_clock::now();
        // solve the (theta_1=0)-line of wave functions
        int theta_1=0;
	// multi-threaded version
        for(int id = 0; id < num_threads; id++) {
	    params[id].theta_1 = theta_1;
	    params[id].n_phi = n_phi;
	    params[id].n_mesh = n_mesh;
	    params[id].dim_m = dim_m;
	    params[id].dim_n = dim_n;
	    params[id].off_head = off_head;
            params[id].theta_2 = thds_theta[id];
            params[id].theta_len = theta_len[id];
            params[id].wave_function = wave_functions;
            params[id].energy = energy_theta;
            params[id].v_mn = v_mn;
            params[id].coeff_jm = coeff_jm;
            params[id].coeff_m_theta = coeff_m_theta;
            pthread_create(&(peer_thds[id]), NULL, peer_solve_projected,  (void* )&params[id]);
        }
        // join all the threads
        for(int id = 0; id < num_threads; id++)
            pthread_join(peer_thds[id], NULL);
	
	//solve_projected(theta_1,n_mesh,n_phi,dim_m,dim_n,off_head,wave_functions,coeff_m_theta,coeff_jm,v_mn,energy_theta);

        _t2=std::chrono::high_resolution_clock::now();
              
        diag_time=std::chrono::duration_cast<chrono::microseconds>(_t2-_t1).count()/1.0E6;
	if(n==0)
	cerr<<setw(40)<<"diagonalization time/k-point: "<<diag_time/(n_mesh+1)<<" s\n";

        // average the energy over the first line of theta_2
        // energy will be divided by (n_mesh+1)*(n_mesh+1) finally
        for(int i=0; i<n_phi;i++){
           for(int j=0; j<=n_mesh;j++)
               energy_levels[i] += energy_theta[j* n_phi +i];
        }

        // solve theta_1-line of wfs from 1 to n_mesh, chern number calculations are also performed
        for(int theta_1=1;theta_1<=n_mesh;theta_1++){
         // solve for another line
	
        for(int id = 0; id < num_threads; id++) {
	    params[id].theta_1 = theta_1;
	    params[id].n_phi = n_phi;
	    params[id].n_mesh = n_mesh;
	    params[id].dim_m = dim_m;
	    params[id].dim_n = dim_n;
	    params[id].off_head = off_head;
            params[id].theta_2 = thds_theta[id];
            params[id].theta_len = theta_len[id];
            params[id].wave_function = wave_functions;
            params[id].energy = energy_theta;
            params[id].v_mn = v_mn;
            params[id].coeff_jm = coeff_jm;
            params[id].coeff_m_theta = coeff_m_theta;
            pthread_create(&(peer_thds[id]), NULL, peer_solve_projected,  (void* )&params[id]);
        }
        // join all the threads
        for(int id = 0; id < num_threads; id++)
            pthread_join(peer_thds[id], NULL);
	
	//solve_projected(theta_1,n_mesh,n_phi,dim_m,dim_n,off_head,wave_functions,coeff_m_theta,coeff_jm,v_mn,energy_theta);

        // energy will be divided by (n_mesh+1)*(n_mesh+1) finally
        for(int i=0; i<n_phi;i++){
           for(int j=0; j<=n_mesh;j++)
               energy_levels[i] += energy_theta[j* n_phi +i];
        }
        

        _t1=std::chrono::high_resolution_clock::now();
       
        for(int id = 0; id < num_threads; id++) {
	    Chern_params[id].n_phi = n_phi;
	    Chern_params[id].dim_vec = n_phi;
	    Chern_params[id].n_mesh = n_mesh;
	    Chern_params[id].theta_1 = theta_1;
            Chern_params[id].theta_2 = thds_ctheta[id];
            Chern_params[id].theta_len = ctheta_len[id];
            Chern_params[id].wave_function = wave_functions;
            Chern_params[id].chern_numbers_theta = chern_numbers_theta;
            pthread_create(&(peer_thds[id]), NULL, peer_cal_Chern,  (void*)&Chern_params[id]);
        }
        // join all the threads
        for(int id = 0; id < num_threads; id++)
            pthread_join(peer_thds[id], NULL);
        
        //cal_Chern(wave_functions, chern_numbers_theta,n_phi,n_phi,n_mesh,theta_1);
        _t2=std::chrono::high_resolution_clock::now();
        chern_time=std::chrono::duration_cast<chrono::microseconds>(_t2-_t1).count()/1.0E6;

	if(theta_1==1 && n==0){
	    total_time=coeff_time+(pot_time+(chern_time*n_mesh+diag_time*(n_mesh+1)))*n_sample;
            cerr<<setw(40)<<"Chern No. calculation time/k-point: "<<chern_time/(n_mesh+1)<<" s\n";
	    cerr<<setw(40)<<"estimated total time:"<<std::right<<setw(3)<<setfill('0')<<(int(total_time))/3600<<setw(1)<<":";
	    cerr<<setw(2)<<setfill('0')<<((int(total_time))%3600)/60<<setw(1)<<":";
	    cerr<<setw(2)<<setfill('0')<<(int(total_time))%60<<endl<<setfill(' ');
	    cerr<<std::left<<setw(40)<<"estimated # of samples/hour:"<<n_sample/(total_time/3600.0)<<endl;
	    cerr<<setw(40)<<"# module"<<"%(time)"<<endl;
	    cerr<<setw(40)<<"coefficients"<<coeff_time/total_time*100.0<<"\n";
	    cerr<<setw(40)<<"potential"<<pot_time*n_sample/total_time*100.0<<"\n";
	    cerr<<setw(40)<<"diagonalization"<<diag_time*(n_mesh+1)*n_sample/total_time*100.0<<endl;
	    cerr<<setw(40)<<"Chern No. calculation"<<chern_time*n_mesh*n_sample/total_time*100.0<<endl;

	}

        // collect the chern number contributions from the theta_1 line
        for(int k=0;k<n_phi;k++){
            for(int i=0;i<n_mesh;i++)
                chern_numbers[k]+=chern_numbers_theta[k*n_mesh+i];
          }
        }
        //average the energy levels by dividing (n_mesh+1)*(n_mesh+1) at the final stage
        for(int k=0;k<n_phi;k++){
             energy_levels[k]/=(n_mesh+1)*(n_mesh+1);
             chern_numbers[k]+=1.0/n_phi;
        }
        // print averaged energy levels
        for(int i = 0; i < n_phi; i++)
            fchern <<std::left<< setw(8) << n << " " << setw(15) << setprecision(8) << energy_levels[i] << " " << chern_numbers[i] << endl;
    }

    fchern.close();
    delete[] chern_numbers;
    delete[] energy_levels;
    delete[] wave_functions;
    delete[] v_mn;
    delete[] i_buffer;
    delete[] impurity_intensity;
    delete[] impurity_y;
    delete[] impurity_x;
    delete[] coeff_jm;
    delete[] coeff_m_theta;
    delete[] coeff_mn;

    return 0;
}

