#include"chern.h"
#include"wfs_file.h"
void cal_Chern( complex<dtype>* wfs_full,dtype* chern_numbers_theta,int n_phi,int dim_vec , int n_mesh, int theta_1) {
    long lki, rki, lkj, rkj, n, m, l;
    unsigned long dim_wfs=n_phi*dim_vec;
    unsigned long wfs_offset_l,wfs_offset_r;
    int inca, incb;
    inca = incb = 1;
    vector< complex<dtype> > prod;
    long dki[4], dkj[4];
    dki[0] = -(dki[2] = (theta_1 % 2 == 0 ? 1 : -1));
    dki[1] = dki[3] = dkj[0] = dkj[2] = 0;
    dkj[1] = -(dkj[3] = -1);
    for(m = 0 ; m < n_mesh; m++) {
        lki = 1 - theta_1 % 2;
        lkj = m;
        complex<dtype> sum=0;
        prod.assign(n_phi,complex<dtype>(1,0));
        for(l = 0; l < 4; l++) {
            rki = lki + dki[l];
            rkj = lkj + dkj[l];
            for(n=0; n<n_phi; n++) {
                wfs_offset_l=(lki*(n_mesh+1)+lkj)*dim_wfs+dim_vec*n;
                wfs_offset_r=(rki*(n_mesh+1)+rkj)*dim_wfs+dim_vec*n;
                DOT(&sum, &dim_vec, wfs_full + wfs_offset_l, &inca, wfs_full+wfs_offset_r, &incb );
                prod[n] *= sum;
            }
            lki = rki;
            lkj = rkj;
        }
        for(n=0; n<n_phi; n++)
            chern_numbers_theta[m+n*n_mesh ]=-arg(prod[n])/(M_PI*2.0);
    }
    prod.clear();
}

void *peer_cal_Chern( void *peer_Chern_params) {
    struct struct_Chern *params=(struct struct_Chern*) peer_Chern_params;
    long lki, rki, lkj, rkj, n, m, l,s, n_mesh_pbc;
    unsigned long wfs_offset_l,wfs_offset_r;
    unsigned long nxn_phi;
    unsigned long dim_wfs=params->n_phi*params->dim_vec;
    int inca, incb;
    inca = incb = 1;
    complex<dtype> sum;
    vector< complex<dtype> > prod;
    long dki[4], dkj[4];
    dki[0] = -(dki[2] = (params->theta_1 % 2 == 0 ? 1 : -1));
    dki[1] = dki[3] = dkj[0] = dkj[2] = 0;
    dkj[1] = -(dkj[3] = -1);

    n_mesh_pbc = params->n_mesh + 1;
    for(m = 0 ; m < params->theta_len; m++) {
        lki = 1 - params->theta_1 %2;
        lkj = params->theta_2[m];
	prod.assign(params->n_phi,1);
        for(l = 0; l < 4; l++) {
            rki = lki + dki[l];
            rkj = lkj + dkj[l];
   	    for(n=0;n<params->n_phi;n++){
               wfs_offset_l=(lki*n_mesh_pbc+lkj)*dim_wfs+params->dim_vec*n;
               wfs_offset_r=(rki*n_mesh_pbc+rkj)*dim_wfs+params->dim_vec*n;
               DOT(&sum, &(params->dim_vec), params->wave_function + wfs_offset_l, &inca, params->wave_function+wfs_offset_r, &incb );
	       //for(s=0;s<params->dim_vec;s++)
		  //sum+=conj(params->wave_function[wfs_offset_l+s])*params->wave_function[wfs_offset_r+s];
              prod[n] *= sum;
  	     }
            lki = rki;
            lkj = rkj;
          }
	for(n=0;n<params->n_phi;n++)
          params-> chern_numbers_theta[params->theta_2[m]+n*params->n_mesh]=-arg(prod[n])/(M_PI*2.0);
    }
    pthread_exit((void*)params);
}

void cal_Chern_wfs_IO(dtype* chern_numbers_theta, int n_phi,int dim_vec,int n_mesh, int theta_1) {
    int lki, rki, lkj, rkj, n, m, l;
    int inca, incb;
    inca = incb = 1;
    long dim_wfs=n_phi*dim_vec;
    complex<dtype>* wfs=new complex<dtype>[4*dim_wfs];
    vector< complex<dtype> > prod;
    complex<dtype> sum;
    int dki[4], dkj[4];
    dki[0] = -(dki[2] = (theta_1 % 2 == 0 ? 1 : -1));
    dki[1] = dki[3] = dkj[0] = dkj[2] = 0;
    dkj[1] = -(dkj[3] = -1);
    for(m = 0 ; m < n_mesh; m++) {
        rki=lki = 1 - theta_1 % 2;
        rkj=lkj = m;
	if(m==0)
        for(l=0;l<4;l++){
         read_wfs(rki,rkj,wfs+l*dim_wfs,dim_wfs);
         rki += dki[l];
         rkj += dkj[l];
        }
	else{
         rki += dki[0]+dki[1];
         rkj += dkj[0]+dkj[1];
	 for(n=0;n<dim_wfs;n++){
	   wfs[n]=wfs[n+3*dim_wfs];
	   wfs[n+dim_wfs]=wfs[n+2*dim_wfs];
	 } 
	 for(l=2;l<4;l++)  {
           read_wfs(rki,rkj,wfs+l*dim_wfs,dim_wfs);
           rki += dki[l];
           rkj += dkj[l];
	 }
	}
        prod.assign(n_phi,1);
        for(l = 0; l < 4; l++) {
            for(n=0; n<n_phi; n++) {
                DOT(&sum, &dim_vec, wfs+((l+1)%4)*dim_wfs + n*dim_vec, &inca, wfs+l*dim_wfs + n*dim_vec, &incb );
                prod[n] *= sum;
            }
        }
        for(n=0; n<n_phi; n++)
            chern_numbers_theta[m+n*n_mesh ]=arg(prod[n])/(M_PI*2.0);
    }
    delete[] wfs;
    prod.clear();
}
