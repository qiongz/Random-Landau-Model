#include"magma_diag.h"
void magma_diag_hamil(magmaFloatComplex* dev_wfs, float*energy,int n_phi){
    magma_init();
    magma_int_t dev=0;
    magma_queue_t queue=NULL;
    magma_queue_create_from_cuda(dev,NULL,NULL,NULL,&queue);
    magmaFloatComplex *h_R,*h_work,aux_work[1];

    float *rwork,aux_rwork[1];
    magma_int_t lrwork,info;
    float *w2;
    magma_int_t *iwork,aux_iwork[1];
    magma_int_t  lwork,liwork,lda,ldda;
    lda=n_phi;
    ldda=magma_roundup(n_phi,32);
    // query for workspace sizes
    magma_cheevd_gpu(MagmaVec,MagmaLower,
		    n_phi,NULL,lda,NULL,
		    NULL,lda,
		    aux_work, -1,
		    aux_rwork,-1,
		    aux_iwork,-1,
		    &info);

    lwork = (magma_int_t) MAGMA_C_REAL(aux_work[0]);
    lrwork = (magma_int_t) aux_rwork[0];
    liwork = aux_iwork[0];

    magma_smalloc_cpu(&w2, n_phi);
    magma_smalloc_cpu(&rwork, lrwork);
    magma_imalloc_cpu(&iwork, liwork);
    magma_cmalloc_pinned(&h_work,lwork);
    magma_cmalloc_pinned(&h_R,lwork);
    
    magma_cheevd_gpu(MagmaVec, MagmaLower,
	  	    n_phi,dev_wfs, ldda, energy,
		    h_R,lda,h_work,lwork,
		    rwork,lrwork,
		    iwork,liwork,
		    &info);

    magma_free_cpu(w2);
    magma_free_cpu(rwork);
    magma_free_cpu(iwork);
    magma_free_pinned(h_R);
    magma_free_pinned(h_work);
    magma_queue_destroy(queue);
    magma_finalize();
}
