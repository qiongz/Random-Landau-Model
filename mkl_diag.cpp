#include"mkl_diag.h"
void mkl_cheevd(complex<float> *hamiltonian, float * energy, int l) {
    int  lda;
    char jobz, uplo;
    jobz = 'V';
    uplo = 'U';
    lda = l;
    LAPACKE_cheevd(LAPACK_COL_MAJOR, jobz, uplo,  l , hamiltonian, lda, energy );
}
