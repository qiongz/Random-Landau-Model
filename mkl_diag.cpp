#include"mkl_diag.h"
void mkl_heevd(complex<dtype> *hamiltonian, dtype * energy, int l) {
    int  lda;
    char jobz, uplo;
    jobz = 'V';
    uplo = 'U';
    lda = l;
    HEEVD(LAPACK_COL_MAJOR, jobz, uplo,  l , hamiltonian, lda, energy );
}
