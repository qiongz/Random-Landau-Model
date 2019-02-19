#include"matrix_coefficients.h"

void prepare_coeff(complex<dtype> *coeff_mn, complex<dtype> *coeff_m_theta, complex<dtype> *coeff_jm, int n_phi, int off_head, int dim_m, int dim_n, int n_mesh, dtype L1, dtype L2) {
    int m, n, j;
    int theta;
    dtype dtheta;
    // coeff_mn = exp(-\pi*(m^2*L2/L1+n^2*L1/L2)/2N_phi) * exp(-i\pi*mn/N_phi)
    // divided by n_phi=L1*L2
    for(m = 0; m < dim_m; m++)
        for(n = 0; n < dim_n; n++)
            coeff_mn[m * dim_n + n] = complex<dtype>(cos(M_PI * (m - off_head) * (n - off_head) / n_phi), -sin(M_PI * (m - off_head) * (n - off_head) / n_phi)) * expf(-M_PI * ((m - off_head) * (m - off_head) * L2 / L1 + (n - off_head) * (n - off_head) * L1 / L2) / 2.0 / n_phi) *(1.0f/ n_phi);

    // coeff_m_theta =  exp(i\pi * m /N_phi)
    for(theta = 0; theta <= n_mesh; theta++) {
        dtheta = theta * PI2 / n_mesh;
        for(m = 0; m < dim_m; m++) {
            coeff_m_theta[theta * dim_m + m] = complex<dtype>(cos((m - off_head) * dtheta / n_phi), sin((m - off_head) * dtheta   / n_phi));
        }
    }
    // coeff_jm = exp(-i2\pi j*m/N_phi)
    for(j = 0; j < n_phi; j++) {
        for(m = 0; m < dim_m; m++)
            coeff_jm[j * dim_m + m] = complex<dtype>(cos(PI2 * (m - off_head) * (j + 1) / n_phi), -sin(PI2 * (m - off_head) * (j + 1) / n_phi));
    }
}

void get_potential_coeff(dtype* impurity_x, dtype * impurity_y, dtype *impurity_intensity, complex<dtype> * Vmn, complex<dtype>* coeff_mn, int dim_m, int dim_n, int off_head,int impurity_num) {
    int m, n, i;
    dtype kx, ky,phase ;
    // change impurity number to one
    for(m = 0; m < dim_m; m++) {
        kx = PI2 * (m - off_head);
        for(n = 0; n < dim_n; n++) {
            Vmn[m * dim_n + n] = 0;
            //  L1, L2 is not divided, since impurity_x and impurity_y in (0,1]
            ky = PI2 * (n - off_head);
            for(i = 0; i < impurity_num; i++){
		phase = kx * impurity_x[i] + ky * impurity_y[i];
                Vmn[m * dim_n + n] += complex<dtype>(cos(phase), -sin(phase)) * impurity_intensity[i];
		}
            Vmn[m*dim_n+n]*=coeff_mn[m*dim_n+n];
        }
    }
}
