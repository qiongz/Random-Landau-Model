#include "rlm/rlm.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace rlm {
namespace {

// MKL zdotc conjugates the first argument: sum_i conj(a[i]) * b[i].
cx overlap_zdotc(const cx* a, const cx* b, int vec_len) {
    cx sum = 0.0;
    for (int i = 0; i < vec_len; ++i) {
        sum += std::conj(a[i]) * b[i];
    }
    return sum;
}

const cx* wave_at(const cx* wfs, int ki, int kj, int state, int n_mesh_pbc, int n_state, int vec_len) {
    const std::size_t dim_wfs = static_cast<std::size_t>(n_state) * static_cast<std::size_t>(vec_len);
    const std::size_t offset =
        (static_cast<std::size_t>(ki) * static_cast<std::size_t>(n_mesh_pbc) + static_cast<std::size_t>(kj)) * dim_wfs +
        static_cast<std::size_t>(vec_len) * static_cast<std::size_t>(state);
    return wfs + offset;
}

}  // namespace

void chern_strip(const cx* wfs, double* chern_theta, int n_state, int vec_len, int n_mesh, int theta1) {
    const int n_mesh_pbc = n_mesh + 1;
    int dki[4];
    int dkj[4];
    dki[0] = -(dki[2] = (theta1 % 2 == 0 ? 1 : -1));
    dki[1] = dki[3] = dkj[0] = dkj[2] = 0;
    dkj[3] = -1;
    dkj[1] = 1;

    std::vector<cx> prod(static_cast<std::size_t>(n_state));
    for (int m = 0; m < n_mesh; ++m) {
        int lki = 1 - (theta1 % 2);
        int lkj = m;
        std::fill(prod.begin(), prod.end(), cx(1.0, 0.0));
        for (int l = 0; l < 4; ++l) {
            const int rki = lki + dki[l];
            const int rkj = lkj + dkj[l];
            for (int state = 0; state < n_state; ++state) {
                const cx* left = wave_at(wfs, lki, lkj, state, n_mesh_pbc, n_state, vec_len);
                const cx* right = wave_at(wfs, rki, rkj, state, n_mesh_pbc, n_state, vec_len);
                prod[static_cast<std::size_t>(state)] *= overlap_zdotc(left, right, vec_len);
            }
            lki = rki;
            lkj = rkj;
        }
        for (int state = 0; state < n_state; ++state) {
            chern_theta[m + state * n_mesh] = -std::arg(prod[static_cast<std::size_t>(state)]) / (2.0 * M_PI);
        }
    }
}

}  // namespace rlm
