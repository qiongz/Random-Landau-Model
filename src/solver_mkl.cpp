#include "rlm/rlm.hpp"

#include <complex>
#include <memory>
#include <stdexcept>
#include <string>

#define MKL_Complex16 std::complex<double>
#include <mkl.h>

namespace rlm {
namespace {

class MklSolver final : public Solver {
public:
    bool on_device() const override { return false; }

    void diagonalize(cx* H, double* eval, int n) override {
        if (n <= 0) {
            throw std::runtime_error("diagonalize: n <= 0");
        }
        // One MKL thread per call. The theta2 sweep is the OpenMP parallelism,
        // and this local count is intentionally not restored.
        mkl_set_num_threads_local(1);
        const lapack_int n_lapack = static_cast<lapack_int>(n);
        const lapack_int info =
            LAPACKE_zheevd(LAPACK_COL_MAJOR, 'V', 'U', n_lapack, H, n_lapack, eval);
        if (info != 0) {
            throw std::runtime_error("LAPACKE_zheevd failed, info=" + std::to_string(info));
        }
    }
};

}  // namespace

std::unique_ptr<Solver> make_mkl_solver() {
    return std::make_unique<MklSolver>();
}

}  // namespace rlm
