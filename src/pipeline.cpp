#include "rlm/rlm.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include <omp.h>

namespace rlm {
namespace {

void solve_theta2_point(const Model& model, Solver& solver, cx* wfs, double* eval, int theta1, int theta2,
                        int parity, int t2n, int ns, int vl) {
    const std::size_t block = static_cast<std::size_t>(ns) * static_cast<std::size_t>(vl);
    const std::size_t slot =
        (static_cast<std::size_t>(parity) * static_cast<std::size_t>(t2n) + static_cast<std::size_t>(theta2)) * block;
    model.solve_point(theta1, theta2, wfs + slot, eval + static_cast<std::size_t>(theta2) * static_cast<std::size_t>(ns),
                      solver);
}

}  // namespace

SampleResult run_sample(const Model& model, Solver& solver, int threads) {
    const int mesh = model.n_mesh();
    const int ns = model.n_state();
    const int vl = model.vec_len();
    const int t2n = mesh + 1;

    const std::size_t block = static_cast<std::size_t>(ns) * static_cast<std::size_t>(vl);
    std::vector<cx> wfs(static_cast<std::size_t>(2) * static_cast<std::size_t>(t2n) * block);
    std::vector<double> eval(static_cast<std::size_t>(t2n) * static_cast<std::size_t>(ns));
    std::vector<double> energy(static_cast<std::size_t>(ns), 0.0);
    std::vector<double> chern(static_cast<std::size_t>(ns), 0.0);
    std::vector<double> chern_theta(static_cast<std::size_t>(mesh) * static_cast<std::size_t>(ns));

    for (int theta1 = 0; theta1 <= mesh; ++theta1) {
        const int parity = theta1 % 2;
        const int t2_end = model.theta2_end();
        // Device solvers are not thread-safe. threads == 1 stays serial; threads <= 0
        // keeps the OpenMP default team size (no num_threads clause).
        if (solver.on_device() || threads == 1) {
            for (int theta2 = 0; theta2 < t2_end; ++theta2) {
                solve_theta2_point(model, solver, wfs.data(), eval.data(), theta1, theta2, parity, t2n, ns, vl);
            }
        } else if (threads <= 0) {
#pragma omp parallel for schedule(static)
            for (int theta2 = 0; theta2 < t2_end; ++theta2) {
                solve_theta2_point(model, solver, wfs.data(), eval.data(), theta1, theta2, parity, t2n, ns, vl);
            }
        } else {
#pragma omp parallel for schedule(static) num_threads(threads)
            for (int theta2 = 0; theta2 < t2_end; ++theta2) {
                solve_theta2_point(model, solver, wfs.data(), eval.data(), theta1, theta2, parity, t2n, ns, vl);
            }
        }

        if (model.close_theta2_with_copy()) {
            cx* const line = wfs.data() + static_cast<std::size_t>(parity) * static_cast<std::size_t>(t2n) * block;
            std::copy_n(line, block, line + static_cast<std::size_t>(mesh) * block);
            std::copy_n(eval.data(), static_cast<std::size_t>(ns),
                        eval.data() + static_cast<std::size_t>(mesh) * static_cast<std::size_t>(ns));
        }

        for (int theta2 = 0; theta2 <= mesh; ++theta2) {
            const double* ev = eval.data() + static_cast<std::size_t>(theta2) * static_cast<std::size_t>(ns);
            for (int state = 0; state < ns; ++state) {
                energy[static_cast<std::size_t>(state)] += ev[state];
            }
        }

        if (theta1 >= 1) {
            chern_strip(wfs.data(), chern_theta.data(), ns, vl, mesh, theta1);
            for (int state = 0; state < ns; ++state) {
                for (int m = 0; m < mesh; ++m) {
                    chern[static_cast<std::size_t>(state)] += chern_theta[m + state * mesh];
                }
            }
        }
    }

    model.finalize(energy, chern);
    return {energy, chern};
}

}  // namespace rlm
