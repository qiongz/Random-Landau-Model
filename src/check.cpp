#include "rlm/rlm.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct CheckLog {
    int fails = 0;

    void operator()(bool ok, const std::string& message) {
        std::cout << (ok ? "PASS " : "FAIL ") << message << '\n';
        if (!ok) {
            ++fails;
        }
    }
};

std::string sci(double value) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(6) << value;
    return out.str();
}

// ||H psi - E psi|| / ||E psi||, psi = column of the diagonalized matrix.
// H is the assembled matrix; diagonalize runs on a copy.
double max_relative_residual(const rlm::Model& model, rlm::Solver& solver, int theta1, int theta2) {
    const int n = model.n_state();
    if (n <= 0) {
        throw std::runtime_error("residual: n_state <= 0");
    }
    const std::size_t nn = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
    std::vector<rlm::cx> H(nn);
    model.assemble(theta1, theta2, H.data());
    std::vector<rlm::cx> work = H;
    std::vector<double> eval(static_cast<std::size_t>(n));
    solver.diagonalize(work.data(), eval.data(), n);

    double max_res = 0.0;
    for (int state = 0; state < n; ++state) {
        const rlm::cx* psi = work.data() + static_cast<std::size_t>(state) * static_cast<std::size_t>(n);
        const double energy = eval[static_cast<std::size_t>(state)];
        double num2 = 0.0;
        double den2 = 0.0;
        for (int row = 0; row < n; ++row) {
            rlm::cx hpsi = 0.0;
            for (int col = 0; col < n; ++col) {
                hpsi += H[static_cast<std::size_t>(row) + static_cast<std::size_t>(col) * static_cast<std::size_t>(n)] *
                        psi[col];
            }
            const rlm::cx diff = hpsi - energy * psi[row];
            num2 += std::norm(diff);
            den2 += std::norm(energy * psi[row]);
        }
        const double num = std::sqrt(num2);
        const double den = std::sqrt(den2);
        const double residual =
            den == 0.0 ? (num == 0.0 ? 0.0 : std::numeric_limits<double>::infinity()) : num / den;
        max_res = std::max(max_res, residual);
    }
    return max_res;
}

bool all_finite(const std::vector<double>& values) {
    if (values.empty()) {
        return false;
    }
    for (double value : values) {
        if (!std::isfinite(value)) {
            return false;
        }
    }
    return true;
}

void check_mkl_case(CheckLog& log, const std::string& tag, std::unique_ptr<rlm::Model> model) {
    auto solver = rlm::make_solver(rlm::Backend::Mkl);
    const double residual = max_relative_residual(*model, *solver, 0, 0);
    log(residual < 1e-8, tag + " residual theta=(0,0) max=" + sci(residual) + " < 1e-8");

    const rlm::SampleResult result = rlm::run_sample(*model, *solver, 4);
    const int ns = model->n_state();
    const bool energy_ok = static_cast<int>(result.energy.size()) == ns && all_finite(result.energy);
    log(energy_ok, tag + " energies finite");

    double sum = 0.0;
    const bool chern_ok = static_cast<int>(result.chern.size()) == ns && all_finite(result.chern);
    if (chern_ok) {
        for (double value : result.chern) {
            sum += value;
        }
    }
    const double err = std::abs(sum - 1.0);
    log(chern_ok && err < 1e-3, tag + " chern |sum-1|=" + sci(err) + " sum=" + sci(sum) + " < 1e-3");
}

rlm::ContinuumParams continuum_case() {
    rlm::ContinuumParams params;
    params.n_phi = 16;
    params.n_mesh = 16;
    params.quanta_concentration = 1.0;
    params.impurity_concentration = 16.0;
    params.off_head = -1;
    params.seed = 1;
    params.sample = 0;
    return params;
}

rlm::LatticeParams lattice_case() {
    rlm::LatticeParams params;
    params.lx = 3;
    params.ly = 3;
    params.Q = 3;
    params.band = 0;
    params.n_mesh = 16;
    params.seed = 1;
    params.sample = 0;
    return params;
}

#ifdef RLM_HAS_GPU
double max_abs_eig_diff(const rlm::Model& model, rlm::Solver& a, rlm::Solver& b, int theta1, int theta2) {
    const int n = model.n_state();
    if (n <= 0) {
        throw std::runtime_error("eigenvalue compare: n_state <= 0");
    }
    const std::size_t nn = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
    std::vector<rlm::cx> H(nn);
    model.assemble(theta1, theta2, H.data());
    std::vector<rlm::cx> Ha = H;
    std::vector<rlm::cx> Hb = H;
    std::vector<double> ea(static_cast<std::size_t>(n));
    std::vector<double> eb(static_cast<std::size_t>(n));
    a.diagonalize(Ha.data(), ea.data(), n);
    b.diagonalize(Hb.data(), eb.data(), n);
    double max_diff = 0.0;
    for (int i = 0; i < n; ++i) {
        max_diff = std::max(max_diff, std::abs(ea[static_cast<std::size_t>(i)] - eb[static_cast<std::size_t>(i)]));
    }
    return max_diff;
}

bool near_integer(double value) {
    return std::abs(value - std::round(value)) <= 0.05;
}

void check_gpu(CheckLog& log) {
    auto mkl = rlm::make_solver(rlm::Backend::Mkl);
    auto gpu = rlm::make_solver(rlm::Backend::Gpu);

    {
        auto model = rlm::make_continuum(continuum_case());
        const double diff = max_abs_eig_diff(*model, *mkl, *gpu, 0, 0);
        log(diff < 1e-8, "gpu/continuum eigenvalues theta=(0,0) max|dE|=" + sci(diff) + " < 1e-8");
    }
    {
        auto model = rlm::make_lattice(lattice_case());
        const double diff = max_abs_eig_diff(*model, *mkl, *gpu, 0, 0);
        log(diff < 1e-8, "gpu/lattice eigenvalues theta=(0,0) max|dE|=" + sci(diff) + " < 1e-8");
    }

    // mesh=4 keeps the paired run short. Per-state Chern on this mesh is not a fail
    // gate: individual plaquettes are noisy, and values near integers are reported only.
    rlm::ContinuumParams coarse = continuum_case();
    coarse.n_mesh = 4;
    auto model = rlm::make_continuum(coarse);
    const rlm::SampleResult mkl_result = rlm::run_sample(*model, *mkl, 4);
    const rlm::SampleResult gpu_result = rlm::run_sample(*model, *gpu, 4);
    const bool finite = all_finite(mkl_result.energy) && all_finite(mkl_result.chern) &&
                        all_finite(gpu_result.energy) && all_finite(gpu_result.chern) &&
                        mkl_result.energy.size() == gpu_result.energy.size() &&
                        mkl_result.chern.size() == gpu_result.chern.size() && !mkl_result.chern.empty();
    log(finite, "gpu/continuum mesh=4 run_sample finite");

    double sum_mkl = 0.0;
    double sum_gpu = 0.0;
    double max_dc = 0.0;
    bool abs_ok = finite;
    bool integers_ok = finite;
    if (finite) {
        for (std::size_t i = 0; i < mkl_result.chern.size(); ++i) {
            const double a = mkl_result.chern[i];
            const double b = gpu_result.chern[i];
            sum_mkl += a;
            sum_gpu += b;
            max_dc = std::max(max_dc, std::abs(a - b));
            if (!(std::abs(a - b) < 1e-6)) {
                abs_ok = false;
            }
            if (!(near_integer(a) && near_integer(b) && std::lround(a) == std::lround(b))) {
                integers_ok = false;
            }
        }
    }
    std::cout << "INFO gpu/continuum mesh=4 chern_sum_mkl=" << sci(sum_mkl) << " chern_sum_gpu=" << sci(sum_gpu)
              << " |sum_mkl-1|=" << sci(std::abs(sum_mkl - 1.0)) << " |sum_gpu-1|=" << sci(std::abs(sum_gpu - 1.0))
              << " max|dC|=" << sci(max_dc) << " abs<1e-6=" << (abs_ok ? "yes" : "no")
              << " integers=" << (integers_ok ? "yes" : "no") << " (mesh=4 Chern is not a fail gate)\n";
}
#endif

template <typename Fn>
void guard(CheckLog& log, const char* name, Fn&& fn) {
    try {
        fn();
    } catch (const std::exception& ex) {
        log(false, std::string(name) + " exception: " + ex.what());
    }
}

}  // namespace

int rlm_check() {
    CheckLog log;
    guard(log, "continuum", [&] { check_mkl_case(log, "continuum", rlm::make_continuum(continuum_case())); });
    guard(log, "lattice", [&] { check_mkl_case(log, "lattice", rlm::make_lattice(lattice_case())); });
#ifdef RLM_HAS_GPU
    guard(log, "gpu", [&] { check_gpu(log); });
#else
    std::cout << "INFO gpu cross-check skipped (rebuild with -DRLM_GPU=HIP or -DRLM_GPU=CUDA)\n";
#endif
    if (log.fails != 0) {
        std::cout << "FAIL " << log.fails << "\n";
        return 1;
    }
    std::cout << "PASS\n";
    return 0;
}
