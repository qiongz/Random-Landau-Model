#pragma once

// fp64 rewrite of the random Landau (continuum) and Hofstadter (lattice) models.
// One theta-sweep pipeline. The eigensolver is the only backend:
//   MKL LAPACKE_zheevd, or GPU rocSOLVER / cuSOLVER zheevd.
// Hamiltonian assembly stays on the host and is the CPU formula from the
// historical sources. Do not re-derive it.

#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace rlm {

using cx = std::complex<double>;

enum class Backend { Mkl, Gpu };

// Column-major Hermitian n×n matrix, both triangles filled, lda = n.
// On success eval[0..n) is ascending and column j of H is eigenvector j
// (LAPACK COL_MAJOR layout: element (row, col) at row + col * n).
struct Solver {
    virtual ~Solver() = default;
    virtual bool on_device() const { return false; }
    virtual void diagonalize(cx* H, double* eval, int n) = 0;
};

std::unique_ptr<Solver> make_mkl_solver();
std::unique_ptr<Solver> make_gpu_solver();
std::unique_ptr<Solver> make_solver(Backend backend);

struct Model {
    virtual ~Model() = default;
    virtual std::string name() const = 0;
    virtual int n_mesh() const = 0;
    virtual int n_state() const = 0;
    // Length of one stored eigenvector. Continuum: n_phi. Lattice: n_sites.
    virtual int vec_len() const = 0;
    // Exclusive end of the independent θ2 solves.
    // Continuum: n_mesh + 1 (θ2 = 0..n_mesh are all diagonalized).
    // Lattice:   n_mesh     (θ2 = n_mesh is a periodic copy of θ2 = 0).
    virtual int theta2_end() const = 0;
    virtual bool close_theta2_with_copy() const = 0;
    // Projected Hamiltonian only, n_state × n_state, column-major, both
    // triangles filled. Used by the residual check. Thread-safe.
    virtual void assemble(int theta1, int theta2, cx* H) const = 0;
    // Write eigenvectors into `vectors`, length n_state * vec_len.
    // State j starts at vectors + j * vec_len.
    // Write eigenvalues into eval[0..n_state), ascending.
    // Must be safe to call concurrently for different θ2 on the same model:
    // disorder and coefficients are read-only, scratch is local.
    virtual void solve_point(int theta1, int theta2, cx* vectors, double* eval, Solver& solver) const = 0;
    // energy_sum: sum of eigenvalues over θ1 = 0..n_mesh and θ2 = 0..n_mesh
    //   (θ2 = n_mesh included after the periodic copy, when that copy is used).
    // chern_sum: sum of chern_strip plaquette values. Each plaquette value is
    //   -arg(overlap product) / (2π), the historical peer_cal_Chern sign.
    // Overwrite both vectors with the reported per-state energy and Chern number.
    virtual void finalize(std::vector<double>& energy_sum, std::vector<double>& chern_sum) const = 0;
};

struct SampleResult {
    std::vector<double> energy;
    std::vector<double> chern;
};

// threads <= 0 selects the OpenMP default. Device solvers run the θ2 loop serially.
SampleResult run_sample(const Model& model, Solver& solver, int threads);

// wfs layout, length 2 * (n_mesh + 1) * n_state * vec_len:
//   ((theta1_parity) * (n_mesh + 1) + theta2) * n_state * vec_len + state * vec_len
// Writes n_mesh * n_state values:
//   chern_theta[theta2 + state * n_mesh] = -arg(prod) / (2π)
void chern_strip(const cx* wfs, double* chern_theta, int n_state, int vec_len, int n_mesh, int theta1);

struct ContinuumParams {
    int n_phi = 16;
    int n_mesh = 16;
    double quanta_concentration = 1.0;
    double impurity_concentration = 16.0;
    int off_head = -1;  // < 0: floor(sqrt(n_imp) / 2 + 1), matching the historical cast
    unsigned seed = 1;
    int sample = 0;
};

std::unique_ptr<Model> make_continuum(const ContinuumParams& params);

struct LatticeParams {
    int lx = 3;
    int ly = 3;
    int Q = 3;
    int band = 0;
    int n_mesh = 16;
    unsigned seed = 1;
    int sample = 0;
};

std::unique_ptr<Model> make_lattice(const LatticeParams& params);

}  // namespace rlm
