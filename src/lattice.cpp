#include "rlm/rlm.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace rlm {
namespace {

void validate_params(const LatticeParams& params) {
    if (params.lx < 1 || params.ly < 1 || params.Q < 1 || params.n_mesh < 1 ||
        params.lx % params.Q != 0 || params.band < 0 || params.band >= params.Q) {
        throw std::runtime_error(
            "lattice: require lx, ly, Q, n_mesh >= 1, lx divisible by Q, and band in [0, Q)");
    }
    if (params.sample < 0) {
        throw std::runtime_error("lattice: sample must be >= 0");
    }
}

// sample k burns k * n_sites draws first, so it does not depend on earlier
// samples having been constructed. std::mt19937, not MKL VSL.
std::vector<double> draw_potential(const LatticeParams& params, int n_sites) {
    std::mt19937 rng(params.seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    const long long skip =
        static_cast<long long>(params.sample) * static_cast<long long>(n_sites);
    for (long long i = 0; i < skip; ++i) {
        (void)u(rng);
    }
    std::vector<double> potential(static_cast<std::size_t>(n_sites));
    for (int j = 0; j < n_sites; ++j) {
        potential[static_cast<std::size_t>(j)] = u(rng);
    }
    return potential;
}

struct ValidatedParams {
    LatticeParams params;
};

ValidatedParams vet_params(const LatticeParams& params) {
    validate_params(params);
    return ValidatedParams{params};
}

class LatticeModel final : public Model {
public:
    explicit LatticeModel(const LatticeParams& params)
        : LatticeModel(vet_params(params)) {}

private:
    explicit LatticeModel(ValidatedParams in)
        : lx_(in.params.lx),
          ly_(in.params.ly),
          Q_(in.params.Q),
          band_(in.params.band),
          n_mesh_(in.params.n_mesh),
          n_sites_(in.params.lx * in.params.ly),
          n_phi_(n_sites_ / in.params.Q),
          potential_(draw_potential(in.params, n_sites_)) {}

public:

    std::string name() const override { return "lattice"; }
    int n_mesh() const override { return n_mesh_; }
    int n_state() const override { return n_phi_; }
    int vec_len() const override { return n_sites_; }
    int theta2_end() const override { return n_mesh_; }
    bool close_theta2_with_copy() const override { return true; }

    void assemble(int theta1, int theta2, cx* H) const override {
        require_theta(theta1, theta2);
        // assemble has no Solver&. The Q×Q blocks still go through Solver
        // (fp64 zheevd). The projected matrix itself is not diagonalized here.
        std::unique_ptr<Solver> blocks = make_mkl_solver();
        std::vector<cx> psi(static_cast<std::size_t>(n_sites_) * static_cast<std::size_t>(n_phi_));
        clean_orbitals(theta1, theta2, psi.data(), *blocks);
        project(psi.data(), H);
    }

    void solve_point(int theta1, int theta2, cx* vectors, double* eval, Solver& solver) const override {
        require_theta(theta1, theta2);
        std::vector<cx> psi(static_cast<std::size_t>(n_sites_) * static_cast<std::size_t>(n_phi_));
        // Same Solver& for the Q×Q clean blocks and the projected matrix.
        clean_orbitals(theta1, theta2, psi.data(), solver);

        std::vector<cx> H(static_cast<std::size_t>(n_phi_) * static_cast<std::size_t>(n_phi_));
        project(psi.data(), H.data());
        solver.diagonalize(H.data(), eval, n_phi_);

        const std::size_t ns = static_cast<std::size_t>(n_sites_);
        const std::size_t np = static_cast<std::size_t>(n_phi_);
        for (int state = 0; state < n_phi_; ++state) {
            for (int site = 0; site < n_sites_; ++site) {
                cx sum = 0.0;
                for (int p = 0; p < n_phi_; ++p) {
                    sum += psi[static_cast<std::size_t>(site) + static_cast<std::size_t>(p) * ns] *
                           H[static_cast<std::size_t>(p) + static_cast<std::size_t>(state) * np];
                }
                vectors[static_cast<std::size_t>(site) + static_cast<std::size_t>(state) * ns] = sum;
            }
        }
    }

    void finalize(std::vector<double>& energy, std::vector<double>& chern) const override {
        if (energy.size() < static_cast<std::size_t>(n_phi_) ||
            chern.size() < static_cast<std::size_t>(n_phi_)) {
            throw std::runtime_error("lattice: finalize vectors are shorter than n_phi");
        }
        const double cells =
            static_cast<double>(n_mesh_ + 1) * static_cast<double>(n_mesh_ + 1);
        for (int k = 0; k < n_phi_; ++k) {
            energy[static_cast<std::size_t>(k)] /= cells;
            chern[static_cast<std::size_t>(k)] = -chern[static_cast<std::size_t>(k)];
        }
    }

private:
    void require_theta(int theta1, int theta2) const {
        if (theta1 < 0 || theta1 > n_mesh_ || theta2 < 0 || theta2 > n_mesh_) {
            throw std::runtime_error("lattice: theta is outside 0 .. n_mesh");
        }
    }

    // Clean subband orbitals at one twist. Transcription of get_wfs_ktheta
    // (the theta_1, theta_2 overload). psi is column-major, ld = n_sites.
    // Column-major Q×Q block, element (row, col) at row + col * Q:
    //   old hamiltonian[p*Q + q] is row q, column p.
    void clean_orbitals(int theta1, int theta2, cx* psi, Solver& solver) const {
        const double inv_mesh = 1.0 / static_cast<double>(n_mesh_);
        const double pi2 = 2.0 * std::acos(-1.0);
        const double pi2_lx = pi2 / static_cast<double>(lx_);
        const double pi2_ly = pi2 / static_cast<double>(ly_);
        const double inv_Q = 1.0 / static_cast<double>(Q_);

        std::vector<cx> h(static_cast<std::size_t>(Q_) * static_cast<std::size_t>(Q_));
        std::vector<double> energy(static_cast<std::size_t>(Q_));

        const int n_kx = lx_ / Q_;
        for (int kx = 0; kx < n_kx; ++kx) {
            const double kx_shifted = static_cast<double>(kx) + static_cast<double>(theta1) * inv_mesh;
            for (int ky = 0; ky < ly_; ++ky) {
                const double ky_shifted =
                    static_cast<double>(ky) + static_cast<double>(theta2) * inv_mesh;
                const int orbital = ky + kx * ly_;

                std::fill(h.begin(), h.end(), cx{0.0, 0.0});
                for (int p = 0; p < Q_; ++p) {
                    h[p + p * Q_] =
                        -2.0 * std::cos(ky_shifted * pi2_ly + pi2 * static_cast<double>(p) * inv_Q);
                    if (p + 1 < Q_) {
                        h[(p + 1) + p * Q_] = -1.0;
                        h[p + (p + 1) * Q_] = -1.0;
                    }
                }
                const double ang = kx_shifted * static_cast<double>(Q_) * pi2_lx;
                const cx corner(-std::cos(ang), -std::sin(ang));
                h[0 + (Q_ - 1) * Q_] = corner;
                h[(Q_ - 1) + 0 * Q_] = std::conj(corner);

                solver.diagonalize(h.data(), energy.data(), Q_);

                cx* col = psi + static_cast<std::size_t>(orbital) * static_cast<std::size_t>(n_sites_);
                for (int j = 0; j < ly_; ++j) {
                    for (int m = 0; m < lx_; m += Q_) {
                        const double phase = kx_shifted * pi2_lx * static_cast<double>(m) +
                                              ky_shifted * pi2_ly * static_cast<double>(j);
                        const cx phase_factor(std::cos(phase), -std::sin(phase));
                        for (int q = 0; q < Q_; ++q) {
                            const int site = j * lx_ + m + q;
                            col[site] = phase_factor * h[q + band_ * Q_];
                        }
                    }
                }

                double norm2 = 0.0;
                for (int site = 0; site < n_sites_; ++site) {
                    norm2 += std::norm(col[site]);
                }
                if (!(norm2 > 0.0)) {
                    throw std::runtime_error("lattice: clean orbital has zero norm");
                }
                const double inv_norm = 1.0 / std::sqrt(norm2);
                for (int site = 0; site < n_sites_; ++site) {
                    col[site] *= inv_norm;
                }
            }
        }
    }

    // H = Psi^H * diag(potential) * Psi, both triangles, column-major.
    void project(const cx* psi, cx* H) const {
        const std::size_t ns = static_cast<std::size_t>(n_sites_);
        const std::size_t np = static_cast<std::size_t>(n_phi_);
        for (int a = 0; a < n_phi_; ++a) {
            for (int b = a; b < n_phi_; ++b) {
                cx s = 0.0;
                for (int site = 0; site < n_sites_; ++site) {
                    s += std::conj(psi[static_cast<std::size_t>(site) + static_cast<std::size_t>(a) * ns]) *
                         potential_[static_cast<std::size_t>(site)] *
                         psi[static_cast<std::size_t>(site) + static_cast<std::size_t>(b) * ns];
                }
                if (a == b) {
                    s = cx(std::real(s), 0.0);
                }
                H[static_cast<std::size_t>(a) + static_cast<std::size_t>(b) * np] = s;
                H[static_cast<std::size_t>(b) + static_cast<std::size_t>(a) * np] = std::conj(s);
            }
        }
    }

    int lx_;
    int ly_;
    int Q_;
    int band_;
    int n_mesh_;
    int n_sites_;
    int n_phi_;
    std::vector<double> potential_;
};

}  // namespace

std::unique_ptr<Model> make_lattice(const LatticeParams& params) {
    return std::make_unique<LatticeModel>(params);
}

}  // namespace rlm
