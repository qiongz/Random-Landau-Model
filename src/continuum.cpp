// Continuum random-Landau model.
// Coefficient and Hamiltonian formulas are transcribed from the CPU sources
// on the legacy branch.
// Disorder positions and signs use std::mt19937, not MKL VSL.

#include "rlm/rlm.hpp"

#include <cmath>
#include <complex>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace rlm {
namespace {

// glibc M_PI, and PI2 = M_PI * 2 from matrix_coefficients.h.
constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = kPi * 2.0;

class ContinuumModel final : public Model {
public:
    explicit ContinuumModel(const ContinuumParams& params);

    std::string name() const override { return "continuum"; }
    int n_mesh() const override { return n_mesh_; }
    int n_state() const override { return n_phi_; }
    int vec_len() const override { return n_phi_; }
    int theta2_end() const override { return n_mesh_ + 1; }
    bool close_theta2_with_copy() const override { return false; }

    void assemble(int theta1, int theta2, cx* H) const override;
    void solve_point(int theta1, int theta2, cx* vectors, double* eval, Solver& solver) const override;
    void finalize(std::vector<double>& energy_sum, std::vector<double>& chern_sum) const override;

private:
    void check_theta(int theta1, int theta2) const;

    int n_phi_ = 0;
    int n_mesh_ = 0;
    int off_head_ = 0;
    int dim_m_ = 0;
    int dim_n_ = 0;
    std::vector<cx> coeff_mn_;
    std::vector<cx> coeff_m_theta_;
    std::vector<cx> coeff_jm_;
    std::vector<cx> vmn_;
};

ContinuumModel::ContinuumModel(const ContinuumParams& params)
    : n_phi_(params.n_phi), n_mesh_(params.n_mesh) {
    if (n_phi_ < 1) {
        throw std::runtime_error("continuum: n_phi < 1");
    }
    if (n_mesh_ < 1) {
        throw std::runtime_error("continuum: n_mesh < 1");
    }
    if (params.quanta_concentration <= 0.0) {
        throw std::runtime_error("continuum: quanta_concentration <= 0");
    }

    // Truncate toward 0, matching the historical assignment to a signed integer.
    const int impurity_num = static_cast<int>(params.impurity_concentration * n_phi_ / params.quanta_concentration);
    if (impurity_num < 1) {
        throw std::runtime_error("continuum: impurity_num < 1");
    }

    if (params.off_head >= 0) {
        off_head_ = params.off_head;
    } else {
        off_head_ = static_cast<int>(std::sqrt(static_cast<double>(impurity_num)) / 2.0 + 1.0);
    }
    if (off_head_ < 0) {
        throw std::runtime_error("continuum: off_head < 0");
    }
    if (off_head_ >= n_phi_) {
        throw std::runtime_error("continuum: off_head >= n_phi");
    }

    dim_m_ = off_head_ * 2 + 1;
    dim_n_ = dim_m_;

    const double L1 = std::sqrt(static_cast<double>(n_phi_) * params.quanta_concentration);
    const double L2 = std::sqrt(static_cast<double>(n_phi_) * params.quanta_concentration);

    coeff_mn_.resize(static_cast<std::size_t>(dim_m_) * static_cast<std::size_t>(dim_n_));
    coeff_m_theta_.resize(static_cast<std::size_t>(n_mesh_ + 1) * static_cast<std::size_t>(dim_m_));
    coeff_jm_.resize(static_cast<std::size_t>(n_phi_) * static_cast<std::size_t>(dim_m_));
    vmn_.resize(coeff_mn_.size());

    for (int m = 0; m < dim_m_; ++m) {
        for (int n = 0; n < dim_n_; ++n) {
            const double ang = kPi * (m - off_head_) * (n - off_head_) / n_phi_;
            const double env = std::exp(-kPi * ((m - off_head_) * (m - off_head_) * L2 / L1
                                                 + (n - off_head_) * (n - off_head_) * L1 / L2)
                                         / 2.0 / n_phi_);
            coeff_mn_[static_cast<std::size_t>(m * dim_n_ + n)] =
                cx(std::cos(ang), -std::sin(ang)) * env * (1.0 / n_phi_);
        }
    }

    for (int theta = 0; theta <= n_mesh_; ++theta) {
        const double dtheta = theta * kTwoPi / n_mesh_;
        for (int m = 0; m < dim_m_; ++m) {
            const double ang = (m - off_head_) * dtheta / n_phi_;
            coeff_m_theta_[static_cast<std::size_t>(theta * dim_m_ + m)] = cx(std::cos(ang), std::sin(ang));
        }
    }

    for (int j = 0; j < n_phi_; ++j) {
        for (int m = 0; m < dim_m_; ++m) {
            const double ang = kTwoPi * (m - off_head_) * (j + 1) / n_phi_;
            coeff_jm_[static_cast<std::size_t>(j * dim_m_ + m)] = cx(std::cos(ang), -std::sin(ang));
        }
    }

    // Sample k draws impurity_num x, then y, then signs. Earlier samples are
    // drawn and discarded so k does not depend on which samples were solved.
    std::mt19937 rng(params.seed);
    std::uniform_real_distribution<double> ux(0.0, L1);
    std::uniform_real_distribution<double> uy(0.0, L2);
    std::uniform_int_distribution<int> bit(0, 1);

    for (int sample = 0; sample < params.sample; ++sample) {
        for (int i = 0; i < impurity_num; ++i) {
            (void)ux(rng);
        }
        for (int i = 0; i < impurity_num; ++i) {
            (void)uy(rng);
        }
        for (int i = 0; i < impurity_num; ++i) {
            (void)bit(rng);
        }
    }

    std::vector<double> impurity_x(static_cast<std::size_t>(impurity_num));
    std::vector<double> impurity_y(static_cast<std::size_t>(impurity_num));
    std::vector<double> impurity_intensity(static_cast<std::size_t>(impurity_num));
    for (int i = 0; i < impurity_num; ++i) {
        impurity_x[static_cast<std::size_t>(i)] = ux(rng);
    }
    for (int i = 0; i < impurity_num; ++i) {
        impurity_y[static_cast<std::size_t>(i)] = uy(rng);
    }
    for (int i = 0; i < impurity_num; ++i) {
        impurity_intensity[static_cast<std::size_t>(i)] = (bit(rng) == 0) ? 1.0 : -1.0;
    }

    for (int m = 0; m < dim_m_; ++m) {
        const double kx = kTwoPi * (m - off_head_);
        for (int n = 0; n < dim_n_; ++n) {
            const double ky = kTwoPi * (n - off_head_);
            cx sum = 0.0;
            for (int i = 0; i < impurity_num; ++i) {
                const double phase = kx * impurity_x[static_cast<std::size_t>(i)]
                    + ky * impurity_y[static_cast<std::size_t>(i)];
                sum += cx(std::cos(phase), -std::sin(phase)) * impurity_intensity[static_cast<std::size_t>(i)];
            }
            const std::size_t idx = static_cast<std::size_t>(m * dim_n_ + n);
            vmn_[idx] = sum * coeff_mn_[idx];
        }
    }
}

void ContinuumModel::check_theta(int theta1, int theta2) const {
    if (theta1 < 0 || theta1 > n_mesh_ || theta2 < 0 || theta2 > n_mesh_) {
        throw std::runtime_error("continuum: theta out of range");
    }
}

void ContinuumModel::assemble(int theta1, int theta2, cx* H) const {
    check_theta(theta1, theta2);
    // Historical index arithmetic from hamiltonian_rlm_cpu.cpp, both triangles.
    // theta1 uses coeff_m_theta unconjugated. Do not swap the triangle LAPACK reads.
    const int idx_kx = theta1 * dim_n_;
    for (int j = 0; j < n_phi_; ++j) {
        for (int k = j; k < n_phi_; ++k) {
            H[j * n_phi_ + k] = 0.0;
            for (int m = 0; m < dim_m_; ++m) {
                const cx coeff_theta_jm =
                    std::conj(coeff_m_theta_[static_cast<std::size_t>(theta2 * dim_m_ + m)])
                    * coeff_jm_[static_cast<std::size_t>(j * dim_m_ + m)];
                cx sum = 0.0;
                if (k - j <= off_head_) {
                    const int n = off_head_ + k - j;
                    sum += vmn_[static_cast<std::size_t>(m * dim_n_ + n)]
                        * coeff_m_theta_[static_cast<std::size_t>(idx_kx + n)];
                }
                if ((n_phi_ + j - k) <= off_head_) {
                    const int n = off_head_ + k - j - n_phi_;
                    sum += vmn_[static_cast<std::size_t>(m * dim_n_ + n)]
                        * coeff_m_theta_[static_cast<std::size_t>(idx_kx + n)];
                }
                H[j * n_phi_ + k] += sum * coeff_theta_jm;
            }
            H[k * n_phi_ + j] = std::conj(H[j * n_phi_ + k]);
        }
    }
}

void ContinuumModel::solve_point(int theta1, int theta2, cx* vectors, double* eval, Solver& solver) const {
    check_theta(theta1, theta2);
    assemble(theta1, theta2, vectors);
    solver.diagonalize(vectors, eval, n_phi_);
}

void ContinuumModel::finalize(std::vector<double>& energy_sum, std::vector<double>& chern_sum) const {
    const double cells = static_cast<double>(n_mesh_ + 1) * static_cast<double>(n_mesh_ + 1);
    const double chern_shift = 1.0 / static_cast<double>(n_phi_);
    for (int k = 0; k < n_phi_; ++k) {
        energy_sum[static_cast<std::size_t>(k)] /= cells;
        chern_sum[static_cast<std::size_t>(k)] += chern_shift;
    }
}

}  // namespace

std::unique_ptr<Model> make_continuum(const ContinuumParams& params) {
    return std::make_unique<ContinuumModel>(params);
}

}  // namespace rlm
