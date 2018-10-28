#include "wdm_sino.h"
#include <cassert>
#include <complex>
#include <iostream>
#include "admm_regression.h"

namespace focss {
WdmPPE::WdmPPE() : c_radius_(0), d_radius_(0), trained_(false) {}

WdmPPE::WdmPPE(const int& radii)
    : c_radius_(radii), d_radius_(radii), trained_(false) {}

WdmPPE::WdmPPE(const int& c_radius, const int& d_radius)
    : c_radius_(c_radius), d_radius_(d_radius), trained_(false) {}

WdmPPE& WdmPPE::set_c_radius(const int& c_radius) {
    c_radius_ = c_radius;
    trained_ = false;
    return *this;
}

WdmPPE& WdmPPE::set_d_radius(const int& d_radius) {
    d_radius_ = d_radius;
    trained_ = false;
    return *this;
}

void WdmPPE::train(const arma::cx_mat& tx, const arma::cx_mat& rx) {
    using namespace arma;

    int n_samples = tx.n_cols;
    int n_feat_spm = (c_radius_ < 0) ? 0 : c_radius_ * (c_radius_ + 2) / 2 + 1;
    int n_feat_xpm = (d_radius_ < 0) ? 0 : 2 * d_radius_ * d_radius_ + 2 * d_radius_ + 1;
    int n_features = n_feat_spm + n_feat_xpm + 1;
    int radius = std::max(c_radius_, d_radius_);

    std::cout << "> preparing responses" << std::endl;
    cx_vec responses =
        vectorise(rx(span(2, 3), span(radius, n_samples - radius - 1)), 1).st();

    std::cout << "> preparing features" << std::endl;
    cx_cube feat_cube(n_samples - 2 * radius, n_features, 2, fill::zeros);
    for (int pol = 0; pol < 2; ++pol) {
        for (int sample = radius; sample < n_samples - radius; ++sample) {
            uword index = 0;

            int m_border = c_radius_ / 2;
            for (int m = -m_border; m <= m_border; ++m) {
                int n_left = std::abs(m);
                int n_right = c_radius_ - std::abs(m);

                for (int n = n_left; n <= n_right; ++n) {
                    if (m == 0 && n == 0) {
                        feat_cube(sample - radius, index++, pol) =
                            T(tx, pol, sample, 0, 0);
                    } else if (m == n) {
                        feat_cube(sample - radius, index++, pol) =
                            T(tx, pol, sample, m, n) +
                            T(tx, pol, sample, -m, -n);
                    } else if (m == -n) {
                        feat_cube(sample - radius, index++, pol) =
                            T(tx, pol, sample, m, n) + T(tx, pol, sample, n, m);
                    } else {
                        feat_cube(sample - radius, index++, pol) =
                            T(tx, pol, sample, m, n) +
                            T(tx, pol, sample, n, m) +
                            T(tx, pol, sample, -m, -n) +
                            T(tx, pol, sample, -n, -m);
                    }
                }
            }
            for (int m = -d_radius_; m <= d_radius_; ++m) {
                int n_border = d_radius_ - std::abs(m);
                for (int n = -n_border; n <= n_border; ++n) {
                    feat_cube(sample - radius, index++, pol) =
                        T1(tx, pol, sample, m, n) + T3(tx, pol, sample, m, -n);
                }
            }
            feat_cube(sample - radius, index, pol) = tx(2 + pol, sample);
        }
    }
    cx_mat features = join_vert(feat_cube.slice(0), feat_cube.slice(1));
    feat_cube.clear();

    std::cout << "> normalizing matrix" << std::endl;
    vec x_stddev = stddev(features, 1, 0).st();
    cx_vec x_mean = mean(features, 0).st();
    complex_t y_mean = mean(responses);
    responses = responses - y_mean;

    for (uword i = 0; i < features.n_rows; ++i)
        for (uword j = 0; j < features.n_cols; ++j)
            features(i, j) = (features(i, j) - x_mean(j)) / x_stddev(j);

    std::cout << "> solving system (" << features.n_rows;
    std::cout << ',' << features.n_cols << ')' << std::endl;
    cx_vec raw_weights = solve(features, responses);

    std::cout << "> denormalizing matrix" << std::endl;
    raw_weights = raw_weights / x_stddev;
    intercept_ = y_mean - sum(raw_weights % x_mean);
    C_ = cx_mat(2 * c_radius_ + 1, 2 * c_radius_ + 1, fill::zeros);
    D_ = cx_mat(2 * d_radius_ + 1, 2 * d_radius_ + 1, fill::zeros);

    uword index = 0;
    int m_border = c_radius_ / 2;
    for (int m = -m_border; m <= m_border; ++m) {
        int n_left = std::abs(m);
        int n_right = c_radius_ - std::abs(m);

        for (int n = n_left; n <= n_right; ++n) {
            if (m == 0 && n == 0) {
                C_(c_radius_, c_radius_) = raw_weights(index);
            } else if (m == n) {
                C_(c_radius_ + m, c_radius_ + n) = raw_weights(index);
                C_(c_radius_ - m, c_radius_ - n) = raw_weights(index);
            } else if (m == -n) {
                C_(c_radius_ + m, c_radius_ + n) = raw_weights(index);
                C_(c_radius_ + n, c_radius_ + m) = raw_weights(index);
            } else {
                C_(c_radius_ + m, c_radius_ + n) = raw_weights(index);
                C_(c_radius_ + n, c_radius_ + m) = raw_weights(index);
                C_(c_radius_ - m, c_radius_ - n) = raw_weights(index);
                C_(c_radius_ - n, c_radius_ - m) = raw_weights(index);
            }
            ++index;
        }
    }

    for (int m = -d_radius_; m <= d_radius_; ++m) {
        int n_border = d_radius_ - std::abs(m);
        for (int n = -n_border; n <= n_border; ++n) {
            D_(d_radius_ + m, d_radius_ + n) = raw_weights(index++);
        }
    }

    linear_weight_ = raw_weights(index);

    trained_ = true;
    std::cout << "> wdm-ppe has been trained (" << c_radius_;
    std::cout << ',' << d_radius_ << ')' << std::endl;
}

complex_t WdmPPE::TT(const arma::cx_mat& A,
                     const int& pol,
                     const int& k,
                     const int& m,
                     const int& n) {
    complex_t term = 0.0;
    term += A(2, k + n) * std::conj(A(2, k + m + n + 1));
    term += A(3, k + n) * std::conj(A(3, k + m + n + 1));
    return A(2 + pol, k + m) * term;
}

complex_t WdmPPE::T(const arma::cx_mat& A,
                    const int& pol,
                    const int& k,
                    const int& m,
                    const int& n) {
    complex_t term = 0.0;
    term += A(2, k + n) * std::conj(A(2, k + m + n));
    term += A(3, k + n) * std::conj(A(3, k + m + n));
    return A(2 + pol, k + m) * term;
}

complex_t WdmPPE::T1(const arma::cx_mat& A,
                     const int& pol,
                     const int& k,
                     const int& m,
                     const int& n) {
    complex_t term = 0.0;
    if (pol == 0) {
        term += 2.0 * A(2, k + m) * A(0, k + n) * std::conj(A(0, k + m + n));
        term += A(2, k + m) * A(1, k + n) * std::conj(A(1, k + m + n));
        term += A(3, k + m) * A(0, k + n) * std::conj(A(1, k + m + n));
    } else {
        term += 2.0 * A(3, k + m) * A(1, k + n) * std::conj(A(1, k + m + n));
        term += A(3, k + m) * A(0, k + n) * std::conj(A(0, k + m + n));
        term += A(2, k + m) * A(1, k + n) * std::conj(A(0, k + m + n));
    }
    return term;
}

complex_t WdmPPE::T3(const arma::cx_mat& A,
                     const int& pol,
                     const int& k,
                     const int& m,
                     const int& n) {
    complex_t term = 0.0;
    if (pol == 0) {
        term += 2.0 * A(2, k + m) * A(4, k + n) * std::conj(A(4, k + m + n));
        term += A(2, k + m) * A(5, k + n) * std::conj(A(5, k + m + n));
        term += A(3, k + m) * A(4, k + n) * std::conj(A(5, k + m + n));
    } else {
        term += 2.0 * A(3, k + m) * A(5, k + n) * std::conj(A(5, k + m + n));
        term += A(3, k + m) * A(4, k + n) * std::conj(A(4, k + m + n));
        term += A(2, k + m) * A(5, k + n) * std::conj(A(4, k + m + n));
    }
    return term;
}

complex_t WdmPPE::T13(const arma::cx_mat& A,
                      const int& pol,
                      const int& k,
                      const int& m,
                      const int& n) {
    complex_t term = 0.0;
    if (pol == 0) {
        term += 2.0 * A(0, k + m) * A(4, k + n) * std::conj(A(2, k + m + n));
        term += A(1, k + m) * A(4, k + n) * std::conj(A(3, k + m + n));
        term += A(0, k + m) * A(5, k + n) * std::conj(A(3, k + m + n));
    } else {
        term += 2.0 * A(1, k + m) * A(5, k + n) * std::conj(A(3, k + m + n));
        term += A(0, k + m) * A(5, k + n) * std::conj(A(2, k + m + n));
        term += A(1, k + m) * A(4, k + n) * std::conj(A(2, k + m + n));
    }
    return term;
}

complex_t WdmPPE::get_intercept() const {
    assert(trained_);
    return intercept_;
}

complex_t WdmPPE::get_linear_weight() const {
    assert(trained_);
    return linear_weight_;
}

complex_t WdmPPE::get_coo_weight() const {
    assert(trained_);
    return C_(c_radius_, c_radius_);
}

arma::cx_mat WdmPPE::model(const arma::cx_mat& field) const {
    assert(trained_);
    return linear(field) + cubic(field);
}

arma::cx_mat WdmPPE::equalize(const arma::cx_mat& tx,
                              const arma::cx_mat& rx) const {
    return (rx.rows(2, 3) - cubic(tx)) / linear_weight_;
}

arma::cx_mat WdmPPE::equalize(const arma::cx_mat& tx_spm,
                              const arma::cx_mat& tx_xpm,
                              const arma::cx_mat& rx) const {
    return (rx.rows(2, 3) - spm(tx_spm) - xpm(tx_xpm)) / linear_weight_;
}

arma::cx_mat WdmPPE::linear(const arma::cx_mat& original) const {
    assert(trained_);
    return linear_weight_ * original.rows(2, 3);
}

arma::cx_mat WdmPPE::spm(const arma::cx_mat& field) const {
    assert(trained_);

    std::cout << "> spm estimation" << std::endl;
    arma::cx_mat noise(2, field.n_cols, arma::fill::zeros);
    for (int k = c_radius_; k < noise.n_cols - c_radius_; ++k) {
        for (int pol = 0; pol < 2; ++pol) {
            for (int n = -c_radius_; n <= c_radius_; ++n) {
                int borders = c_radius_ - std::abs(n);
                for (int m = -borders; m <= borders; ++m) {
                    noise(pol, k) += C_(c_radius_ + m, c_radius_ + n) *
                                     T(field, pol, k, m, n);
                }
            }
        }
    }

    return noise;
}

arma::cx_mat WdmPPE::xpm(const arma::cx_mat& field) const {
    assert(trained_);

    std::cout << "> xpm estimation" << std::endl;
    arma::cx_mat noise(2, field.n_cols, arma::fill::zeros);
    for (int k = d_radius_; k < noise.n_cols - d_radius_; ++k) {
        for (int pol = 0; pol < 2; ++pol) {
            for (int m = -d_radius_; m <= d_radius_; ++m) {
                int n_border = d_radius_ - std::abs(m);
                for (int n = -n_border; n <= n_border; ++n) {
                    noise(pol, k) +=
                        D_(d_radius_ + m, d_radius_ + n) *
                        (T1(field, pol, k, m, n) + T3(field, pol, k, m, -n));
                }
            }
        }
    }

    return noise;
}

arma::cx_mat WdmPPE::cubic(const arma::cx_mat& field) const {
    assert(trained_);
    return spm(field) + xpm(field);
}
}  // namespace focss
