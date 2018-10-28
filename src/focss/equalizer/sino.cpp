#include "sino.h"
#include <cassert>
#include <iostream>
#include "admm_regression.h"

namespace focss {
SINO::SINO() : radius_(0), lambda_(1), type_(ORDINARY), trained_(false) {}

SINO::SINO(const int& radius)
    : radius_(radius), lambda_(1), type_(ORDINARY), trained_(false) {}

SINO::SINO(const int& radius, const double& parameter)
    : radius_(radius), lambda_(parameter), type_(ORDINARY), trained_(false) {}

SINO::SINO(const int& radius, const double& parameter, const SinoType& type)
    : radius_(radius), lambda_(parameter), type_(type), trained_(false) {}

SINO& SINO::set_radius(const int& radius) {
    radius_ = radius;
    trained_ = false;
    return *this;
}

SINO& SINO::set_lambda(const double& lambda) {
    lambda_ = lambda;
    trained_ = false;
    return *this;
}

SINO& SINO::set_type(const SinoType& type) {
    type_ = type;
    trained_ = false;
    return *this;
}

SINO& SINO::set_linear_weight(const complex_t& linear_weight) {
    linear_weight_ = linear_weight;
    trained_ = false;
    return *this;
}

void SINO::train(const Field& transmitted, const Field& received) {
    assert(transmitted.domain() == received.domain());
    using namespace arma;

    int n_samples_per_mode = transmitted.samples() - 2 * radius_;
    int n_modes = transmitted.modes();
    int n_samples = n_samples_per_mode * n_modes;
    int n_features = radius_ * (radius_ + 2) / 2 + 1;
    if (type_ == ORDINARY) n_features++;
    int shift = radius_;

    cx_vec responses(n_samples);
    if (type_ == ORDINARY || type_ == HOMOGENEOUS) {
        for (int mode = 0; mode < n_modes; ++mode)
            for (int sample = 0; sample < n_samples_per_mode; ++sample)
                responses(sample + mode * n_samples_per_mode) =
                    received(mode, sample + shift);
    } else {  // type_ == FIXED
        for (int mode = 0; mode < n_modes; ++mode)
            for (int sample = 0; sample < n_samples_per_mode; ++sample)
                responses(sample + mode * n_samples_per_mode) =
                    received(mode, sample + shift) -
                    transmitted(mode, sample + shift) * linear_weight_;
    }

    cx_mat features(n_samples, n_features, fill::zeros);
    for (int mode = 0; mode < n_modes; ++mode) {
        for (int sample = 0; sample < n_samples_per_mode; ++sample) {
            int index = 0;

            int m_border = radius_ / 2;
            for (int m = -m_border; m <= m_border; ++m) {
                int n_left = std::abs(m);
                int n_right = radius_ - std::abs(m);

                for (int n = n_left; n <= n_right; ++n) {
                    if (m == 0 && n == 0) {
                        features(sample + mode * n_samples_per_mode, index++) =
                            term(transmitted, mode, sample + shift, 0, 0);
                    } else if (m == n) {
                        features(sample + mode * n_samples_per_mode, index++) =
                            term(transmitted, mode, sample + shift, m, n) +
                            term(transmitted, mode, sample + shift, -m, -n);
                    } else if (m == -n) {
                        features(sample + mode * n_samples_per_mode, index++) =
                            term(transmitted, mode, sample + shift, m, n) +
                            term(transmitted, mode, sample + shift, n, m);
                    } else {
                        features(sample + mode * n_samples_per_mode, index++) =
                            term(transmitted, mode, sample + shift, m, n) +
                            term(transmitted, mode, sample + shift, n, m) +
                            term(transmitted, mode, sample + shift, -m, -n) +
                            term(transmitted, mode, sample + shift, -n, -m);
                    }
                }
            }

            if (type_ == ORDINARY)
                features(sample + mode * n_samples_per_mode, index) =
                    transmitted(mode, sample + shift);
        }
    }

    vec x_stddev = stddev(features, 1, 0).st();
    cx_vec x_mean = mean(features, 0).st();
    complex_t y_mean = mean(responses);
    responses = responses - y_mean;

    for (uword i = 0; i < n_samples; ++i)
        for (uword j = 0; j < n_features; ++j)
            features(i, j) = (features(i, j) - x_mean(j)) / x_stddev(j);

    cx_vec raw_weights = solve(features, responses);
    // cx_vec raw_weights = cx_lad_admm(features, responses);
    // cx_vec raw_weights = cx_ridge_cholesky(features, responses, lambda_);
    // cx_vec raw_weights = cx_sag(features, responses, lambda_);
    // cx_vec raw_weights = cx_modulus_cholesky(features, responses,
    // x_stddev(n_features - 1));
    std::cout << arma::norm(features * raw_weights - responses) << std::endl;

    raw_weights = raw_weights / x_stddev;
    intercept_ = y_mean - sum(raw_weights % x_mean);
    weights_ = cx_mat(2 * radius_ + 1, 2 * radius_ + 1, fill::zeros);

    int index = 0;
    int m_border = radius_ / 2;
    for (int m = -m_border; m <= m_border; ++m) {
        int n_left = std::abs(m);
        int n_right = radius_ - std::abs(m);

        for (int n = n_left; n <= n_right; ++n) {
            if (m == 0 && n == 0) {
                weights_(radius_, radius_) = raw_weights(index++);
            } else if (m == n) {
                weights_(radius_ + m, radius_ + n) = raw_weights(index);
                weights_(radius_ - m, radius_ - n) = raw_weights(index++);
            } else if (m == -n) {
                weights_(radius_ + m, radius_ + n) = raw_weights(index);
                weights_(radius_ + n, radius_ + m) = raw_weights(index++);
            } else {
                weights_(radius_ + m, radius_ + n) = raw_weights(index);
                weights_(radius_ + n, radius_ + m) = raw_weights(index);
                weights_(radius_ - m, radius_ - n) = raw_weights(index);
                weights_(radius_ - n, radius_ - m) = raw_weights(index++);
            }
        }
    }

    if (type_ == ORDINARY)
        linear_weight_ = raw_weights(index);
    else if (type_ == HOMOGENEOUS)
        linear_weight_ = 0.0;

    trained_ = true;
    std::cout << "> sino has been trained (";
    std::cout << radius_ << ")" << std::endl;
    save_matrix("data/mmm.csv");
}

complex_t SINO::term(const Field& field,
                     const int& mode,
                     const int& k,
                     const int& m,
                     const int& n) {
    complex_t term = 0.0;
    for (int p = 0; p < field.modes(); ++p)
        term += field(p, k + n) * conj(field(p, k + m + n));

    return field(mode, k + m) * term;
}

complex_t SINO::get_intercept() const {
    assert(trained_);
    return intercept_;
}

complex_t SINO::get_linear_weight() const {
    assert(trained_);
    return linear_weight_;
}

complex_t SINO::get_coo_weight() const {
    assert(trained_);
    return weights_(radius_, radius_);
}

Field SINO::model(const Field& original) const {
    assert(trained_);
    return linear(original) + cubic(original);
}

Field SINO::equalize(const Field& rx, const Field& tx) const {
    assert(trained_);
    return (rx - cubic(tx)) / linear_weight_;
}

Field SINO::linear(const Field& original) const {
    assert(trained_);
    return linear_weight_ * original;
}

Field SINO::cubic(const Field& original) const {
    assert(trained_);

    Field equalized(original.domain());
    for (int mode = 0; mode < original.modes(); ++mode) {
        for (int k = radius_; k < original.samples() - radius_; ++k) {
            for (int m = -radius_; m <= radius_; ++m) {
                int borders = radius_ - std::abs(m);
                for (int n = -borders; n <= borders; ++n) {
                    equalized(mode, k) += weights_(m + radius_, n + radius_) *
                                          term(original, mode, k, m, n);
                }
            }
        }
    }

    return equalized;
}

void SINO::save_matrix(const std::string& filename) {
    std::ofstream fout(filename);
    save_matrix(fout);
    fout.close();
}

void SINO::save_matrix(std::ostream& output) {
    assert(trained_);

    for (int m = 0; m < 2 * radius_ + 1; ++m) {
        for (int n = 0; n < 2 * radius_ + 1; ++n) {
            output << weights_(m, n).real() << ' ';
        }
        output << '\n';
    }

    output << "\n\n";
    for (int m = 0; m < 2 * radius_ + 1; ++m) {
        for (int n = 0; n < 2 * radius_ + 1; ++n) {
            output << weights_(m, n).imag() << ' ';
        }
        output << '\n';
    }

    output << "\n\n";
    for (int m = 0; m < 2 * radius_ + 1; ++m) {
        for (int n = 0; n < 2 * radius_ + 1; ++n) {
            output << std::abs(weights_(m, n)) << ' ';
        }
        output << '\n';
    }

    output << "\n\n";
    for (int m = 0; m < 2 * radius_ + 1; ++m) {
        for (int n = 0; n < 2 * radius_ + 1; ++n) {
            output << std::arg(weights_(m, n)) << ' ';
        }
        output << '\n';
    }
}
}  // namespace focss
