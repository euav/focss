#include "focss/equalizer/sino.h"
#include "lasso_admm.h"
#include <armadillo>

SINO::SINO() : radius(0), lambda(1), trained(false) {}

SINO::SINO(const unsigned long& symbol_radius)
    : radius(symbol_radius), lambda(1), trained(false) {}

SINO::SINO(const unsigned long& symbol_radius, const double& parameter)
    : radius(symbol_radius), lambda(parameter), trained(false) {}

void SINO::setSymbolRadius(const unsigned long& symbol_radius) {
    radius = symbol_radius;
    trained = false;
}

void SINO::setRegularizationParameter(const double& parameter) {
    lambda = parameter;
    trained = false;
}

void SINO::train(const Field& desired, const Field& actual) {
    if (desired.size() != actual.size() || actual.size() < 4 * radius + 1)
        return;

    unsigned long n_samples = actual.size() - 4 * radius;
    unsigned long n_features = (2 * radius + 1) * (2 * radius + 1);
    arma::cx_mat features(n_samples, n_features + 1);
    arma::cx_vec responses = desired.chomp(2 * radius, 2 * radius);

    for (unsigned long i = 2 * radius; i < actual.size() - 2 * radius; ++i) {
        for (unsigned long j = 0; j < n_features; ++j) {
            features(i - 2 * radius, j) = term(actual, i, j);
        }
        features(i - 2 * radius, n_features) = actual[i];
    }

    arma::vec x_stddev = stddev(features, 1, 0).st();
    arma::cx_vec x_mean = mean(features, 0).st();
    Complex y_mean = mean(responses);
    responses = responses - y_mean;

    for (unsigned long i = 0; i < n_samples; ++i)
        for (unsigned long j = 0; j < n_features + 1; ++j)
            features(i, j) = (features(i, j) - x_mean(j)) / x_stddev(j);

    arma::cx_vec raw_weights = cx_lasso_admm(features, responses, lambda);
    raw_weights = raw_weights / x_stddev;
    intercept = y_mean - sum(raw_weights % x_mean / x_stddev);

    weights = arma::conv_to<ComplexVector>::from(raw_weights);
    trained = true;
}

Complex SINO::term(const Field& symbols,
                   const unsigned long& center_index,
                   const unsigned long& term_index) const {
    unsigned long n = term_index % (2 * radius + 1) - radius;
    unsigned long m = term_index / (2 * radius + 1) - radius;

    Complex term = symbols[center_index + n];
    term *= symbols[center_index + m];
    term *= conj(symbols[center_index + n + m]);

    return term;
}

ComplexVector SINO::getWeights() const {
    if (trained)
        return weights;
    else
        return ComplexVector(1, 1);
}

Field SINO::equalize(const Field& original) const {
    if (!trained || original.size() < 4 * radius + 1) return original;
    unsigned long n_features = (2 * radius + 1) * (2 * radius + 1);

    Field equalized(original.size());
    for (unsigned long i = 2 * radius; i < original.size() - 2 * radius; ++i) {
        equalized[i] = weights[n_features] * original[i];
        for (unsigned long j = 0; j < n_features; ++j)
            equalized[i] += weights[j] * term(original, i, j);
    }

    return equalized;
}
