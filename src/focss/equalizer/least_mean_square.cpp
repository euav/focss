#include "focss/equalizer/least_mean_square.h"

const double stability_epsilon = 2.2204460492503131e-016;

LeastMeanSquare::LeastMeanSquare() : trained(false), radius(0) {}

LeastMeanSquare::LeastMeanSquare(const unsigned long& symbol_radius)
    : trained(false), radius(symbol_radius) {}

void LeastMeanSquare::setSymbolRadius(const unsigned long& symbol_radius) {
    radius = symbol_radius;
    trained = false;
}

void LeastMeanSquare::train(const Field& desired, const Field& actual) {
    if (desired.size() != actual.size() || actual.size() < 2 * radius + 1) {
        trained = false;
        return;
    }

    weights.assign(2 * radius + 1, 0);

    mu = 0.16;
    alpha = 0.004;
    Complex residual, normalization, step;
    for (unsigned long i = radius; i < actual.size() - radius; ++i) {
        residual = desired[i];
        for (unsigned long j = 0; j <= 2 * radius; ++j)
            residual -= weights[j] * actual[i + j - radius];

        normalization = stability_epsilon;
        for (unsigned long j = 0; j <= 2 * radius; ++j)
            normalization += norm(actual[i + j - radius]);

        step = mu / exp(i * alpha) * residual / normalization;
        for (unsigned long j = 0; j <= 2 * radius; ++j)
            weights[j] += step * conj(actual[i + j - radius]);
    }

    trained = true;
}

Field LeastMeanSquare::getWeights() const {
    if (trained)
        return weights;
    else
        return Field(1, 1);
}

Field LeastMeanSquare::equalize(const Field& original) const {
    if (!trained || original.size() < 2 * radius + 1) return original;

    Field equalized(original.size(), 0);
    for (unsigned long i = radius; i < original.size() - radius; ++i) {
        for (unsigned long j = 0; j <= 2 * radius; ++j) {
            equalized[i] += weights[j] * original[i + j - radius];
        }
    }

    return equalized;
}
