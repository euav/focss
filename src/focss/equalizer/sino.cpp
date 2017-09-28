#include "focss/equalizer/sino.h"

const double stability_epsilon = 2.2204460492503131e-016;

SINO::SINO() : trained(false), radius(0) {}

SINO::SINO(const unsigned long& symbol_radius)
    : trained(false), radius(symbol_radius) {}

void SINO::setSymbolRadius(const unsigned long& symbol_radius) {
    radius = symbol_radius;
    trained = false;
}

Complex SINO::term(const Field& field,
                   const unsigned long& center_index,
                   const unsigned long& term_index) const {
    unsigned long n = term_index % (2 * radius + 1);
    unsigned long m = term_index / (2 * radius + 1);

    Complex term = field[center_index + n - radius];
    term *= field[center_index + m - radius];
    term *= field[center_index + n + m - 2 * radius];

    return term;
}

void SINO::train(const Field& desired, const Field& actual) {
    if (desired.size() != actual.size() || actual.size() < 2 * radius + 1) {
        trained = false;
        return;
    }

    unsigned long size = (2 * radius + 1) * (2 * radius + 1);
    weights.assign(size, 0);

    double initial_step = 0.16, degradation = 0.004;
    Complex residual, normalization, step;
    for (unsigned long i = 2 * radius; i < actual.size() - 2 * radius; ++i) {
        residual = desired[i] - scalar_weight * actual[i];
        for (unsigned long j = 0; j < size; ++j) {
            residual -= weights[j] * term(actual, i, j);
        }

        normalization = stability_epsilon + norm(actual[i]);
        for (unsigned long j = 0; j < size; ++j) {
            normalization += norm(term(actual, i, j));
        }

        step = initial_step / exp(i * degradation) * residual / normalization;
        scalar_weight += step * conj(actual[i]);
        for (unsigned long j = 0; j < size; ++j) {
            weights[j] += step * conj(term(actual, i, j));
        }
    }

    trained = true;
}

Field SINO::getWeights() const {
    if (trained)
        return weights;
    else
        return Field(1, 0);
}

Field SINO::equalize(const Field& original) const {
    unsigned long size = (2 * radius + 1) * (2 * radius + 1);
    if (!trained || original.size() < size) return original;

    Field equalized = original;
    for (unsigned long i = 2 * radius; i < original.size() - 2 * radius; ++i) {
        equalized[i] *= scalar_weight;
        for (unsigned long j = 0; j < size; ++j) {
            equalized[i] += weights[j] * term(original, i, j);
        }
    }

    return equalized;
}
