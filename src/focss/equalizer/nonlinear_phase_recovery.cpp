#include "focss/equalizer/nonlinear_phase_recovery.h"

const double stability_epsilon = 2.2204460492503131e-016;

NonlinearPhaseRecovery::NonlinearPhaseRecovery() : trained(false), radius(0) {}

NonlinearPhaseRecovery::NonlinearPhaseRecovery(
    const unsigned long& symbol_radius)
    : trained(false), radius(symbol_radius) {}

void NonlinearPhaseRecovery::setSymbolRadius(
    const unsigned long& symbol_radius) {
    radius = symbol_radius;
    trained = false;
}

void NonlinearPhaseRecovery::train(const Field& desired, const Field& actual) {
    if (desired.size() != actual.size() || actual.size() < 2 * radius + 1) {
        trained = false;
        return;
    }

    // double score, best_score, best_weight;
    // Field test = desired.chomp(radius, radius);

    // weights.assign(radius + 1, 0);
    // for (unsigned long current = 0; current <= radius; ++current) {
    //     weights[current] = -math_pi;
    //     best_weight = -math_pi;

    //     score = q2_factor(test, equalize(actual).chomp(radius, radius));
    //     best_score = score;

    //     for (int degree = -1440; degree < 1440; ++degree) {
    //         weights[current] = degree * math_pi / 5760;
    //         score = q2_factor(test, equalize(actual).chomp(radius, radius));

    //         if (best_score < score) {
    //             best_score = score;
    //             best_weight = weights[current];
    //         }
    //     }

    //     weights[current] = best_weight;
    // }

    weights.assign(2 * radius + 1, 0);
    RealVector powers(actual.size());
    for (unsigned long i = 0; i < actual.size(); ++i)
        powers[i] = norm(actual[i]);

    double initial_step = 0.16, degradation = 0.004;
    double residual, normalization, step;
    for (unsigned long i = radius; i < powers.size() - radius; ++i) {
        residual = arg(desired[i] / actual[i]);
        for (unsigned long j = 0; j <= 2 * radius; ++j)
            residual -= weights[j] * powers[i + j - radius];

        normalization = stability_epsilon;
        for (unsigned long j = 0; j <= 2 * radius; ++j)
            normalization += powers[i + j - radius] * powers[i + j - radius];

        step = initial_step / exp(i * degradation) * residual / normalization;
        for (unsigned long j = 0; j <= 2 * radius; ++j)
            weights[j] += step * powers[i + j - radius];
    }

    trained = true;
}

RealVector NonlinearPhaseRecovery::getWeights() const {
    if (trained)
        return weights;
    else
        return RealVector(1, 0);
}

Field NonlinearPhaseRecovery::equalize(const Field& original) const {
    if (!trained || original.size() < 2 * radius + 1) return original;

    Field equalized = original;
    for (unsigned long i = radius; i < original.size() - radius; ++i) {
        for (unsigned long j = 0; j <= 2 * radius; ++j) {
            equalized[i] *= i_exp(weights[j] * norm(original[i + j - radius]));
        }
    }

    return equalized;
}
