#include "../headers/lms_equalizer.h"

LmsEqualizer::LmsEqualizer() : trained(false), length(0) {}

LmsEqualizer::LmsEqualizer(const unsigned long& filter_length)
    : trained(false), length(filter_length), weights(filter_length, 0) {}

void LmsEqualizer::setFilterLength(const unsigned long& filter_length) {
    length = filter_length;
    weights.assign(length, 0);
    trained = false;
}

void LmsEqualizer::train(const Field& desired, const Field& actual) {
    weights.assign(length, 0);

    if (desired.size() != actual.size() || actual.size() < length) return;

    mu = 0.16; alpha = 0.004;
    Complex residual, normalization, step;
    for (unsigned long i = length - 1; i < actual.size(); ++i) {
        residual = desired[i];
        for (unsigned long j = 0; j < length; ++j)
            residual -= weights[j] * actual[i - j];

        normalization = stability_epsilon;
        for (unsigned long j = 0; j < length; ++j)
            normalization += norm(actual[i - j]);

        step = mu / exp(i * alpha) * residual / normalization;
        for (unsigned long j = 0; j < length; ++j)
            weights[j] += step * conj(actual[i - j]);
    }

    trained = true;
}

Field LmsEqualizer::getWeights() const { return weights; }

Field LmsEqualizer::equalize(const Field& original) const {
    if (!trained || original.size() < length) return original;

    Field equalized(original.size(), 0);
    unsigned long cyclic_index;
    for (unsigned long i = 0; i < original.size(); ++i) {
        for (unsigned long j = 0; j < length; ++j) {
            cyclic_index = (i - j + original.size()) % original.size();
            equalized[i] += weights[j] * original[cyclic_index];
        }
    }

    return equalized;
}
