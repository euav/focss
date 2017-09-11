#include "../headers/lms_equalizer.h"

LmsEqualizer::LmsEqualizer() : trained_flag(false), order(0) {}

LmsEqualizer::LmsEqualizer(const int& filter_order)
    : trained_flag(false), order(filter_order), weights(filter_order, 0) {}

void LmsEqualizer::setFilterOrder(const int& filter_order) {
    order = filter_order;
    weights.assign(order, 0);
    trained_flag = false;
}

void LmsEqualizer::train(const Field& desired, const Field& actual) {
    weights.assign(order, 0);

    if (desired.size() != actual.size() || actual.size() < order) return;

    Complex residual, normalization, step;
    for (int i = order - 1; i < actual.size(); ++i) {
        residual = desired[i];
        for (int j = 0; j < order; ++j)
            residual -= weights[j] * actual[i - j];

        normalization = stability_epsilon;
        for (int j = 0; j < order; ++j)
            normalization += norm(actual[i - j]);

        step = residual / normalization;
        for (int j = 0; j < order; ++j)
            weights[j] += step * conj(actual[i - j]);
    }

    trained_flag = true;
}

Field LmsEqualizer::getWeights() const { return weights; }

Field LmsEqualizer::equalize(const Field& original) const {
    if (!trained_flag || original.size() < order) return original;

    return convolution(original, weights).chomp(order - 1, 0);
}
