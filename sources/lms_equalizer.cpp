#include "../headers/lms_equalizer.h"

LmsEqualizer::LmsEqualizer() : order(1), weights(1, 0) {}

LmsEqualizer::LmsEqualizer(const int& filter_order)
    : order(filter_order), weights(order, 0) {}

void LmsEqualizer::train(const Signal& desired, const Signal& actual) {
    weights.assign(order, 0);

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
}
