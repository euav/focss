#include "../headers/linear_equalizer.h"

LinearEqualizer::LinearEqualizer() {
    estimated_angle = 0;
    angle_steps = 360;
}

Field LinearEqualizer::equalize(const Field& original) const {
    return rotate(original, estimated_angle);
}

void LinearEqualizer::train(const Field& desired, const Field& actual) {
    double alpha, q2;
    double argmax_alpha = 0;
    double max_q2 = q2_factor(desired, actual);
    for (int i = 0; i < angle_steps; ++i) {
        alpha = 2 * M_PI * i / angle_steps;
        q2 = q2_factor(desired, rotate(actual, alpha));

        if (max_q2 < q2) {
            max_q2 = q2;
            argmax_alpha = alpha;
        }
    }

    estimated_angle = argmax_alpha;
}

Field LinearEqualizer::rotate(const Field& original,
                               const double& angle) const {
    return original * i_exp(angle);
}
