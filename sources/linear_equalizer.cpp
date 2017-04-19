#include "../headers/linear_equalizer.h"

LinearEqualizer::LinearEqualizer() {
    estimated_angle = 0;
    angle_steps = 360;
}

Signal LinearEqualizer::equalize(const Signal &original) const {
    return rotate(original, estimated_angle);
}

double LinearEqualizer::train(const Signal &tx, const Signal &rx) {
    double alpha, q2;
    double argmax_alpha = 0;
    double max_q2 = q2_factor(tx, rx);
    for (int i = 0; i < angle_steps; ++i) {
        alpha = 2 * M_PI * i / angle_steps;
        q2 = q2_factor(tx, rotate(rx, alpha));

        if (max_q2 < q2) {
            max_q2 = q2;
            argmax_alpha = alpha;
        }
    }

    estimated_angle = argmax_alpha;
    return max_q2;
}

Signal LinearEqualizer::rotate(const Signal& original, const double& angle) const {
    Signal rotated = original;
    for (int i = 0; i < original.size(); ++i)
        rotated[i] *= std::exp(i_unit * angle);

    return rotated;
}
