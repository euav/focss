#include "phase_shift_equalizer.h"
#include "focss/functions.h"

namespace focss {
PhaseShiftEqualizer::PhaseShiftEqualizer() : steps(360), estimated_angle(0) {}

PhaseShiftEqualizer::PhaseShiftEqualizer(const int& angle_steps)
    : steps(angle_steps), estimated_angle(0) {}

void PhaseShiftEqualizer::setAngleSteps(const int& angle_steps) {
    steps = angle_steps;
}

void PhaseShiftEqualizer::train(const ComplexVector& desired,
                                const ComplexVector& actual) {
    double angle, q2;
    double argmax_angle = 0;
    double max_q2 = q2_factor(desired, actual);
    for (int i = 0; i < steps; ++i) {
        angle = 2.0 * M_PI * double(i) / steps;
        q2 = q2_factor(desired, actual * i_exp(angle));

        if (max_q2 < q2) {
            max_q2 = q2;
            argmax_angle = angle;
        }
    }

    estimated_angle = argmax_angle;
}

double PhaseShiftEqualizer::getAngle() const { return estimated_angle; }

ComplexVector PhaseShiftEqualizer::equalize(
    const ComplexVector& original) const {
    return original * i_exp(estimated_angle);
}
}  // namespace focss
