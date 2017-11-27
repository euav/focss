#include "phase_shift_equalizer.h"
#include "focss/functions.h"

namespace focss {
PhaseShiftEqualizer::PhaseShiftEqualizer()
    : angle_steps_(360), estimated_angle_(0) {}

PhaseShiftEqualizer::PhaseShiftEqualizer(const int& angle_steps)
    : angle_steps_(angle_steps), estimated_angle_(0) {}

void PhaseShiftEqualizer::set_angle_steps(const int& angle_steps) {
    angle_steps_ = angle_steps;
}

void PhaseShiftEqualizer::train(const ComplexVector& desired,
                                const ComplexVector& actual) {
    double angle, q2;
    double argmax_angle = 0;
    double max_q2 = q2_factor(desired, actual);
    for (int i = 0; i < angle_steps_; ++i) {
        angle = 2 * math_pi * double(i) / angle_steps_;
        q2 = q2_factor(desired, actual * i_exp(angle));

        if (max_q2 < q2) {
            max_q2 = q2;
            argmax_angle = angle;
        }
    }

    estimated_angle_ = argmax_angle;
}

void PhaseShiftEqualizer::train(const Field& desired, const Field& actual) {
    double angle, q2;
    double argmax_angle = 0;
    double max_q2 = q2_factor(desired, actual);
    for (int i = 0; i < angle_steps_; ++i) {
        angle = 2 * math_pi * double(i) / angle_steps_;
        q2 = q2_factor(desired, actual * i_exp(angle));

        if (max_q2 < q2) {
            max_q2 = q2;
            argmax_angle = angle;
        }
    }

    estimated_angle_ = argmax_angle;
}

double PhaseShiftEqualizer::get_angle() const { return estimated_angle_; }

Field PhaseShiftEqualizer::equalize(const Field& original) const {
    return original * i_exp(estimated_angle_);
}

ComplexVector PhaseShiftEqualizer::equalize(
    const ComplexVector& original) const {
    return original * i_exp(estimated_angle_);
}
}  // namespace focss
