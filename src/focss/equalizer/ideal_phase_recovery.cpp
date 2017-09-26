#include "equalizer/ideal_phase_recovery.h"

IdealPhaseRecovery::IdealPhaseRecovery()
    : steps(360), estimated_angle(0) {}

IdealPhaseRecovery::IdealPhaseRecovery(const int& angle_steps)
    : steps(angle_steps), estimated_angle(0) {}

void IdealPhaseRecovery::setAngleSteps(const int& angle_steps) {
    steps = angle_steps;
}

void IdealPhaseRecovery::train(const Field& desired, const Field& actual) {
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

double IdealPhaseRecovery::getAngle() const { return estimated_angle; }

Field IdealPhaseRecovery::equalize(const Field& original) const {
    return original * i_exp(estimated_angle);
}
