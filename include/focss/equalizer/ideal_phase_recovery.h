#ifndef IDEAL_PHASE_RECOVERY_H_
#define IDEAL_PHASE_RECOVERY_H_

#include "focss/field.h"
#include "focss/utility.h"

class IdealPhaseRecovery {
    int steps;
    double estimated_angle;

  public:
    IdealPhaseRecovery();
    IdealPhaseRecovery(const int& angle_steps);
    void setAngleSteps(const int& angle_steps);

    void train(const Field& desired, const Field& actual);

    double getAngle() const;
    Field equalize(const Field& original) const;
};

#endif  // IDEAL_PHASE_RECOVERY_H_
