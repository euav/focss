#ifndef PHASE_SHIFT_EQUALIZER_H_
#define PHASE_SHIFT_EQUALIZER_H_

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

#endif  // PHASE_SHIFT_EQUALIZER_H_
