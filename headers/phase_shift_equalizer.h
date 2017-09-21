#ifndef PHASE_SHIFT_EQUALIZER_H_
#define PHASE_SHIFT_EQUALIZER_H_

#include "field.h"
#include "utility.h"

class PhaseShiftEqualizer {
    int steps;
    double estimated_angle;

  public:
    PhaseShiftEqualizer();
    PhaseShiftEqualizer(const int& angle_steps);
    void setAngleSteps(const int& angle_steps);

    void train(const Field& desired, const Field& actual);

    double getAngle() const;
    Field equalize(const Field& original) const;
};

#endif  // PHASE_SHIFT_EQUALIZER_H_
