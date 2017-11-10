#ifndef FOCSS_EQUALIZER_PHASE_SHIFT_EQUALIZER_H_
#define FOCSS_EQUALIZER_PHASE_SHIFT_EQUALIZER_H_

#include "focss/functions.h"

namespace focss {
class PhaseShiftEqualizer {
    int steps;
    double estimated_angle;

  public:
    PhaseShiftEqualizer();
    PhaseShiftEqualizer(const int& angle_steps);
    void setAngleSteps(const int& angle_steps);

    void train(const ComplexVector& desired, const ComplexVector& actual);

    double getAngle() const;
    ComplexVector equalize(const ComplexVector& original) const;
};
}  // namespace focss
#endif  // FOCSS_EQUALIZER_PHASE_SHIFT_EQUALIZER_H_
