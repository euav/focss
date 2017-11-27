#ifndef FOCSS_EQUALIZER_PHASE_SHIFT_EQUALIZER_H_
#define FOCSS_EQUALIZER_PHASE_SHIFT_EQUALIZER_H_

#include "focss/functions.h"
#include "focss/field.h"

namespace focss {
class PhaseShiftEqualizer {
    int angle_steps_;
    double estimated_angle_;

  public:
    PhaseShiftEqualizer();
    PhaseShiftEqualizer(const int& angle_steps);
    void set_angle_steps(const int& angle_steps);

    void train(const Field& desired, const Field& actual);
    void train(const ComplexVector& desired, const ComplexVector& actual);

    double get_angle() const;
    Field equalize(const Field& original) const;
    ComplexVector equalize(const ComplexVector& original) const;
};
}  // namespace focss

#endif  // FOCSS_EQUALIZER_PHASE_SHIFT_EQUALIZER_H_
