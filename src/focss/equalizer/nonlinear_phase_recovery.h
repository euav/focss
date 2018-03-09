#ifndef FOCSS_EQUALIZER_NONLINEAR_PHASE_RECOVERY_H_
#define FOCSS_EQUALIZER_NONLINEAR_PHASE_RECOVERY_H_

#include "focss/field.h"

namespace focss {
class NonlinearPhaseRecovery {
  public:
    NonlinearPhaseRecovery();

    void train(const Field& desired, const Field& actual);
    Field equalize(const Field& original) const;
    RealVector get_weights() const;

  private:
    bool trained_;
    double a_, b_;
};
}  // namespace focss

#endif  // FOCSS_EQUALIZER_NONLINEAR_PHASE_RECOVERY_H_
