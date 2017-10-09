#ifndef NONLINEAR_PHASE_RECOVERY_H_
#define NONLINEAR_PHASE_RECOVERY_H_

#include "focss/field.h"
#include "focss/utility.h"


class NonlinearPhaseRecovery {
    bool trained;
    double a, b;

  public:
    NonlinearPhaseRecovery();
    
    void train(const Field& desired, const Field& actual);
    
    RealVector getWeights() const;
    Field equalize(const Field& original) const;
};

#endif  // NONLINEAR_PHASE_RECOVERY_H_
