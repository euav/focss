#ifndef LINEAR_EQUALIZER_H_
#define LINEAR_EQUALIZER_H_

#include "modulation.h"
#include "utility.h"

class LinearEqualizer {
    double estimated_angle;
    int angle_steps;

  public:
    LinearEqualizer();
    double train(const Signal& tx, const Signal& rx);
    Signal equalize(const Signal& original) const;

  private:
    Signal rotate(const Signal& original, const double& angle) const;
};

#endif  // LINEAR_EQUALIZER_H_
