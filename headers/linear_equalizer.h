#ifndef LINEAR_EQUALIZER_H_
#define LINEAR_EQUALIZER_H_

#include "signal.h"
#include "utility.h"

class LinearEqualizer {
    double estimated_angle;
    int angle_steps;

  public:
    LinearEqualizer();
    void train(const Signal& desired, const Signal& actual);
    Signal equalize(const Signal& original) const;

  private:
    Signal rotate(const Signal& original, const double& angle) const;
};

#endif  // LINEAR_EQUALIZER_H_
