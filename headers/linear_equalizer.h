#ifndef LINEAR_EQUALIZER_H_
#define LINEAR_EQUALIZER_H_

#include "field.h"
#include "utility.h"

class LinearEqualizer {
    double estimated_angle;
    int angle_steps;

  public:
    LinearEqualizer();
    void train(const Field& desired, const Field& actual);
    Field equalize(const Field& original) const;

  private:
    Field rotate(const Field& original, const double& angle) const;
};

#endif  // LINEAR_EQUALIZER_H_
