#ifndef AMPLIFIER_H_
#define AMPLIFIER_H_

#include <random>
#include "focss/field.h"
#include "focss/utility.h"

struct Amplifier {
    double gain;
    double noise_factor;
    double center_wavelength = 1.55e-6;
    double bandwidth;

    enum Scale { LINEAR, DECIBELS };

  public:
    Amplifier();
    Amplifier(const double& gain, Scale = LINEAR);
    Amplifier(const double& db_gain,
                  const double& noise_figure,
                  Scale = LINEAR);

    void amplify(Field& field) const;
    void give_gain(Field& field) const;
    void drop_power(Field& field) const;
};

#endif  // EDFA_H_
