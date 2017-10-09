#ifndef EDFA_H_
#define EDFA_H_

#include <random>
#include "focss/field.h"
#include "focss/utility.h"

struct EDFA {
    double gain = 1;
    double noise_factor = 0;
    double center_wavelength = 1.55e-6;
    double bandwidth;

    enum Scale { LINEAR, DECIBELS };

  public:
    EDFA(const double& gain, Scale = LINEAR);
    EDFA(const double& db_gain,
                  const double& noise_figure,
                  Scale = LINEAR);

    void amplify(Field& field) const;
    void drop_power(Field& field) const;
};

#endif  // EDFA_H_
