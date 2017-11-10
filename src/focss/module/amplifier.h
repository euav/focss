#ifndef FOCSS_MODULE_AMPLIFIER_H_
#define FOCSS_MODULE_AMPLIFIER_H_

#include "focss/field.h"

namespace focss {
struct Amplifier {
    double gain_;
    double noise_factor_;
    double center_wavelength = 1.55e-6;

  public:
    Amplifier();
    Amplifier(const double& gain);
    Amplifier(const double& gain, const double& noise_factor);

  public:
    void amplify(Field& field) const;
    void give_gain(Field& field) const;
    void drop_power(Field& field) const;
};
}  // namespace focss

#endif  // FOCSS_MODULE_AMPLIFIER_H_
