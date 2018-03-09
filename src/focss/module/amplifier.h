#ifndef FOCSS_MODULE_AMPLIFIER_H_
#define FOCSS_MODULE_AMPLIFIER_H_

#include "focss/field.h"

namespace focss {
struct Amplifier {
    double gain;
    double noise_factor;

  public:
    Amplifier();
    Amplifier(const double& gain);
    Amplifier(const double& gain, const double& noise_factor);

  public:
    void amplify(Field& field) const;
    void give_power(Field& field) const;
    void drop_power(Field& field) const;
};
}  // namespace focss

#endif  // FOCSS_MODULE_AMPLIFIER_H_
