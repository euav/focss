#include "amplifier.h"
#include <cmath>
#include <random>
#include "focss/field.h"
#include "focss/functions.h"

namespace focss {
Amplifier::Amplifier() : gain(1), noise_factor(0) {}

Amplifier::Amplifier(const double& gain) : gain(gain), noise_factor(0) {}

Amplifier::Amplifier(const double& gain, const double& noise_factor)
    : gain(gain), noise_factor(noise_factor) {}

void Amplifier::amplify(Field& field) const {
    field *= std::sqrt(gain);

    if (noise_factor != 0) {
        double variance = focss::planck * field.center_frequency();
        variance *= noise_factor * (gain - 1) / 2;
        variance *= field.bandwidth();

        std::normal_distribution<double> awgn(0, std::sqrt(variance / 2));
        for (int i = 0; i < field.size(); ++i) 
            field[i] += complex_t(awgn(global_urng()), awgn(global_urng()));
    }
}

void Amplifier::give_power(Field& field) const {
    if (gain > 1)
        field *= std::sqrt(gain);
    else
        field *= std::sqrt(1 / gain);
}

void Amplifier::drop_power(Field& field) const {
    if (gain < 1)
        field *= std::sqrt(gain);
    else
        field *= std::sqrt(1 / gain);
}
}  // namespace focss
