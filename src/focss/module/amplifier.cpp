#include "amplifier.h"
#include <cmath>
#include <random>
#include "focss/field.h"

namespace focss {
Amplifier::Amplifier() : gain_(1), noise_factor_(0) {}

Amplifier::Amplifier(const double& gain) : gain_(gain) {}

Amplifier::Amplifier(const double& gain, const double& noise_factor)
    : gain_(gain), noise_factor_(noise_factor) {}

void Amplifier::amplify(Field& field) const {
    field *= std::sqrt(gain_);

    if (noise_factor_ != 0) {
        double variance = focss::planck * field.get_center_frequency();
        variance *= noise_factor_ * (gain_ - 1) * field.get_sampling_rate(); 

        std::mt19937 generator(time(0));
        std::normal_distribution<double> awgn(0, std::sqrt(variance / 4));
        for (int i = 0; i < field.samples(); ++i) {
            field.x(i) += Complex(awgn(generator), awgn(generator));
            field.y(i) += Complex(awgn(generator), awgn(generator));
        }
    }
}

void Amplifier::give_gain(Field& field) const {
    if (gain_ > 1)
        field *= std::sqrt(gain_);
    else
        field *= std::sqrt(1 / gain_);
}

void Amplifier::drop_power(Field& field) const {
    if (gain_ < 1)
        field *= std::sqrt(gain_);
    else
        field *= std::sqrt(1 / gain_);
}
}  // namespace focss
