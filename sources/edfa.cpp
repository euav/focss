#include "../headers/edfa.h"

EDFA::EDFA(const double& gain, Scale scale) {
    if (scale == LINEAR)
        this->gain = gain;
    else
        this->gain = db_to_linear(gain);
}

EDFA::EDFA(const double& gain, const double& noise_factor, Scale scale) {
    if (scale == LINEAR) {
        this->gain = gain;
        this->noise_factor = noise_factor;
    } else {
        this->gain = db_to_linear(gain);
        this->noise_factor = db_to_linear(noise_factor);
    }
}

void EDFA::amplify(Field& field) const {
    field *= std::sqrt(gain);

    double variance = 0.5 * noise_factor * (gain - 1);
    variance *= planck * light_speed / center_wavelength;
    variance *= field.getSamplingRate();

    std::mt19937 generator(time(0));
    std::normal_distribution<double> awgn(0, std::sqrt(variance / 2));
    for (unsigned long i = 0; i < field.size(); ++i)
        field[i] += Complex(awgn(generator), awgn(generator));
}

void EDFA::drop_power(Field& field) const { field *= std::sqrt(1 / gain); }