#include "../headers/fiber.h"

Fiber::Fiber(const double& center_wavelength) {
    wavelength = center_wavelength;

    alpha = 0;
    beta2 = 0;
    gamma = 0;
    fiber_length = 0;
    sampling_rate = 0;
    total_steps = 0;
}

void Fiber::setAttenuation(const double& alpha) { this->alpha = alpha; }

void Fiber::setAttenuationDB(const double& alpha_dB) {
    this->alpha = alpha_dB * std::log(10) / 10;
}

void Fiber::setDispersionPhysical(const double& dispersion) {
    double w_square = wavelength * wavelength;
    this->beta2 = -w_square * dispersion / (2 * M_PI * light_speed);
}

void Fiber::setDispersionEngineering(const double& beta2) {
    this->beta2 = beta2;
}

void Fiber::setNonlinearity(const double& gamma) { this->gamma = gamma; }

void Fiber::setFiberLength(const double& length) {
    this->fiber_length = length;
}

void Fiber::setTotalSteps(const int& steps) { this->total_steps = steps; }

void Fiber::setSamplingRate(const double& rate) { this->sampling_rate = rate; }

Signal Fiber::estimateLinearity(const int& samples) const {
    double freq;
    double freq_step = 2 * M_PI * sampling_rate / samples;
    double step = fiber_length / total_steps;

    Signal linearity(samples, 0);
    for (int i = 0; i < samples; ++i) {
        if (i <= samples / 2)
            freq = freq_step * i;
        else
            freq = freq_step * (i - samples);

        linearity[i] = exp(i_unit * beta2 / 2.0 * freq * freq * step);
        linearity[i] *= exp(-alpha * step / 2.0);
    }

    return linearity;
}

void Fiber::propagate(Signal& field) const {
    int samples = field.size();
    double step = fiber_length / total_steps;
    Signal linearity = estimateLinearity(samples);
    Complex nonlinearity = exp(i_unit * gamma * step);

    for (int j = 0; j < samples; ++j)
        field[j] *= pow(nonlinearity, norm(field[j]) / 2);

    for (int i = 0; i < total_steps; ++i) {
        field.fft_inplace();
        field *= linearity;
        field.ifft_inplace();

        for (int j = 0; j < samples; ++j)
            field[j] *= pow(nonlinearity, norm(field[j]));
    }

    for (int j = 0; j < samples; ++j)
        field[j] *= pow(nonlinearity, -norm(field[j]) / 2);
}

void Fiber::compensateCD(Signal& field, const double& times) const {
    int samples = field.size();
    double freq;
    double freq_step = 2 * M_PI * sampling_rate / samples;

    field.fft_inplace();
    for (int i = 0; i < samples; ++i) {
        if (i <= samples / 2)
            freq = freq_step * i;
        else
            freq = freq_step * (i - samples);

        field[i] *= exp(-i_unit * beta2 / 2.0 * freq * freq * times * fiber_length);
    }
    field.ifft_inplace();
}

void Fiber::amplify(Signal& field) const {
    field *= exp(alpha * fiber_length / 2.0);
}
