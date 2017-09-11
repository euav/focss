#include "../headers/fiber.h"
#include <iostream>

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

void Fiber::setTotalSteps(const double& steps) { this->total_steps = steps; }

void Fiber::setSamplingRate(const double& rate) { this->sampling_rate = rate; }

Field Fiber::estimateLinearity(const int& samples) const {
    double freq;
    double freq_step = 2 * M_PI * sampling_rate / samples;
    double step = fiber_length / total_steps;

    Field linearity(samples, exp(-0.5 * alpha * step));
    for (int i = 0; i < samples; ++i) {
        if (i <= samples / 2)
            freq = freq_step * i;
        else
            freq = freq_step * (i - samples);

        linearity[i] *= i_exp(-0.5 * beta2 * freq * freq * step); 
    }

    return linearity;
}

void Fiber::propagate(Field& field) const {
    int samples = field.size();
    double step = fiber_length / total_steps;
    Field linearity = estimateLinearity(samples);

    for (int j = 0; j < samples; ++j)
        field[j] *= i_exp(0.5 * gamma * step * norm(field[j]));

    for (int i = 0; i < total_steps; ++i) {
        field.fft_inplace();
        field *= linearity;
        field.ifft_inplace();

        for (int j = 0; j < samples; ++j)
            field[j] *= i_exp(gamma * step * norm(field[j]));
    }

    for (int j = 0; j < samples; ++j)
        field[j] *= i_exp(-0.5 * gamma * step * norm(field[j]));
}

void Fiber::compensateCD(Field& field, const double& times) const {
    int samples = field.size();
    double freq;
    double freq_step = 2 * M_PI * sampling_rate / samples;

    field.fft_inplace();
    for (int i = 0; i < samples; ++i) {
        if (i <= samples / 2)
            freq = freq_step * i;
        else
            freq = freq_step * (i - samples);

        field[i] *= i_exp(0.5 * beta2 * freq * freq * times * fiber_length);
    }
    field.ifft_inplace();
}

void Fiber::amplify(Field& field) const {
    field *= exp(0.5 * alpha * fiber_length);
}
