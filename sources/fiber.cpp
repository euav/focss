#include "../headers/fiber.h"
#include <iostream>

Fiber::Fiber(const double& center_wavelength) {
    alpha = 0;
    beta2 = 0;
    gamma = 0;
    length = 0;
    wavelength = center_wavelength;

    uniform_grid = true;
    total_steps = 0;
    max_ps = 0;
    max_step = 0;
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

void Fiber::setFiberLength(const double& length) { this->length = length; }

void Fiber::setTotalSteps(const int& steps) {
    uniform_grid = true;
    total_steps = steps;
}

void Fiber::setMaximumShiftAndStep(const double& shift, const double& step) {
    uniform_grid = false;
    max_ps = shift;
    max_step = step;
}

void Fiber::linearStep(Field& field, const double& step) const {
    field *= exp(-0.5 * alpha * step);

    field.fft_inplace();
    for (int i = 0; i < field.size(); ++i)
        field[i] *= i_exp(-0.5 * beta2 * field.w(i) * field.w(i) * step);
    field.ifft_inplace();
}

void Fiber::nonlinearStep(Field& field, const double& step) const {
    for (int j = 0; j < field.size(); ++j)
        field[j] *= i_exp(gamma * step * norm(field[j]));
}

void Fiber::propagate(Field& field) const {
    if (uniform_grid) {
        double step = length / total_steps;

        linearStep(field, 0.5 * step);
        for (int i = 0; i < total_steps; ++i) {
            nonlinearStep(field, step);
            linearStep(field, step);
        }
        linearStep(field, -0.5 * step);
    } else {
        double distance = 0;
        double phase_shift = std::abs(gamma) * field.peak_power();
        double step = std::min(max_step, max_ps / phase_shift);

        while (distance + step < length) {
            linearStep(field, 0.5 * step);
            nonlinearStep(field, step);
            linearStep(field, 0.5 * step);

            distance += step;
            phase_shift = std::abs(gamma) * field.peak_power();
            step = std::min(max_step, max_ps / phase_shift);
        }

        step = length - distance;
        linearStep(field, 0.5 * step);
        nonlinearStep(field, step);
        linearStep(field, 0.5 * step);
    }
}

void Fiber::compensateCD(Field& field, const double& times) const {
    int samples = field.size();
    double freq;
    double freq_step = 2 * M_PI * field.getSamplingRate() / samples;

    field.fft_inplace();
    for (int i = 0; i < samples; ++i) {
        if (i <= samples / 2)
            freq = freq_step * i;
        else
            freq = freq_step * (i - samples);

        field[i] *= i_exp(0.5 * beta2 * freq * freq * times * length);
    }
    field.ifft_inplace();
}

void Fiber::amplify(Field& field) const { field *= exp(0.5 * alpha * length); }
