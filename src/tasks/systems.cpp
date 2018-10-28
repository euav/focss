#include "systems.h"
#include <iostream>
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/module/amplifier.h"
#include "focss/module/fiber.h"
#include "focss/solver/ssfm.h"
#include "parameters.h"

const focss::Fiber fiber = {focss::db_to_natural(::attenuation),
                            focss::disp_to_beta2(::dispersion, ::wavelength),
                            ::nonlinearity,
                            ::fiber_length,
                            ::wavelength,
                            ::polarization_coupling};

void back_to_back(focss::Field& signal) {
    focss::Amplifier amplifier(
        focss::db_to_linear(::attenuation * ::fiber_length),
        focss::db_to_linear(::noise_figure));

    for (int i = 0; i < ::spans; ++i)
        amplifier.add_noise(signal);
}

void forward_propagation(focss::Field& signal) {
    std::cout << ">> forward propagation" << std::endl;

    focss::SSFM ssfm(::fiber);
    ssfm.set_grid(focss::SSFM::ADAPTIVE);
    ssfm.set_maximum_phase_shift(::max_phase_shift);
    ssfm.set_maximum_step_size(::max_forward_step);
    focss::Amplifier amplifier(
        focss::db_to_linear(::attenuation * ::fiber_length),
        focss::db_to_linear(::noise_figure));

    for (int i = 0; i < ::spans; ++i) {
        ssfm.run(signal, focss::SSFM::FORWARD);
        amplifier.amplify(signal);
        std::cout << "    * span " << i + 1 << " has done" << std::endl;
    }
    std::cout << std::endl;
}

void forward_propagation(focss::Field& signal, const int& span) {
    std::cout << ">> forward propagation" << std::endl;

    focss::SSFM ssfm(::fiber);
    ssfm.set_grid(focss::SSFM::ADAPTIVE);
    ssfm.set_maximum_phase_shift(::max_phase_shift);
    ssfm.set_maximum_step_size(::max_forward_step);
    focss::Amplifier amplifier(
            focss::db_to_linear(::attenuation * ::fiber_length),
            focss::db_to_linear(::noise_figure));

    for (int i = 0; i < span; ++i) {
        ssfm.run(signal, focss::SSFM::FORWARD);
        amplifier.amplify(signal);
    }
    std::cout << std::endl;
}

void cd_compensation(focss::Field& signal, const int& n_spans) {
    focss::SSFM ssfm(::fiber.dispersive() * n_spans, 1);
    ssfm.run(signal, focss::SSFM::BACKWARD);
}

void cd_compensation(focss::Field& signal) {
    focss::SSFM ssfm(::fiber.dispersive() * ::spans, 1);
    ssfm.run(signal, focss::SSFM::BACKWARD);
}

void backward_propagation(focss::Field& signal, const int& steps) {
    focss::SSFM ssfm(::fiber);
    focss::Amplifier amplifier(
        focss::db_to_linear(::attenuation * ::fiber_length), 0);

    if (steps == 0) {
        ssfm.set_grid(focss::SSFM::LOGARITHMIC);
        ssfm.set_maximum_phase_shift(::max_phase_shift);
        ssfm.set_maximum_step_size(::max_forward_step);
    } else {
        ssfm.set_grid(focss::SSFM::UNIFORM);
        ssfm.set_total_steps(steps);
    }

    for (int i = 0; i < ::spans; ++i) {
        amplifier.drop_power(signal);
        ssfm.run(signal, focss::SSFM::BACKWARD);
    }
}