#include "ssfm.h"
#include <cmath>
#include <iostream>
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/module/fiber.h"

namespace focss {
const int SSFM::UNIFORM = 0;
const int SSFM::LOGARITHMIC = 1;
const int SSFM::ADAPTIVE = 2;

const double SSFM::FORWARD = 1.0;
const double SSFM::BACKWARD = -1.0;

SSFM::SSFM() {}

SSFM::SSFM(const Fiber& fiber) : fiber_(fiber) {}

SSFM::SSFM(const Fiber& fiber, const int& total_steps)
    : fiber_(fiber), grid_(UNIFORM), total_steps_(total_steps) {}

void SSFM::set_fiber(const Fiber& fiber) { fiber_ = fiber; }

void SSFM::set_grid(const Grid& grid) { grid_ = grid; }

void SSFM::set_total_steps(const int& total_steps) {
    total_steps_ = total_steps;
}

void SSFM::set_maximum_phase_shift(const double& maximum_phase_shift) {
    maximum_phase_shift_ = maximum_phase_shift;
}

void SSFM::set_maximum_step_size(const double& maximum_step_size) {
    maximum_step_size_ = maximum_step_size;
}

void SSFM::run(Field& field, const Direction& direction) {
    int freq_step = field.dw();
    int cache_size = field.samples();

    disp_cache_ = new double[cache_size];
    for (int i = 0; i <= cache_size / 2; ++i)
        disp_cache_[i] = static_cast<double>(i) * freq_step;
    for (int i = cache_size / 2 + 1; i < cache_size; ++i)
        disp_cache_[i] = static_cast<double>(i - cache_size) * freq_step;
    for (int i = 0; i < cache_size; ++i)
        disp_cache_[i] *= disp_cache_[i];

    if (grid_ == UNIFORM) {
        uniform_solve(field, direction);
    } else if (grid_ == LOGARITHMIC) {
        logarithmic_solve(field, direction);
    } else {  // grid_ == ADAPTIVE
        adaptive_solve(field, direction);
    }

    delete[] disp_cache_;
}

void SSFM::uniform_solve(Field& field, const Direction& direction) const {
    double step = fiber_.length / total_steps_;
    step = direction * step;

    linear_step(field, 0.5 * step);
    for (int i = 0; i < total_steps_; ++i) {
        nonlinear_step(field, step);
        linear_step(field, step);
    }
    linear_step(field, -0.5 * step);
}

void SSFM::logarithmic_solve(Field& field, const Direction& direction) const {
    auto logarithmic_steps = estimate_logarithmic_steps(field);
    auto inv_alpha = -1.0 / fiber_.alpha * direction;
    auto inv_gain = std::exp(-fiber_.alpha * fiber_.length) - 1;
    inv_gain /= logarithmic_steps;

    double step;
    field.fft_inplace();
    for (int i = 0; i < logarithmic_steps; ++i) {
        step = inv_alpha * std::log(1 + inv_gain / (1 + inv_gain * i));

        nofft_linear_step(field, 0.5 * step);
        fft_nonlinear_step(field, step);
        nofft_linear_step(field, 0.5 * step);
    }
    field.ifft_inplace();
}

void SSFM::adaptive_solve(Field& field, const Direction& direction) const {
    double distance = 0;
    double next_step = estimate_next_adaptive_step(field);
    double step = next_step;
    int steps_done = 1;

    field.fft_inplace();
    while (distance + step < fiber_.length) {
        nofft_linear_step(field, 0.5 * step * direction);
        field.ifft_inplace();
        nonlinear_step(field, step * direction);
        next_step = estimate_next_adaptive_step(field);
        field.fft_inplace();
        nofft_linear_step(field, 0.5 * step * direction);

        distance += step;
        step = next_step;
        steps_done++;
    }

    step = fiber_.length - distance;
    nofft_linear_step(field, 0.5 * step * direction);
    fft_nonlinear_step(field, step * direction);
    nofft_linear_step(field, 0.5 * step * direction);
    field.ifft_inplace();
    std::cout << "    * SSFM done " << steps_done << " steps" << std::endl;
}

int SSFM::estimate_logarithmic_steps(const Field& field) const {
    auto alpha = fiber_.alpha;
    auto L_nl = 1.0 / fiber_.gamma / field.average_power();
    auto CL_nl_kappa = maximum_phase_shift_ * L_nl / fiber_.kappa;

    auto L_eff = (1.0 - std::exp(-alpha * fiber_.length)) / alpha;
    auto L_nl_eff = (1.0 - std::exp(-alpha * CL_nl_kappa)) / alpha;

    auto estimated_steps = static_cast<int>(L_eff / L_nl_eff);
    auto minimum_steps = static_cast<int>(fiber_.length / maximum_step_size_);

    std::cout << "    * predicting  " << estimated_steps;
    std::cout << " SSFM steps" << std::endl;
    return std::max(minimum_steps, estimated_steps);
}

double SSFM::estimate_next_adaptive_step(const Field& field) const {
    double argument = fiber_.kappa * std::abs(fiber_.gamma);
    double phase_shift = argument * field.peak_power();
    return std::min(maximum_step_size_, maximum_phase_shift_ / phase_shift);
}

void SSFM::linear_step(Field& field, const double& step) const {
    field *= exp(-0.5 * fiber_.alpha * step);

    field.fft_inplace();
    double argument = 0.5 * fiber_.beta2 * step;
    for (int index = 0, mode = 0; mode < field.modes(); ++mode)
        for (int sample = 0; sample < field.samples(); ++sample, ++index)
            field[index] *= i_exp(argument * disp_cache_[sample]);
    field.ifft_inplace();
}

void SSFM::nonlinear_step(Field& field, const double& step) const {
    double argument = fiber_.kappa * fiber_.gamma * step;
    for (int index = 0, mode = 0; mode < field.modes(); ++mode)
        for (int sample = 0; sample < field.samples(); ++sample, ++index)
            field[index] *= i_exp(argument * field.power(sample));
}

void SSFM::nofft_linear_step(Field& field, const double& step) const {
    field *= exp(-0.5 * fiber_.alpha * step);

    double argument = 0.5 * fiber_.beta2 * step;
    for (int index = 0, mode = 0; mode < field.modes(); ++mode)
        for (int sample = 0; sample < field.samples(); ++sample, ++index)
            field[index] *= i_exp(argument * disp_cache_[sample]);
}

void SSFM::fft_nonlinear_step(Field& field, const double& step) const {
    field.ifft_inplace();
    double argument = fiber_.kappa * fiber_.gamma * step;
    for (int index = 0, mode = 0; mode < field.modes(); ++mode)
        for (int sample = 0; sample < field.samples(); ++sample, ++index)
            field[index] *= i_exp(argument * field.power(sample));
    field.fft_inplace();
}
}  // namespace focss
