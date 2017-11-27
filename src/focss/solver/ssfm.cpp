#include "ssfm.h"
#include <cmath>
#include <iostream>
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/module/fiber.h"

namespace focss {
const double SSFM::FORWARD = 1;
const double SSFM::BACKWARD = -1;

SSFM::SSFM() {}

SSFM::SSFM(const Fiber& fiber) : fiber_(fiber) {}

void SSFM::set_fiber(const Fiber& fiber) { fiber_ = fiber; }

void SSFM::set_total_steps(const int& total_steps) {
    uniform_grid_ = true;
    total_steps_ = total_steps;
}

void SSFM::set_maximum_shift_and_step(const double& maximum_phase_shift,
                                      const double& maximum_step) {
    uniform_grid_ = false;
    maximum_phase_shift_ = maximum_phase_shift;
    maximum_step_ = maximum_step;
}

void SSFM::run(Field& field, Direction direction) const {
    if (uniform_grid_) {
        uniform_solve(field, direction);
    } else {
        adaptive_solve(field, direction);
    }
}

void SSFM::uniform_solve(Field& field, Direction direction) const {
    double step = fiber_.length / total_steps_;
    step = direction * step;

    linear_step(field, 0.5 * step);
    for (int i = 0; i < total_steps_; ++i) {
        nonlinear_step(field, step);
        linear_step(field, step);
    }
    linear_step(field, -0.5 * step);
}

void SSFM::adaptive_solve(Field& field, Direction direction) const {
    double distance = 0;
    double next_step = estimate_next_step(field);
    double step = next_step;
    int steps_done = 1;

    field.fft_inplace();
    while (distance + step < fiber_.length) {
        linear_step_without_fft(field, 0.5 * step * direction);
        field.ifft_inplace();
        nonlinear_step(field, step * direction);
        next_step = estimate_next_step(field);
        field.fft_inplace();
        linear_step_without_fft(field, 0.5 * step * direction);

        distance += step;
        step = next_step;
        steps_done++;
    }

    step = fiber_.length - distance;
    linear_step_without_fft(field, 0.5 * step * direction);
    nonlinear_step_with_fft(field, step * direction);
    linear_step_without_fft(field, 0.5 * step * direction);
    field.ifft_inplace();
    std::cout << "    * SSFM done " << steps_done << " steps" << std::endl;
}

double SSFM::estimate_next_step(const Field& field) const {
    double argument = fiber_.kappa * std::abs(fiber_.gamma);
    double phase_shift = argument * field.peak_power();
    return std::min(maximum_step_, maximum_phase_shift_ / phase_shift);
}

void SSFM::linear_step(Field& field, const double& step) const {
    field.fft_inplace();
    linear_step_without_fft(field, step);
    field.ifft_inplace();
}

void SSFM::nonlinear_step(Field& field, const double& step) const {
    double argument = fiber_.kappa * fiber_.gamma * step;

    for (int i = 0; i < field.samples(); ++i) {
        field(i) *= i_exp(argument * field.power(i));
    }
}

void SSFM::linear_step_without_fft(Field& field, const double& step) const {
    field *= exp(-0.5 * fiber_.alpha * step);

    double argument = 0.5 * fiber_.beta2 * step;
    for (int i = 0; i < field.samples(); ++i)
        field(i) *= i_exp(argument * field.w(i) * field.w(i));
}

void SSFM::nonlinear_step_with_fft(Field& field, const double& step) const {
    field.ifft_inplace();
    nonlinear_step(field, step);
    field.fft_inplace();
}
}  // namespace focss
