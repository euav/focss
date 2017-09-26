#include "ssfm.h"

SSFM::SSFM() {}

SSFM::SSFM(const Fiber& fiber) : fiber(fiber) {}

void SSFM::setFiber(const Fiber& fiber) {
    this->fiber = fiber;
}

void SSFM::setTotalSteps(const int& steps) {
    uniform_grid = true;
    total_steps = steps;
}

void SSFM::setMaximumShiftAndStep(const double& shift, const double& step) {
    uniform_grid = false;
    max_ps = shift;
    max_step = step;
}

void SSFM::linearStep(Field& field, const double& step) const {
    field *= exp(-0.5 * fiber.alpha * step);

    field.fft_inplace();
    for (unsigned long i = 0; i < field.size(); ++i)
        field[i] *= i_exp(-0.5 * fiber.beta2 * field.w(i) * field.w(i) * step);
    field.ifft_inplace();
}

void SSFM::nonlinearStep(Field& field, const double& step) const {
    for (unsigned long j = 0; j < field.size(); ++j)
        field[j] *= i_exp(fiber.gamma * step * norm(field[j]));
}

void SSFM::run(Field& field) const {
    if (uniform_grid) {
        double step = fiber.length / total_steps;

        linearStep(field, 0.5 * step);
        for (int i = 0; i < total_steps; ++i) {
            nonlinearStep(field, step);
            linearStep(field, step);
        }
        linearStep(field, -0.5 * step);
    } else {
        double distance = 0;
        double phase_shift = std::abs(fiber.gamma) * field.peak_power();
        double step = std::min(max_step, max_ps / phase_shift);

        while (distance + step < fiber.length) {
            linearStep(field, 0.5 * step);
            nonlinearStep(field, step);
            linearStep(field, 0.5 * step);

            distance += step;
            phase_shift = std::abs(fiber.gamma) * field.peak_power();
            step = std::min(max_step, max_ps / phase_shift);
        }

        step = fiber.length - distance;
        linearStep(field, 0.5 * step);
        nonlinearStep(field, step);
        linearStep(field, 0.5 * step);
    }
}
