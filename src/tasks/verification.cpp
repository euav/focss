#include "focss/field.h"
#include "focss/functions.h"
#include "focss/solver/ssfm.h"
using namespace focss;

Fiber fiber = {0.0, 1.0, 1.0, 1.0, 0.0, 1.0};
double time_domain = 10.0;
double power = 1.0;

Field sech_pulse(const int& points, const double& z) {
    double dt = time_domain / static_cast<double>(points);
    Domain domain{1, points, dt, 0.0};
    Field field(domain);

    for (int i = 0; i < field.size(); ++i) {
        double t = static_cast<double>(i - field.size() / 2) * dt;
        double a = std::sqrt(power);
        field[i] = a * i_exp(0.5 * power * z) / std::cosh(a * t);
    }

    return field;
}

Field dark_pulse(const int& points, const double& z) {
    double dt = time_domain / static_cast<double>(points);
    Domain domain{1, points, dt, 0.0};
    Field field(domain);

    for (int i = 0; i < field.size(); ++i) {
        double t = static_cast<double>(i - field.size() / 2) * dt;
        double a = std::sqrt(power);
        field[i] = a * i_exp(power * z) * std::tanh(a * t);
    }

    return field;
}

void verification(const int& n_samples, const int& n_steps) {
    Field field = dark_pulse(n_samples, 0);
    Field solution = dark_pulse(n_samples, fiber.length);

    SSFM ssfm;
    ssfm.set_fiber(fiber);
    ssfm.set_grid(SSFM::UNIFORM);
    ssfm.set_total_steps(n_steps);
    ssfm.run(field, SSFM::FORWARD);

    double error = std::sqrt((field - solution).average_power() /
                             solution.average_power());

    std::cout << error << std::endl;
}
