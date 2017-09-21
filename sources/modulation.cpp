#include "../headers/modulation.h"
#include "../headers/utility.h"

Field sech_pulse(const int& nodes_quantity, const double& width) {
    Field sech_pulse(nodes_quantity, 0);
    double argument;
    for (int i = 0; i < nodes_quantity; i++) {
        argument = double(i) / (nodes_quantity - 1.0);
        argument = 20 * argument - 20 / 2.0;
        sech_pulse[i] = 1.0 / std::cosh(argument / width);
        sech_pulse[i] /= width;
    }

    return sech_pulse;
}

Field rrc_filter(const double& roll_off, const int& width, const int& osf) {
    Field filter(width * osf, 0);

    for (int i = 0; i < width * osf; ++i) {
        double t = double(i) / osf - width / 2.0;

        if (t == 0) {
            filter[i] = 1 - roll_off + 4 * roll_off / M_PI;
        } else if (std::abs(t) == 0.25 / roll_off) {
            filter[i] += (1 + 2 / M_PI) * std::sin(0.25 * M_PI / roll_off);
            filter[i] += (1 - 2 / M_PI) * std::cos(0.25 * M_PI / roll_off);
            filter[i] *= roll_off / std::sqrt(2.0);
        } else {
            filter[i] += std::cos(M_PI * t * (1 + roll_off));
            filter[i] *= 4 * roll_off * t;
            filter[i] += std::sin(M_PI * t * (1 - roll_off));
            filter[i] /= 1 - 16 * roll_off * roll_off * t * t;
            filter[i] /= M_PI * t;
        }
    }

    return filter *= 1.0 / filter.peak_power();
}
