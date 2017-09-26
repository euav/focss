#include "modulation.h"

Field random_16qam_symbols(const unsigned long& length,
                           const double& baudrate) {
    Field symbols(length);
    for (unsigned long i = 0; i < length; ++i)
        symbols[i] = gray_symbols_16qam[i % 16];

    std::shuffle(symbols.begin(), symbols.end(), std::mt19937(time(0)));
    symbols.setSampligRate(baudrate);
    return symbols;
}

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

Field rrc_filter(const double& roll_off,
                 const unsigned long& width,
                 const unsigned long& osf) {
    Field filter(width * osf, 0);

    for (unsigned long i = 0; i < width * osf; ++i) {
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

Field rrc_modulate(const Field& symbols,
                   const unsigned long& osf,
                   const double& power) {
    Field rrc = rrc_filter(rrc_roll_off, symbols.size(), osf);
    Field signal = symbols.upsample(osf);
    signal.apply_filter(rrc);
    signal *= std::sqrt(power / signal.average_power());

    return signal;
}

Field rrc_demodulate(const Field& signal, const unsigned long& osf) {
    Field rrc = rrc_filter(rrc_roll_off, signal.size() / osf, osf);
    Field symbols = signal;
    symbols = symbols.apply_filter(rrc).downsample(osf);
    symbols *= std::sqrt(1.0 / symbols.average_power());

    return symbols;
}

void lowpass_inplace(Field& field, const double& cutoff_frequency) {
    double time_domain = field.size() * field.dt();
    unsigned long cutoff_index = cutoff_frequency * time_domain;

    field.fft_inplace();
    for (unsigned long i = cutoff_index + 1; i < field.size() - cutoff_index;
         ++i)
        field[i] = 0;
    field.ifft_inplace();
}
