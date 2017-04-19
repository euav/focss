#include "../headers/modulation.h"
#include "../headers/utility.h"

Signal sech_pulse(const int& nodes_quantity, const double& width) {
    Signal sech_pulse(nodes_quantity, 0);
    double argument;
    for (int i = 0; i < nodes_quantity; i++) {
        argument = double(i) / (nodes_quantity - 1.0);
        argument = 20 * argument - 20 / 2.0;
        sech_pulse[i] = 1.0 / std::cosh(argument / width);
        sech_pulse[i] /= width;
    }

    return sech_pulse;
}

Signal rrc_filter(const int& samples, const double& roll_off, const int& N) {
    Signal filter(N * samples, 0);

    for (int i = 0; i < N * samples; ++i) {
        double t = double(i) / samples - N / 2.0;

        if (t == 0) {
            filter[i] = 1 - roll_off + 4 * roll_off / M_PI;
        } else if (std::abs(t) == 0.25 / roll_off) {
            filter[i] += (1 + 2 / M_PI) * std::sin(0.25 * M_PI / roll_off);
            filter[i] += (1 - 2 / M_PI) * std::cos(0.25 * M_PI / roll_off);
            filter[i] *= roll_off / sqrt(2);
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

Signal modulate_16qam(const Information& data) {
    Signal symbols(data.size(), 0);
    for (int i = 0; i < data.size(); ++i) {
        InformationType code = data[i];
        symbols[i] = gray_symbols_16qam[code & 0x0F];
    }

    return symbols;
}

Information demodulate_16qam(const Signal& signal) {
    Information data(signal.size());

    InformationType argmin;
    double distance, min_distance;
    for (int i = 0; i < signal.size(); ++i) {
        argmin = 0;
        min_distance = norm(gray_symbols_16qam[0] - signal[i]);

        for (InformationType code = 1; code < 16; ++code) {
            distance = norm(gray_symbols_16qam[code] - signal[i]);
            if (min_distance > distance) {
                min_distance = distance;
                argmin = code;
            }
        }

        data[i] = argmin;
    }

    return data;
}

double bit_error_rate(const Information& data_tx, const Information& data_rx) {
    double error = 0;
    for (int i = 0; i < data_tx.size(); ++i) {
        // Brian Kernighan algorithm
        InformationType difference = data_tx[i] ^ data_rx[i];
        while (difference) {
            ++error;
            difference &= difference - 1;
        }
    }

    return error / data_tx.size();
}
