#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include "headers/fiber.h"
#include "headers/linear_equalizer.h"
#include "headers/lms_equalizer.h"
#include "headers/modulation.h"
#include "headers/Field.h"
#include "headers/utility.h"

// reference units [m], [s], [W]
const double bandwidth = 32e9;       // [baud]
const double distance = 1e5;         // [m]
const double spans = 10;             // [1]
const double noise_figure = 4.5;     // [dB]
const double attenuation = 2e-4;     // [dB/m]
const double dispersion = 1.7e-5;    // [s/m/m]
const double wavelength = 1.55e-6;   // [m]
const double nonlinearity = (16.0 / 9.0) * 1.4e-3;  // [1/W/m]
const double rrc_roll_off = 0.01;    // [1]
const int pulse_width = 1024;        // [1]
const int oversampling = 16;         // [1]
const double forward_steps = 1000;   // [1]
const double backward_steps = 100;    // [1]

Field random_16qam_symbols(const int& constellations) {
    Field symbols(16 * constellations);
    for (int i = 0; i < constellations; ++i) {
        for (int j = 0; j < 16; ++j) {
            symbols[16 * i + j] = gray_symbols_16qam[j];
        }
    }
    std::srand(time(0));
    std::random_shuffle(symbols.begin(), symbols.end());

    return symbols;
}

Field modulate_rrc(const Field& symbols, const double& launch_power) {
    Field rrc = rrc_filter(oversampling, rrc_roll_off, pulse_width);
    Field field = symbols.upsample(oversampling);
    field.apply_filter(rrc);
    field *= std::sqrt(launch_power / field.average_power());

    return field;
}

Field match_rrc(const Field& field, const int& osf) {
    Field rrc = rrc_filter(osf, rrc_roll_off, pulse_width);
    Field matched = field;
    matched = matched.apply_filter(rrc).downsample(osf);
    matched *= std::sqrt(1.0 / matched.average_power());

    return matched;
}

void forward_propagation(Field& field) {
    std::cout << ">> forward propagation" << std::endl;

    double noise_variance = 0.5 * pow(10, noise_figure / 10);
    noise_variance *= pow(10, attenuation * distance / 10) - 1;
    noise_variance *= planck * light_speed / wavelength;
    noise_variance *= oversampling * bandwidth;

    Fiber fiber(wavelength);
    fiber.setAttenuationDB(attenuation);
    fiber.setDispersionPhysical(dispersion);
    fiber.setNonlinearity(nonlinearity);
    fiber.setFiberLength(distance);
    fiber.setTotalSteps(forward_steps);
    fiber.setSamplingRate(oversampling * bandwidth);
    for (int i = 0; i < spans; ++i) {
        fiber.propagate(field);
        fiber.amplify(field);
        awgn_generator(field, noise_variance);
        std::cout << "    * span " << i + 1 << " has done" << std::endl;
    }
    std::cout << std::endl;
}

void cd_compensation(Field& field) {
    Fiber fiber(wavelength);
    fiber.setAttenuationDB(attenuation);
    fiber.setDispersionPhysical(dispersion);
    fiber.setNonlinearity(nonlinearity);
    fiber.setFiberLength(distance);
    fiber.setTotalSteps(forward_steps);
    fiber.setSamplingRate(oversampling * bandwidth);
    fiber.compensateCD(field, spans);
}

void backward_propagation(Field& field) {
    std::cout << ">> backward propagation" << std::endl;

    Fiber fiber(wavelength);
    fiber.setAttenuationDB(-attenuation);
    fiber.setDispersionPhysical(-dispersion);
    fiber.setNonlinearity(-nonlinearity);
    fiber.setFiberLength(distance);
    fiber.setTotalSteps(backward_steps);
    fiber.setSamplingRate(2 * bandwidth);
    for (int i = 0; i < spans; ++i) {
        fiber.amplify(field);
        fiber.propagate(field);
        std::cout << "    * span " << spans - i << " has done" << std::endl;
    }
    std::cout << std::endl;
}

void tx_rx(const int& constellations, const double& launch_power) {
    std::cout.precision(15);
    Field symbols = random_16qam_symbols(constellations);
    Field field = modulate_rrc(symbols, 0.5 * launch_power);

    // propagation
    forward_propagation(field);
    Field dbp_field = field.downsample(8);
    backward_propagation(dbp_field);
    cd_compensation(field);

    // Field matching
    field = match_rrc(field, oversampling);
    dbp_field = match_rrc(dbp_field, 2);
    std::cout << ">> Fields has matched" << std::endl;

    // linear equalization
    LinearEqualizer leq;
    leq.train(symbols, field);
    Field leq_field = leq.equalize(field);
    std::cout << ">> linear equalization has done" << std::endl;

    // saving results
    std::string filename = "results/output.txt";
    std::ofstream fout(filename);
    fout.precision(15);

    fout << "# [ AVG POWER, CD, LEQ, DBP ] = " << launch_power
         << " ";
    fout << q2_factor(symbols, field) << " ";
    fout << q2_factor(symbols, leq_field) << " ";
    fout << q2_factor(symbols, dbp_field) << std::endl;

    for (int i = 0; i < symbols.size(); ++i) {
        fout << symbols[i].real() << " ";
        fout << symbols[i].imag() << " ";
        fout << field[i].real() << " ";
        fout << field[i].imag() << " ";
        fout << dbp_field[i].real() << " ";
        fout << dbp_field[i].imag() << " ";
        fout << leq_field[i].real() << " ";
        fout << leq_field[i].imag() << "\n";
    }
    fout.close();
    std::cout << ">> results has saved" << std::endl;
}

void check_pulse_and_match_sampling() {
    // generate random Field
    int constellations = 256;
    Field symbols(16 * constellations);
    for (int i = 0; i < constellations; ++i) {
        for (int j = 0; j < 16; ++j) {
            symbols[16 * i + j] = gray_symbols_16qam[j];
        }
    }
    std::srand(time(0));
    std::random_shuffle(symbols.begin(), symbols.end());

    Field field = symbols.upsample(oversampling);
    Field rrc16 = rrc_filter(oversampling, rrc_roll_off, pulse_width);
    field.apply_filter(rrc16);

    // match Field and output constellations
    field = field.downsample(8);
    Field rrc2 = rrc_filter(2, rrc_roll_off, pulse_width);
    field.apply_filter(rrc2);
    field = field.downsample(2);
    field *= std::sqrt(1 / field.average_power());
    for (int i = 0; i < field.size(); ++i) {
        std::cout << field[i].real() << " ";
        std::cout << field[i].imag() << "\n";
    }
}

void check_different_sps_filter_spectra(const int& sps) {
    Field rrc = rrc_filter(sps, rrc_roll_off, pulse_width);
    rrc.fft_shift().fft_inplace().fft_shift();

    double T = pulse_width / bandwidth;
    for (int i = 0; i < rrc.size(); ++i) {
        std::cout << i / T - 0.5 * rrc.size() / T << " ";
        std::cout << rrc[i].real() << "\n";
    }
}

void check_different_sps_spectra_broadering(const int& sps) {
    int constellations = 64;
    Field symbols(16 * constellations);
    for (int i = 0; i < constellations; ++i) {
        for (int j = 0; j < 16; ++j) {
            symbols[16 * i + j] = gray_symbols_16qam[j];
        }
    }
    std::srand(time(0));
    std::random_shuffle(symbols.begin(), symbols.end());

    Field field = symbols.upsample(sps);
    Field rrc = rrc_filter(sps, rrc_roll_off, pulse_width);
    field.apply_filter(rrc);
    field *= std::sqrt(dbm_to_watts(0) / field.average_power());

    Fiber fiber(wavelength);
    fiber.setAttenuationDB(attenuation);
    fiber.setDispersionPhysical(dispersion);
    fiber.setNonlinearity(nonlinearity);
    fiber.setFiberLength(distance);
    fiber.setTotalSteps(forward_steps);
    fiber.setSamplingRate(sps * bandwidth);
    for (int i = 0; i < spans; ++i) {
        fiber.propagate(field);
        fiber.amplify(field);
    }

    field.fft_inplace().fft_shift();
    double T = pulse_width / bandwidth;
    for (int i = 0; i < rrc.size(); ++i) {
        std::cout << i / T - 0.5 * field.size() / T << " ";
        std::cout << std::sqrt(norm(field[i])) << "\n";
    }
}

void check_noise_power(const int& N) {
    Field zero(N, 0);
    awgn_generator(zero, 4 * N);
    std::cout << zero.average_power() << std::endl;
    zero.fft_inplace();
    std::cout << zero.average_power() << std::endl;
}

void check_noise_variance() {
    int constellations = 256;
    Field symbols(16 * constellations);
    for (int i = 0; i < constellations; ++i) {
        for (int j = 0; j < 16; ++j) {
            symbols[16 * i + j] = gray_symbols_16qam[j];
        }
    }
    std::srand(time(0));
    std::random_shuffle(symbols.begin(), symbols.end());

    Field field = symbols.upsample(oversampling);
    Field rrc = rrc_filter(oversampling, rrc_roll_off, pulse_width);
    field.apply_filter(rrc);
    field *= std::sqrt(dbm_to_watts(0) / field.average_power());

    // adding gaussian noise
    double noise_variance = oversampling * pow(10, noise_figure / 10);
    noise_variance *= pow(10, attenuation * std::log(10) * distance / 100) - 1;
    noise_variance *= bandwidth * planck * light_speed / wavelength;

    // field.fft_inplace();
    awgn_generator(field, 10 * noise_variance);
    // field.ifft_inplace();

    field = field.apply_filter(rrc).downsample(16);
    field *= std::sqrt(1 / field.average_power());
    std::cout << "# Q2 = " << q2_factor(symbols, field) << std::endl;
    for (int i = 0; i < field.size(); ++i) {
        std::cout << field[i].real() << " ";
        std::cout << field[i].imag() << "\n";
    }
}

void check_width() {
    int constellations = 64;
    Field symbols(16 * constellations);
    for (int i = 0; i < constellations; ++i) {
        for (int j = 0; j < 16; ++j) {
            symbols[16 * i + j] = gray_symbols_16qam[j];
        }
    }
    std::srand(time(0));
    std::random_shuffle(symbols.begin(), symbols.end());

    Field field = symbols.upsample(oversampling);
    Field rrc = rrc_filter(oversampling, rrc_roll_off, pulse_width);
    field.apply_filter(rrc);
    field *= std::sqrt(dbm_to_watts(0) / field.average_power());
    field.apply_filter(rrc);
    field *= std::sqrt(1 / field.average_power());
    field.fft_inplace().fft_shift();

    double freq_step = oversampling * bandwidth / field.size();
    for (int i = 0; i < field.size(); ++i) {
        std::cout <<  (i * freq_step - 0.5 * field.size() * freq_step) / 1e9 << " ";
        std::cout << norm(field[i]) << "\n";
    }
}

void check_lms() {
    Field symbols = random_16qam_symbols(64);
    Field field = modulate_rrc(symbols, 0.5 * dbm_to_watts(0));
    forward_propagation(field);
    cd_compensation(field);
    field = match_rrc(field, oversampling);
    std::cout << ">> Fields has matched" << std::endl;

    // linear equalization
    LinearEqualizer leq;
    leq.train(symbols, field);
    Field leq_field = leq.equalize(field);
    std::cout << ">> linear equalization has done" << std::endl;

    // lms equalization
    LmsEqualizer lms(3);
    lms.train(symbols, field);
    Field lms_field = lms.equalize(field);
    std::cout << ">> least-mean square equalization has done" << std::endl;

    std::cout << "Q2 of [CD, LEQ, LMS]" << std::endl;
    std::cout << q2_factor(symbols, field) << std::endl;
    std::cout << q2_factor(symbols, leq_field) << std::endl;
    std::cout << q2_factor(symbols, lms_field) << std::endl;
}

int main() {
    tx_rx(64, dbm_to_watts(0));
    //check_lms();
    return 0;
}
