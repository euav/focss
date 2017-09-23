#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include "headers/fiber.h"
#include "headers/field.h"
#include "headers/lms_equalizer.h"
#include "headers/modulation.h"
#include "headers/phase_shift_equalizer.h"
#include "headers/utility.h"

// reference units [m], [s], [W]
const double bandwidth = 32e9;        // [baud]
const double fiber_length = 1e5;      // [m]
const double spans = 10;              // [1]
const double noise_figure = 4.5;      // [dB]
const double attenuation = 2e-4;      // [dB/m]
const double dispersion = 1.7e-5;     // [s/m/m]
const double wavelength = 1.55e-6;    // [m]
const double nonlinearity = 1.4e-3;   // [1/W/m]
const double rrc_roll_off = 0.01;     // [1]
const int pulse_width = 1024;         // [1]
const int oversampling = 16;          // [1]
const double max_phase_shift = 1e-3;  // [1]
const double max_forward_step = 1e3;  // [m]
const int training_length = 4096;     // [1]
const int IDEAL = 0;

void lowpass_inplace(Field& field, const double& cutoff_frequency) {
    double time_domain = field.size() / field.getSamplingRate();
    int cutoff_index = int(cutoff_frequency * time_domain);
    std::cout << "+ cutoff: " << cutoff_index << "\n";

    field.fft_inplace();
    for (int i = cutoff_index + 1; i < field.size() - cutoff_index; ++i)
        field[i] = 0;
    field.ifft_inplace();
}

Field random_16qam_symbols(const int& length) {
    Field symbols(length);
    for (int i = 0; i < length; ++i)
        symbols[i] = gray_symbols_16qam[i % 16];

    std::shuffle(symbols.begin(), symbols.end(), std::mt19937(time(0)));
    return symbols;
}

Field modulate_rrc(const Field& symbols, const double& launch_power) {
    Field rrc = rrc_filter(rrc_roll_off, pulse_width, oversampling);
    Field signal = symbols.upsample(oversampling);
    signal.apply_filter(rrc);
    signal *= std::sqrt(launch_power / signal.average_power());
    signal.setSampligRate(bandwidth * oversampling);

    return signal;
}

Field demodulate_rrc(const Field& signal, const int& osf) {
    Field symbols = signal;
    Field rrc = rrc_filter(rrc_roll_off, pulse_width, osf);
    symbols = symbols.apply_filter(rrc).downsample(osf);
    symbols *= std::sqrt(1.0 / symbols.average_power());

    return symbols;
}

void forward_propagation(Field& signal) {
    std::cout << ">> forward propagation" << std::endl;

    double noise_variance = 0.5 * pow(10, noise_figure / 10);
    noise_variance *= pow(10, attenuation * fiber_length / 10) - 1;
    noise_variance *= planck * light_speed / wavelength;
    noise_variance *= signal.getSamplingRate();

    Fiber fiber(wavelength);
    fiber.setAttenuationDB(attenuation);
    fiber.setDispersionPhysical(dispersion);
    fiber.setNonlinearity(nonlinearity);
    fiber.setFiberLength(fiber_length);
    fiber.setMaximumShiftAndStep(max_phase_shift, max_forward_step);
    for (int i = 0; i < spans; ++i) {
        fiber.propagate(signal);
        fiber.amplify(signal);
        awgn_generator(signal, noise_variance);
        std::cout << "    * span " << i + 1 << " has done" << std::endl;
    }
    std::cout << std::endl;
}

void cd_compensation(Field& signal) {
    Fiber fiber(wavelength);
    fiber.setAttenuationDB(attenuation);
    fiber.setDispersionPhysical(dispersion);
    fiber.setNonlinearity(nonlinearity);
    fiber.setFiberLength(fiber_length);
    fiber.compensateCD(signal, spans);
}

void backward_propagation(Field& signal, const int& steps = 0) {
    Fiber fiber(wavelength);
    fiber.setAttenuationDB(-attenuation);
    fiber.setDispersionPhysical(-dispersion);
    fiber.setNonlinearity(-nonlinearity);
    fiber.setFiberLength(fiber_length);
    if (steps > 0)
        fiber.setTotalSteps(steps);
    else
        fiber.setMaximumShiftAndStep(max_phase_shift, max_forward_step);

    for (int i = 0; i < spans; ++i) {
        fiber.amplify(signal);
        fiber.propagate(signal);
    }
}

void experiment(const int& symbols_length,
                const double& average_power_dbm,
                std::ostream& output_stream) {
    output_stream << average_power_dbm << ',';

    // propagation
    Field symbols = random_16qam_symbols(symbols_length);
    Field signal = modulate_rrc(symbols, dbm_to_watts(average_power_dbm));
    forward_propagation(signal);

    std::cout << ">> experimenting with dbp" << std::endl;
    for (int i = 1; i <= 11; ++i) {
        Field dbp_signal = signal;
        lowpass_inplace(dbp_signal, 0.5 * bandwidth);
        dbp_signal = dbp_signal.downsample(8);

        // sps parameter
        if (i <= 10)
            backward_propagation(dbp_signal, i);
        else
            backward_propagation(dbp_signal, IDEAL);

        dbp_signal = demodulate_rrc(dbp_signal, 2);

        Field train_x = symbols.chomp(0, symbols_length - training_length);
        Field train_y = dbp_signal.chomp(0, symbols_length - training_length);
        PhaseShiftEqualizer pse(720);
        pse.train(train_x, train_y);
        Field valid_x = symbols.chomp(training_length, 0);
        Field valid_y = pse.equalize(dbp_signal).chomp(training_length, 0);

        output_stream << q2_factor(valid_x, valid_y) << ',';
    }

    std::cout << ">> experimenting with lms and pse" << std::endl;
    {
        Field cd_signal = signal;
        cd_compensation(cd_signal);
        cd_signal = demodulate_rrc(cd_signal, oversampling);

        Field train_x = symbols.chomp(0, symbols_length - training_length);
        Field train_y = cd_signal.chomp(0, symbols_length - training_length);
        PhaseShiftEqualizer pse(720);
        LmsEqualizer lms(3);

        pse.train(train_x, train_y);
        lms.train(train_x, train_y);

        Field valid_x = symbols.chomp(training_length, 0);
        Field valid_pse = pse.equalize(cd_signal).chomp(training_length, 0);
        Field valid_lms = lms.equalize(cd_signal).chomp(training_length, 0);

        output_stream << q2_factor(valid_x, valid_pse) << ',';
        output_stream << q2_factor(valid_x, valid_lms) << std::endl;
    }
    std::cout << ">> experiment has done" << std::endl;
}

int main() {
    std::ofstream fout("experiment/right.txt");
    fout << "#AVG(dBm),DBP(1-10,ideal),PSE,LMS(3)\n";
    fout.precision(15);

    for (int i = 31; i <= 40; ++i) {
        double average_power_dbm = double(i) / 2.0 - 10.0;
        std::cout << "\nstarting experiment at power ";
        std::cout << average_power_dbm << " dBm" << std::endl;

        experiment(16384, average_power_dbm, fout);
    }
    fout.close();

    return 0;
}
