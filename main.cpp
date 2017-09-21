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
const int backward_steps = 10;        // [1]
const int training_length = 100;      // [1]

Field random_16qam_symbols(const int& constellations) {
    Field symbols(16 * constellations);
    for (int i = 0; i < constellations; ++i) {
        for (int j = 0; j < 16; ++j) {
            symbols[16 * i + j] = gray_symbols_16qam[j];
        }
    }
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
    Field rrc = rrc_filter(rrc_roll_off, pulse_width, osf);
    Field matched = signal;
    matched = matched.apply_filter(rrc).downsample(osf);
    matched *= std::sqrt(1.0 / matched.average_power());

    return matched;
}

void bandpass_inplace(Field& field, const double& cutoff_frequency) {
    double time_domain = field.size() / field.getSamplingRate();
    int cutoff_index = 1 + int(cutoff_frequency * time_domain);

    field.fft_inplace();
    for (int i = cutoff_index; i < field.size() - cutoff_index; ++i)
        field[i] = 0;
    field.ifft_inplace();
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
        fiber.setTotalSteps(backward_steps);

    for (int i = 0; i < spans; ++i) {
        fiber.amplify(signal);
        fiber.propagate(signal);
    }
}

void experiment(const int& constellations,
                const double& average_power,
                std::ostream& os) {
    Field symbols = random_16qam_symbols(constellations);
    Field signal = modulate_rrc(symbols, dbm_to_watts(average_power));
    os << average_power << ',';

    // propagation
    forward_propagation(signal);

    std::cout << ">> experimenting with dbp step per span" << std::endl;
    for (int i = 1; i <= 10; ++i) {
        Field dbp_signal = signal;
        bandpass_inplace(dbp_signal, 0.5 * bandwidth);
        // dbp_signal = dbp_signal.downsample(8);
        backward_propagation(dbp_signal, i);
        dbp_signal = demodulate_rrc(dbp_signal, oversampling);

        Field training_x = symbols.chomp(0, symbols.size() - training_length);
        Field training_y = dbp_signal.chomp(0, dbp_signal.size() - training_length);
        LmsEqualizer lms(3);
        lms.train(training_x, training_y);
        Field validation_x = symbols.chomp(training_length, 0);
        Field validation_y = lms.equalize(dbp_signal).chomp(training_length, 0);

        os << q2_factor(validation_x, validation_y) << ',';
    }

    std::cout << ">> experimenting with dbp BPF and PSE drop" << std::endl;
    {
        Field dbp_signal = signal;  //.downsample(8);
        Field bpf_signal = signal;
        bandpass_inplace(bpf_signal, 0.5 * bandwidth);
        // bpf_signal = bpf_signal.downsample(8);
        backward_propagation(dbp_signal, 10);
        backward_propagation(bpf_signal, 10);
        dbp_signal = demodulate_rrc(dbp_signal, oversampling);
        bpf_signal = demodulate_rrc(bpf_signal, oversampling);

        os << q2_factor(symbols, dbp_signal) << ',';
        os << q2_factor(symbols, bpf_signal) << ',';

        Field training_x = symbols.chomp(0, symbols.size() - training_length);
        Field training_y = dbp_signal.chomp(0, dbp_signal.size() - training_length);
        LmsEqualizer lms(3);
        lms.train(training_x, training_y);
        Field validation_x = symbols.chomp(training_length, 0);
        Field validation_y = lms.equalize(dbp_signal).chomp(training_length, 0);


        os << q2_factor(validation_x, validation_y) << ',';
    }

    std::cout << ">> experimenting with lms and pse" << std::endl;
    {
        Field cd_signal = signal;
        cd_compensation(cd_signal);
        cd_signal = demodulate_rrc(cd_signal, oversampling);

        for (int i = 1; i <= 3; ++i) {
            Field training_x = symbols.chomp(0, symbols.size() - training_length);
            Field training_y = cd_signal.chomp(0, cd_signal.size() - training_length);
            LmsEqualizer lms(i);
            lms.train(training_x, training_y);
            Field validation_x = symbols.chomp(training_length, 0);
            Field validation_y = lms.equalize(cd_signal).chomp(training_length, 0);

            os << q2_factor(validation_x, validation_y) << ',';
        }

        Field training_x = symbols.chomp(0, symbols.size() - training_length);
        Field training_y = cd_signal.chomp(0, cd_signal.size() - training_length);
        PhaseShiftEqualizer pse(720);
        pse.train(training_x, training_y);
        Field validation_x = symbols.chomp(training_length, 0);
        Field validation_y = pse.equalize(cd_signal).chomp(training_length, 0);

        os << q2_factor(validation_x, validation_y) << std::endl;
    }
    std::cout << ">> experiment has done" << std::endl;
}

int main() {
    /*
    std::ofstream fout("experiment/right.txt");
    fout << "#AVG(dBm),DBP(1-10),DBP(BPF-PSE),LMS(1-3),PSE\n";
    fout.precision(15);

    for (int i = 21; i <= 41; ++i) {
        double average_power = double(i) / 2.0 - 10.0;
        std::cout << "\nstarting experiment at power ";
        std::cout << average_power << " dBm" << std::endl;

        experiment(256, average_power, fout);
    }
    fout.close();
    */

    experiment(64, -2, std::cout);
    return 0;
}
