#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include "headers/fiber.h"
#include "headers/linear_equalizer.h"
#include "headers/modulation.h"
#include "headers/signal.h"
#include "headers/utility.h"

// reference units [m], [s], [W]
const double bandwidth = 32e9;       // [baud]
const double distance = 1e5;         // [m]
const double spans = 10;             // [1]
const double noise_figure = 4.5;     // [dB]
const double attenuation = 2e-4;     // [dB/m]
const double dispersion = 1.7e-5;    // [s/m/m]
const double wavelength = 1.55e-6;   // [m]
const double nonlinearity = 1.4e-3;  // [1/W/m]
const double rrc_roll_off = 0.01;    // [1]
const int pulse_width = 1024;        // [1]
const int oversampling = 16;         // [1]
const int forward_steps = 1000;      // [1]
const int backward_steps = 1000;     // [1]

void forward_propagation(Signal& field) {
    std::cout << ">> forward propagation" << std::endl;

    double noise_variance = 2 * pow(10, noise_figure / 10);
    noise_variance *= pow(10, attenuation * distance / 10) - 1;
    noise_variance *= bandwidth * planck * light_speed / wavelength;

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

void cd_compensation(Signal& field) {
    Fiber fiber(wavelength);
    fiber.setAttenuationDB(attenuation);
    fiber.setDispersionPhysical(dispersion);
    fiber.setNonlinearity(nonlinearity);
    fiber.setFiberLength(distance);
    fiber.setTotalSteps(forward_steps);
    fiber.setSamplingRate(oversampling * bandwidth);
    fiber.compensateCD(field, spans);
}

void backward_propagation(Signal& field) {
    std::cout << ">> backward propagation" << std::endl;

    Fiber fiber(wavelength);
    fiber.setAttenuationDB(-attenuation);
    fiber.setDispersionPhysical(-dispersion);
    fiber.setNonlinearity(-nonlinearity);
    fiber.setFiberLength(distance);
    fiber.setTotalSteps(backward_steps);
    fiber.setSamplingRate(oversampling * bandwidth);
    for (int i = 0; i < spans; ++i) {
        fiber.amplify(field);
        fiber.propagate(field);
        std::cout << "    * span " << spans - i << " has done" << std::endl;
    }
    std::cout << std::endl;
}

void tx_rx(const int& constellations, const int& point) {
    std::cout.precision(15);
    std::cout << "Trasmission number " << point << std::endl;

    // symbols generation
    Signal symbols(16 * constellations);
    for (int i = 0; i < constellations; ++i) {
        for (int j = 0; j < 16; ++j) {
            symbols[16 * i + j] = gray_symbols_16qam[j];
        }
    }
    std::srand(time(0));
    std::random_shuffle(symbols.begin(), symbols.end());

    // signal forming
    Signal filter = rrc_filter(oversampling, rrc_roll_off, pulse_width);
    Signal field = convolution(symbols.upsample(oversampling), filter);
    double launch_power =
        dbm_to_watts(6 * std::cos((2.0 * point + 1.0) * M_PI / 42.0));
    double center_power = field.average_power(
        filter.size() / 2, field.size() - filter.size() / 2 - 1);
    field *= std::sqrt(launch_power / center_power);
    std::cout << "[ AVG = " << launch_power << ", PEAK = " << field.peak_power()
              << " ]" << std::endl;
    std::cout << ">> signal has formed" << std::endl;

    // propagation
    forward_propagation(field);
    Signal dbp_field = field;
    backward_propagation(dbp_field);
    cd_compensation(field);

    // signal matching
    field = convolution(field, filter);
    field = field.downsample(oversampling).chomp(pulse_width, pulse_width - 1);
    field *= std::sqrt(1.0 / field.average_power());
    dbp_field = convolution(dbp_field, filter);
    dbp_field =
        dbp_field.downsample(oversampling).chomp(pulse_width, pulse_width - 1);
    dbp_field *= std::sqrt(1.0 / dbp_field.average_power());
    std::cout << ">> signals has matched" << std::endl;

    // linear equalization
    LinearEqualizer leq;
    leq.train(symbols, field);
    Signal leq_field = leq.equalize(field);
    std::cout << ">> linear equalization has done" << std::endl;

    // saving results
    std::string filename = "results/" + std::to_string(point) + ".txt";
    std::ofstream fout(filename);
    fout.precision(15);

    fout << "# [ NUM, AVG, DBP, LEQ ] = " << point << " ";
    fout << launch_power << " ";
    fout << q2_factor(symbols, dbp_field) << " ";
    fout << q2_factor(symbols, leq_field) << std::endl;

    for (int i = 0; i < symbols.size(); ++i) {
        fout << symbols[i].real() << ",";
        fout << symbols[i].imag() << ",";
        fout << field[i].real() << ",";
        fout << field[i].imag() << ",";
        fout << dbp_field[i].real() << ",";
        fout << dbp_field[i].imag() << ",";
        fout << leq_field[i].real() << ",";
        fout << leq_field[i].imag() << "\n";
    }
    fout.close();
    std::cout << ">> results has saved" << std::endl;
}

int main() {
    for (int k = 0; k < 21; ++k)
        tx_rx(64, k);

    return 0;
}
