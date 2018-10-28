#include "functions.h"
#include <fftw3.h>
#include <cassert>
#include <cmath>
#include <fstream>
#include <random>
#include "focss/definitions.h"
#include "focss/field.h"

namespace focss {
void setup() {
    // fftw_init_threads();
    // fftw_plan_with_nthreads(2);
    focss::reseed_global_urng();
}

void reseed_global_urng() {
    static std::random_device rdev;
    global_urng().seed(rdev());
}

std::default_random_engine& global_urng() {
    static std::default_random_engine urng;
    return urng;
}

double evm2_factor(const ComplexVector& tx, const ComplexVector& rx) {
    assert(tx.size() == rx.size());

    double numerator = 0;
    double denominator = 0;
    for (int i = 0; i < tx.size(); ++i) {
        numerator += norm(tx[i] - rx[i]);
        denominator += norm(tx[i]);
    }

    return numerator / denominator;
}

double evm2_factor(const Field& tx, const Field& rx) {
    assert(tx.size() == rx.size());

    double numerator = 0;
    double denominator = 0;
    for (int i = 0; i < tx.size(); ++i) {
        numerator += norm(tx[i] - rx[i]);
        denominator += norm(tx[i]);
    }

    return numerator / denominator;
}

double evm2_factor(const Field& tx, const Field& rx, const int& cut) {
    assert(tx.modes() == rx.modes());
    assert(tx.samples() == rx.samples());

    double numerator = 0;
    double denominator = 0;
    for (int mode = 0; mode < tx.modes(); ++mode) {
        for (int sample = cut; sample < tx.samples() - cut; ++sample) {
            numerator += norm(tx(mode, sample) - rx(mode, sample));
            denominator += norm(tx(mode, sample));
        }
    }

    return numerator / denominator;
}

double q2_factor(const ComplexVector& tx, const ComplexVector& rx) {
    return -10 * std::log10(evm2_factor(tx, rx));
}

double q2_factor(const Field& tx, const Field& rx) {
    return -10 * std::log10(evm2_factor(tx, rx));
}

double q2_factor(const Field& tx, const Field& rx, const int& cut) {
    return -10 * std::log10(evm2_factor(tx, rx, cut));
    // return  3.0 / 8.0 * std::erfc(sqrt(0.1 / evm2_factor(tx, rx, cut)));
}

double db_to_linear(const double& db_value) {
    return std::pow(10, db_value / 10);
}

double dbm_to_watts(const double& dbm_power) {
    return 1e-3 * std::pow(10, dbm_power / 10);
}

double db_to_natural(const double& db_value) {
    return db_value * std::log(10) / 10;
}

double disp_to_beta2(const double& dispersion, const double& wavelength) {
    return -wavelength * wavelength * dispersion / (2 * math_pi * light_speed);
}

double sinc(const double& x) {
    if (x == 0) return 1;
    return std::sin(math_pi * x) / math_pi / x;
}

complex_t i_exp(const double& x) { return complex_t(std::cos(x), std::sin(x)); }

// ----------------------------------------------------------------------
// -------------------- complex_t* fast fourier transform (fft)
// ----------------------------------------------------------------------
void fft_inplace(const int& size, complex_t* input_output) {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(input_output),
                         reinterpret_cast<fftw_complex*>(input_output),
                         FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);

    for (int i = 0; i < size; ++i)
        input_output[i] /= size;
}

void ifft_inplace(const int& size, complex_t* input_output) {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(input_output),
                         reinterpret_cast<fftw_complex*>(input_output),
                         FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);
}

void fft(const int& size, complex_t* input, complex_t* output) {
    fftw_plan complex_inplace = fftw_plan_dft_1d(
        size, reinterpret_cast<fftw_complex*>(input),
        reinterpret_cast<fftw_complex*>(output), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);

    for (int i = 0; i < size; ++i)
        output[i] /= size;
}

void ifft(const int& size, complex_t* input, complex_t* output) {
    fftw_plan complex_inplace = fftw_plan_dft_1d(
        size, reinterpret_cast<fftw_complex*>(input),
        reinterpret_cast<fftw_complex*>(output), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);
}

// ----------------------------------------------------------------------
// -------------------- ComplexVector fast fourier transform (fft)
// ----------------------------------------------------------------------
void fft_inplace(ComplexVector* data) {
    fft_inplace(data->size(), data->raw());
}

void ifft_inplace(ComplexVector* data) {
    ifft_inplace(data->size(), data->raw());
}

ComplexVector fft(const ComplexVector& data) {
    ComplexVector copy(data);
    fft_inplace(&copy);
    return copy;
}

ComplexVector ifft(const ComplexVector& data) {
    ComplexVector copy(data);
    ifft_inplace(&copy);
    return copy;
}

// ----------------------------------------------------------------------
// -------------------- File import/export of Field
// ----------------------------------------------------------------------
void save_transmission(std::ostream& output, const Field& tx, const Field& rx) {
    assert(tx.modes() == rx.modes());
    assert(tx.samples() == rx.samples());

    for (int i = 0; i < tx.samples(); ++i) {
        for (int mode = 0; mode < tx.modes(); ++mode)
            output << tx(mode, i).real() << ',' << tx(mode, i).imag() << ',';
        for (int mode = 0; mode < rx.modes() - 1; ++mode)
            output << rx(mode, i).real() << ',' << rx(mode, i).imag() << ',';

        output << rx(rx.modes() - 1, i).real() << ',';
        output << rx(rx.modes() - 1, i).imag() << '\n';
    }
}

void load_transmission(std::istream& input, Field* tx, Field* rx) {
    assert(tx->modes() == rx->modes());
    assert(tx->samples() == rx->samples());

    double re, im;
    char delimeter = ',';
    for (int i = 0; i < tx->samples(); ++i) {
        input >> re >> delimeter >> im;
        tx->operator()(0, i) = complex_t(re, im);

        for (int mode = 1; mode < tx->modes(); ++mode) {
            input >> delimeter >> re >> delimeter >> im;
            tx->operator()(mode, i) = complex_t(re, im);
        }
        for (int mode = 0; mode < tx->modes(); ++mode) {
            input >> delimeter >> re >> delimeter >> im;
            rx->operator()(mode, i) = complex_t(re, im);
        }
    }
}

void save_transmission(const std::string& filename, const Field& tx, const Field& rx) {
    std::ofstream file(filename);
    file.precision(15);
    save_transmission(file, tx, rx);
    file.close();
}

void load_transmission(const std::string& filename, Field* tx, Field* rx) {
    std::ifstream file(filename);
    load_transmission(file, tx, rx);
    file.close();
}
}  // namespace focss
