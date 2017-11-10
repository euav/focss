#include "functions.h"
#include <fftw3.h>
#include <cmath>
#include <complex>
#include "focss/vector.h"

namespace focss {
void relax_max(double& a, const double& b) {
    if (a < b) a = b;
}

void relax_mix(double& a, const double& b) {
    if (a > b) a = b;
}

double sinc(const double& x) {
    if (x == 0) return 1;
    return std::sin(math_pi * x) / math_pi / x;
}

double evm2_factor(const ComplexVector& tx, const ComplexVector& rx) {
    double numerator = 0;
    double denominator = 0;

    unsigned long transmission_size = std::min(tx.size(), rx.size());
    for (unsigned long i = 0; i < transmission_size; ++i) {
        numerator += norm(tx[i] - rx[i]);
        denominator += norm(tx[i]);
    }

    return numerator / denominator;
}

double q2_factor(const ComplexVector& tx, const ComplexVector& rx) {
    return -10 * std::log10(evm2_factor(tx, rx));
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

Complex i_exp(const double& x) { return Complex(std::cos(x), std::sin(x)); }

// ----------------------------------------------------------------------
// -------------------- Complex* fast fourier transform (fft)
// ----------------------------------------------------------------------
void fft_inplace(const int& size, Complex* input_output) {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(size,
                         reinterpret_cast<fftw_complex*>(input_output),
                         reinterpret_cast<fftw_complex*>(input_output),
                         FFTW_FORWARD,
                         FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);

    for (int i = 0; i < size; ++i)
        input_output[i] /= size;
}

void ifft_inplace(const int& size, Complex* input_output) {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(size,
                         reinterpret_cast<fftw_complex*>(input_output),
                         reinterpret_cast<fftw_complex*>(input_output),
                         FFTW_BACKWARD,
                         FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);
}

void fft(const int& size, Complex* input, Complex* output) {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(size,
                         reinterpret_cast<fftw_complex*>(input),
                         reinterpret_cast<fftw_complex*>(output),
                         FFTW_FORWARD,
                         FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);

    for (int i = 0; i < size; ++i)
        output[i] /= size;
}

void ifft(const int& size, Complex* input, Complex* output) {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(size,
                         reinterpret_cast<fftw_complex*>(input),
                         reinterpret_cast<fftw_complex*>(output),
                         FFTW_FORWARD,
                         FFTW_ESTIMATE);
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

// void save_transmission(const ComplexVector& tx,
//                        const ComplexVector& rx,
//                        const char* filename) {
//     std::ofstream file(filename);
//     file.precision(16);
//     unsigned long n_symbols = std::min(tx.size(), rx.size());
//     for (unsigned long i = 0; i < n_symbols; ++i) {
//         file << tx[i].real() << ',' << tx[i].imag() << ',';
//         file << rx[i].real() << ',' << rx[i].imag() << '\n';
//     }
//     file.close();
// }

// void load_transmission(ComplexVector* tx,
//                        ComplexVector* rx,
//                        const char* filename) {
//     std::ifstream file(filename);
//     char delimeter = ',';
//     double tx_re, tx_im, rx_re, rx_im;
//     while (file.good()) {
//         file >> tx_re >> delimeter >> tx_im >> delimeter;
//         file >> rx_re >> delimeter >> rx_im >> std::ws;

//         tx->push_back(Complex(tx_re, tx_im));
//         rx->push_back(Complex(rx_re, rx_im));
//     }
//     file.close();
// }
}  // namespace focss
