#ifndef FIELD_H_
#define FIELD_H_

#include <fftw3.h>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <stdexcept>

typedef std::complex<double> Complex;
typedef std::vector<Complex> ComplexVector;
typedef std::vector<double> RealVector;
const Complex i_unit = Complex(0, 1);

class Field : public ComplexVector {
    double sampling_rate;
    RealVector omega;

  public:
    using ComplexVector::vector;
    // Field();
    // Field(const int& size);
    // Field(const int& size, const int& default_value);

    double peak_power() const;
    double average_power() const;

    Field upsample(const unsigned long& factor) const;
    Field downsample(const unsigned long& factor) const;
    Field chomp(const unsigned long& at_begin,
                const unsigned long& at_end) const;
    Field decimate(const unsigned long& factor) const;

    Field operator+(const Complex& summand) const;
    Field operator+(const Field& summands) const;
    Field& operator+=(const Complex& summand);
    Field& operator+=(const Field& summands);
    
    Field operator*(const Complex& multiplier) const;
    Field operator*(const Field& multipliers) const;
    Field& operator*=(const Complex& multiplier);
    Field& operator*=(const Field& multipliers);

    void setSampligRate(const double& rate);
    double getSamplingRate() const;
    double dt() const;
    double df() const;
    double dw() const;
    double f(const unsigned long& i) const;
    double w(const unsigned long& i) const;
    RealVector temporal_power() const;
    RealVector spectral_power() const;

    Field fft() const;
    Field ifft() const;
    Field& fft_inplace();
    Field& ifft_inplace();
    Field& fft_shift();
    Field& apply_filter(const Field& filter);
};

struct Polarizations {
    Field x, y;

    Polarizations operator*(const Complex& multiplier) const;
    Polarizations operator*(const Field& multipliers) const;
    Polarizations& operator*=(const Complex& multiplier);
    Polarizations& operator*=(const Field& multipliers);
};

Field convolution(const Field& x, const Field& y);

#endif  // FIELD_H_
