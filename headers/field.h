#ifndef FIELD_H_
#define FIELD_H_

#include <fftw3.h>
#include <cmath>
#include <complex>
#include <vector>

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

    Field upsample(const int& factor) const;
    Field downsample(const int& factor) const;
    Field chomp(const int& at_begin, const int& at_end) const;

    Field operator*(const Complex& multiplier) const;
    Field operator*(const Field& multipliers) const;
    Field& operator*=(const Complex& multiplier);
    Field& operator*=(const Field& multipliers);

    void setSampligRate(const double& rate);
    double getSamplingRate() const;
    double f(const int& i) const;
    double w(const int& i) const;
    RealVector temporal_power() const;
    RealVector spectral_power() const;

    Field& fft_inplace();
    Field& ifft_inplace();
    Field& fft_shift();
    Field& apply_filter(const Field& filter);
};

Field convolution(const Field& x, const Field& y);

#endif  // FIELD_H_
