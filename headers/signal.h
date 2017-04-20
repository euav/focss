#ifndef SIGNAL_H_
#define SIGNAL_H_

#include <fftw3.h>
#include <cmath>
#include <complex>
#include <vector>

typedef std::complex<double> Complex;
typedef std::vector<Complex> ComplexVector;
const Complex i_unit = Complex(0, 1);

class Signal : public ComplexVector {
  public:
    using ComplexVector::vector;

	double null_power() const;
    double peak_power() const;
    double average_power() const;
    double peak_power(const int& from, const int& to) const;
    double average_power(const int& from, const int& to) const;

    Signal& fft_inplace();
    Signal& ifft_inplace();

    Signal upsample(const int& factor);
    Signal downsample(const int& factor);
    Signal chomp(const int& at_begin, const int& at_end);
    Signal& operator*=(const Complex& multiplier);
    Signal& operator*=(const Signal& multipliers);
};

Signal convolution(const Signal& x, const Signal& y);

#endif  // SIGNAL_H_
