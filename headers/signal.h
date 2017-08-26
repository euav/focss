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

    double peak_power() const;
    double average_power() const;
    
    Signal upsample(const int& factor) const;
    Signal downsample(const int& factor) const;
    Signal chomp(const int& at_begin, const int& at_end) const;
    //TODO: Signal& upsample_inplace(const int& factor);
    //TODO: Signal& downsample_inplace(const int& factor);
    //TODO: Signal& chomp_inplace(const int& at_begin, const int& at_end);

    Signal operator*(const Complex& multiplier) const;
    Signal operator*(const Signal& multipliers) const;
    Signal& operator*=(const Complex& multiplier);
    Signal& operator*=(const Signal& multipliers);

    Signal& fft_inplace();
    Signal& ifft_inplace();
    Signal& fft_shift();
    Signal& apply_filter(const Signal& filter);
};

Signal convolution(const Signal& x, const Signal& y);

#endif  // SIGNAL_H_
