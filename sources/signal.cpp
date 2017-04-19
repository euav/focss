#include "../headers/signal.h"

double Signal::peak_power() const {
    double power = 0;
    for (int i = 0; i < size(); ++i)
        if (power < norm(at(i))) power = norm(at(i));

    return power;
}

double Signal::average_power() const {
    double power = 0;
    for (int i = 0; i < size(); ++i) power += norm(at(i));

    return power / size();
}

double Signal::peak_power(const int& from, const int& to) const {
    double power = 0;
    for (int i = from; i < to; ++i)
        if (power < norm(at(i))) power = norm(at(i));

    return power;
}

double Signal::average_power(const int& from, const int& to) const {
    double power = 0;
    for (int i = from; i < to; ++i) power += norm(at(i));

    return power / (to - from);
}

Signal& Signal::fft_inplace() {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(size(),
                         reinterpret_cast<fftw_complex*>(data()),
                         reinterpret_cast<fftw_complex*>(data()),
                         FFTW_FORWARD,
                         FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);

    for (int i = 0; i < size(); ++i) at(i) /= size();
    return *this;
}

Signal& Signal::ifft_inplace() {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(size(),
                         reinterpret_cast<fftw_complex*>(data()),
                         reinterpret_cast<fftw_complex*>(data()),
                         FFTW_BACKWARD,
                         FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);
    return *this;
}

Signal& Signal::operator*=(const Complex& multiplier) {
    for (int i = 0; i < size(); ++i) at(i) *= multiplier;
    return *this;
}

Signal& Signal::operator*=(const Signal& multipliers) {
    if (size() == multipliers.size())
        for (int i = 0; i < size(); ++i) at(i) *= multipliers[i];
    return *this;
}

Signal Signal::upsample(const int& factor) {
    Signal upsampled(factor * size(), 0);
    for (int i = 0; i < size(); ++i) {
        upsampled[factor * i] = at(i);
    }

    return upsampled;
}

Signal Signal::downsample(const int& factor) {
    Signal downsampled(size() / factor);
    for (int i = 0; i < downsampled.size(); ++i) {
        downsampled[i] = at(i * factor);
    }

    return downsampled;
}

Signal Signal::chomp(const int& at_begin, const int& at_end) {
    Signal chomped(size() - at_begin - at_end);
    for (int i = 0; i < chomped.size(); ++i) {
        chomped[i] = at(i + at_begin);
    }

    return chomped;
}

Signal convolution(const Signal& x, const Signal& y) {
    Signal z(x.size() + y.size() - 1);
    for (int i = 0; i < x.size(); ++i)
        for (int j = 0; j < y.size(); ++j) z[i + j] += x[i] * y[j];

    return z;
}
