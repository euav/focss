#include "../headers/field.h"

double Field::peak_power() const {
    double power = 0;
    for (int i = 0; i < size(); ++i)
        if (power < norm(at(i))) power = norm(at(i));

    return power;
}

double Field::average_power() const {
    double power = 0;
    for (int i = 0; i < size(); ++i)
        power += norm(at(i));

    return power / size();
}

Field Field::upsample(const int& factor) const {
    Field upsampled(factor * size(), 0);
    upsampled.setSampligRate(sampling_rate * factor);
    for (int i = 0; i < size(); ++i) {
        upsampled[factor * i] = at(i);
    }

    return upsampled;
}

Field Field::downsample(const int& factor) const {
    Field downsampled(size() / factor);
    downsampled.setSampligRate(sampling_rate / factor);
    for (int i = 0; i < downsampled.size(); ++i) {
        downsampled[i] = at(i * factor);
    }

    return downsampled;
}

Field Field::chomp(const int& at_begin, const int& at_end) const {
    Field chomped(size() - at_begin - at_end);
    for (int i = 0; i < chomped.size(); ++i) {
        chomped[i] = at(i + at_begin);
    }

    return chomped;
}

Field Field::operator*(const Complex& multiplier) const {
    Field copy(*this);
    for (int i = 0; i < copy.size(); ++i)
        copy[i] *= multiplier;
    return copy;
}

Field Field::operator*(const Field& multipliers) const {
    Field copy(*this);
    if (size() == multipliers.size())
        for (int i = 0; i < copy.size(); ++i)
            copy[i] *= multipliers[i];

    return copy;
}

Field& Field::operator*=(const Complex& multiplier) {
    for (int i = 0; i < size(); ++i)
        at(i) *= multiplier;
    return *this;
}

Field& Field::operator*=(const Field& multipliers) {
    if (size() == multipliers.size())
        for (int i = 0; i < size(); ++i)
            at(i) *= multipliers[i];
    return *this;
}

void Field::setSampligRate(const double& rate) {
    sampling_rate = rate;
    int samples = size();
    omega.assign(samples, 2 * M_PI * rate / samples);

    for (int i = 0; i <= samples / 2; ++i)
        omega[i] *= i;

    for (int i = samples / 2 + 1; i < samples; ++i) {
        omega[i] *= i - samples;
    }
}

double Field::getSamplingRate() const { return sampling_rate; }

double Field::f(const int& i) const { return omega[i] / (2 * M_PI); }

double Field::w(const int& i) const { return omega[i]; }

RealVector Field::temporal_power() const {
    RealVector power(size(), 0);
    for (int i = 0; i < size(); ++i)
        power[i] = norm(at(i));

    return power;
}

RealVector Field::spectral_power() const {
    RealVector power(size(), 0);
    Field copy = *this;
    copy.fft_inplace();
    for (int i = 0; i < size(); ++i)
        power[i] = norm(copy[i]);

    return power;
}

Field& Field::fft_inplace() {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(size(),
                         reinterpret_cast<fftw_complex*>(data()),
                         reinterpret_cast<fftw_complex*>(data()),
                         FFTW_FORWARD,
                         FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);

    for (int i = 0; i < size(); ++i)
        at(i) /= size();
    return *this;
}

Field& Field::ifft_inplace() {
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

Field& Field::fft_shift() {
    Complex buffer;
    int half_size = size() / 2;
    for (int i = 0; i < half_size; ++i) {
        buffer = at(i);
        at(i) = at(i + half_size);
        at(i + half_size) = buffer;
    }

    return *this;
}

Field& Field::apply_filter(const Field& filter) {
    Field fft_filter(size(), 0);
    int half_size = filter.size() / 2;
    int padding_shift = size() - half_size;
    for (int i = 0; i < filter.size() / 2; ++i) {
        fft_filter[i] = filter[i + half_size];
        fft_filter[padding_shift + i] = filter[i];
    }

    fft_inplace() *= fft_filter.fft_inplace();
    return ifft_inplace();
}

Field convolution(const Field& x, const Field& y) {
    Field z(x.size() + y.size() - 1);
    for (int i = 0; i < x.size(); ++i)
        for (int j = 0; j < y.size(); ++j)
            z[i + j] += x[i] * y[j];

    return z;
}
