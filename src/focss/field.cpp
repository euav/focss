#include "field.h"

double Field::peak_power() const {
    double power = 0;
    for (unsigned long i = 0; i < size(); ++i)
        if (power < norm(at(i))) power = norm(at(i));

    return power;
}

double Field::average_power() const {
    double power = 0;
    for (unsigned long i = 0; i < size(); ++i)
        power += norm(at(i));

    return power / size();
}

Field Field::upsample(const unsigned long& factor) const {
    Field upsampled(factor * size(), 0);
    upsampled.setSampligRate(sampling_rate * factor);
    for (unsigned long i = 0; i < size(); ++i) {
        upsampled[factor * i] = at(i);
    }

    return upsampled;
}

Field Field::downsample(const unsigned long& factor) const {
    Field downsampled(size() / factor);
    downsampled.setSampligRate(sampling_rate / factor);
    for (unsigned long i = 0; i < downsampled.size(); ++i) {
        downsampled[i] = at(i * factor);
    }

    return downsampled;
}

Field Field::decimate(const unsigned long& factor) const {
    Field lowpassed = this->fft();
    Field downsampled(size() / factor);
    downsampled.setSampligRate(sampling_rate / factor);

    unsigned long cutoff_index = lowpassed.size() / factor / 2;
    unsigned long shift_index = size() - cutoff_index;
    for (unsigned long i = 0; i < cutoff_index; ++i) {
        downsampled[i] = lowpassed[i];
        downsampled[cutoff_index + i] = lowpassed[shift_index + i];
    }
    downsampled.ifft_inplace();

    return downsampled;
}

Field Field::chomp(const unsigned long& at_begin,
                   const unsigned long& at_end) const {
    Field chomped(size() - at_begin - at_end);
    for (unsigned long i = 0; i < chomped.size(); ++i) {
        chomped[i] = at(i + at_begin);
    }

    return chomped;
}

Field Field::operator+(const Complex& summand) const {
    Field copy(*this);
    for (unsigned long i = 0; i < copy.size(); ++i)
        copy[i] += summand;
    return copy;
}

Field Field::operator+(const Field& summands) const {
    Field copy(*this);
    if (size() == summands.size())
        for (unsigned long i = 0; i < copy.size(); ++i)
            copy[i] += summands[i];
    else
        throw std::logic_error("fields size mismatch");

    return copy;
}

Field& Field::operator+=(const Complex& summand) {
    for (unsigned long i = 0; i < size(); ++i)
        at(i) += summand;
    return *this;
}

Field& Field::operator+=(const Field& summands) {
    if (size() == summands.size())
        for (unsigned long i = 0; i < size(); ++i)
            at(i) += summands[i];
    else
        throw std::logic_error("fields size mismatch");
    return *this;
}

Field Field::operator*(const Complex& multiplier) const {
    Field copy(*this);
    for (unsigned long i = 0; i < copy.size(); ++i)
        copy[i] *= multiplier;
    return copy;
}

Field Field::operator*(const Field& multipliers) const {
    Field copy(*this);
    if (size() == multipliers.size())
        for (unsigned long i = 0; i < copy.size(); ++i)
            copy[i] *= multipliers[i];
    else
        throw std::logic_error("fields size mismatch");

    return copy;
}

Field& Field::operator*=(const Complex& multiplier) {
    for (unsigned long i = 0; i < size(); ++i)
        at(i) *= multiplier;
    return *this;
}

Field& Field::operator*=(const Field& multipliers) {
    if (size() == multipliers.size())
        for (unsigned long i = 0; i < size(); ++i)
            at(i) *= multipliers[i];
    else
        throw std::logic_error("fields size mismatch");

    return *this;
}

void Field::setSampligRate(const double& rate) {
    sampling_rate = rate;
    unsigned long samples = size();
    omega.assign(samples, 2 * M_PI * rate / samples);

    for (unsigned long i = 0; i <= samples / 2; ++i)
        omega[i] *= double(i);

    for (unsigned long i = samples / 2 + 1; i < samples; ++i) {
        omega[i] *= double(i) - double(samples);
    }
}

double Field::getSamplingRate() const { return sampling_rate; }

double Field::dt() const { return 1 / sampling_rate; }

double Field::df() const { return sampling_rate / size(); }

double Field::dw() const { return 2 * M_PI * sampling_rate / size(); }

double Field::f(const unsigned long& i) const { return omega[i] / (2 * M_PI); }

double Field::w(const unsigned long& i) const { return omega[i]; }

RealVector Field::temporal_power() const {
    RealVector power(size(), 0);
    for (unsigned long i = 0; i < size(); ++i)
        power[i] = norm(at(i));

    return power;
}

RealVector Field::spectral_power() const {
    RealVector power(size(), 0);
    Field copy = *this;
    copy.fft_inplace();
    for (unsigned long i = 0; i < size(); ++i)
        power[i] = norm(copy[i]);

    return power;
}

Field Field::fft() const {
    Field transformed = *this;
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(int(size()),
                         reinterpret_cast<fftw_complex*>(transformed.data()),
                         reinterpret_cast<fftw_complex*>(transformed.data()),
                         FFTW_FORWARD,
                         FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);

    transformed *= 1 / double(size());
    return transformed;
}

Field Field::ifft() const {
    Field transformed = *this;
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(int(size()),
                         reinterpret_cast<fftw_complex*>(transformed.data()),
                         reinterpret_cast<fftw_complex*>(transformed.data()),
                         FFTW_BACKWARD,
                         FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);
    return transformed;
}

Field& Field::fft_inplace() {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(int(size()),
                         reinterpret_cast<fftw_complex*>(data()),
                         reinterpret_cast<fftw_complex*>(data()),
                         FFTW_FORWARD,
                         FFTW_ESTIMATE);
    fftw_execute(complex_inplace);
    fftw_destroy_plan(complex_inplace);

    for (unsigned long i = 0; i < size(); ++i)
        at(i) /= size();
    return *this;
}

Field& Field::ifft_inplace() {
    fftw_plan complex_inplace =
        fftw_plan_dft_1d(int(size()),
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
    unsigned long half_size = size() / 2;
    for (unsigned long i = 0; i < half_size; ++i) {
        buffer = at(i);
        at(i) = at(i + half_size);
        at(i + half_size) = buffer;
    }

    return *this;
}

Field& Field::apply_filter(const Field& filter) {
    Field fft_filter(size(), 0);
    unsigned long half_size = filter.size() / 2;
    unsigned long shift_index = size() - half_size;
    for (unsigned long i = 0; i < filter.size() / 2; ++i) {
        fft_filter[i] = filter[i + half_size];
        fft_filter[shift_index + i] = filter[i];
    }

    fft_inplace() *= fft_filter.fft_inplace();
    return ifft_inplace();
}

Polarizations Polarizations::operator*(const Complex& multiplier) const {
    return Polarizations{x * multiplier, y * multiplier};
}

Polarizations Polarizations::operator*(const Field& multipliers) const {
    return Polarizations{x * multipliers, y * multipliers};
}

Polarizations& Polarizations::operator*=(const Complex& multiplier) {
    x *= multiplier;
    y *= multiplier;

    return *this;
}

Polarizations& Polarizations::operator*=(const Field& multipliers) {
    x *= multipliers;
    y *= multipliers;

    return *this;
}

Field convolution(const Field& x, const Field& y) {
    Field z(x.size() + y.size() - 1);
    for (unsigned long i = 0; i < x.size(); ++i)
        for (unsigned long j = 0; j < y.size(); ++j)
            z[i + j] += x[i] * y[j];

    return z;
}
