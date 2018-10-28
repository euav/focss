#include "field.h"
#include <cassert>
#include <cmath>
#include "focss/definitions.h"
#include "focss/functions.h"

namespace focss {
// ----------------------------------------------------------------------
// -------------------- Field constructors and destructor
// ----------------------------------------------------------------------
Field::Field() {
    modes_ = 0;
    samples_ = 0;
    time_step_ = 0;
    center_frequency_ = 0;

    size_ = 0;
    data_ = nullptr;

    forward_inplace_ = nullptr;
    backward_inplace_ = nullptr;
}

Field::Field(const Domain& domain) {
    assert(domain.modes > 0 && domain.samples > 0 && domain.time_step > 0);

    modes_ = domain.modes;
    samples_ = domain.samples;
    time_step_ = domain.time_step;
    center_frequency_ = domain.center_frequency;

    size_ = modes_ * samples_;
    data_ = reinterpret_cast<complex_t*>(fftw_alloc_complex(size_));
    prepare_fft_plans();

    for (int index = 0; index < size_; ++index)
        data_[index] = complex_t(0.0, 0.0);
}

Field::Field(const Domain& domain, const ComplexVector& data) {
    assert(domain.size() == data.size() && domain.time_step > 0);

    modes_ = domain.modes;
    samples_ = domain.samples;
    time_step_ = domain.time_step;
    center_frequency_ = domain.center_frequency;

    size_ = modes_ * samples_;
    data_ = reinterpret_cast<complex_t*>(fftw_alloc_complex(size_));
    prepare_fft_plans();

    for (int index = 0; index < size_; ++index)
        data_[index] = data[index];
}

Field::Field(const Field& other) {
    modes_ = other.modes_;
    samples_ = other.samples_;
    time_step_ = other.time_step_;
    center_frequency_ = other.center_frequency_;

    size_ = modes_ * samples_;
    data_ = reinterpret_cast<complex_t*>(fftw_alloc_complex(size_));
    prepare_fft_plans();

    for (int index = 0; index < size_; ++index)
        data_[index] = other.data_[index];
}

Field::Field(Field&& other) {
    modes_ = other.modes_;
    samples_ = other.samples_;
    time_step_ = other.time_step_;
    center_frequency_ = other.center_frequency_;

    size_ = other.size_;
    data_ = other.data_;

    forward_inplace_ = other.forward_inplace_;
    backward_inplace_ = other.backward_inplace_;

    other.data_ = nullptr;
    other.forward_inplace_ = nullptr;
    other.backward_inplace_ = nullptr;
}

Field& Field::operator=(const Field& other) {
    if (this != &other) {
        release_resources();

        modes_ = other.modes_;
        samples_ = other.samples_;
        time_step_ = other.time_step_;
        center_frequency_ = other.center_frequency_;

        size_ = modes_ * samples_;
        data_ = reinterpret_cast<complex_t*>(fftw_alloc_complex(size_));
        prepare_fft_plans();

        for (int index = 0; index < size_; ++index)
            data_[index] = other.data_[index];
    }

    return *this;
}

Field& Field::operator=(Field&& other) {
    if (this != &other) {
        release_resources();

        modes_ = other.modes_;
        samples_ = other.samples_;
        time_step_ = other.time_step_;
        center_frequency_ = other.center_frequency_;

        size_ = other.size_;
        data_ = other.data_;

        forward_inplace_ = other.forward_inplace_;
        backward_inplace_ = other.backward_inplace_;

        other.data_ = nullptr;
        other.forward_inplace_ = nullptr;
        other.backward_inplace_ = nullptr;
    }

    return *this;
}

Field::~Field() { release_resources(); }

// ----------------------------------------------------------------------
// -------------------- Field accessors
// ----------------------------------------------------------------------
complex_t Field::operator[](const int& index) const {
    assert(0 <= index && index < size_);
    return data_[index];
}

complex_t Field::operator()(const int& mode, const int& sample) const {
    assert(0 <= mode && mode < modes_);
    assert(0 <= sample && sample < samples_);
    return data_[sample + mode * samples_];
}

complex_t& Field::operator[](const int& index) {
    assert(0 <= index && index < size_);
    return data_[index];
}

complex_t& Field::operator()(const int& mode, const int& sample) {
    assert(0 <= mode && mode < modes_);
    assert(0 <= sample && sample < samples_);
    return data_[sample + mode * samples_];
}

ComplexVector Field::operator()(const int& mode) {
    return ComplexVector::proxy(samples_, data_ + samples_ * mode);
}

// ----------------------------------------------------------------------
// -------------------- Field power methods
// ----------------------------------------------------------------------
real_t Field::power(const int& sample) const {
    assert(0 <= sample && sample < samples_);
    real_t power = 0;
    for (int index = sample; index < size_; index += samples_)
        power += norm(data_[index]);

    return power;
}

real_t Field::peak_power() const {
    real_t peak_power = 0;
    for (int i = 0; i < samples_; ++i)
        if (peak_power < power(i)) peak_power = power(i);

    return peak_power;
}

real_t Field::average_power() const {
    real_t average_power = 0;
    for (int index = 0; index < size_; ++index)
        average_power += norm(data_[index]);

    return average_power / samples_;
}

Field& Field::peak_normalize(const real_t& power) {
    return (*this) *= std::sqrt(power / peak_power());
}

Field& Field::average_normalize(const real_t& power) {
    return (*this) *= std::sqrt(power / average_power());
}

// ----------------------------------------------------------------------
// -------------------- Field shape and domain methods
// ----------------------------------------------------------------------
int Field::size() const { return size_; }

int Field::modes() const { return modes_; }

int Field::samples() const { return samples_; }

real_t Field::duration() const { return samples_ * time_step_; }

real_t Field::bandwidth() const { return 1.0 / time_step_; }

real_t Field::time_step() const { return time_step_; }

real_t Field::sampling_rate() const { return 1.0 / time_step_; }

real_t Field::center_frequency() const { return center_frequency_; }

real_t Field::dt() const { return time_step_; }

real_t Field::df() const { return 1.0 / (time_step_ * samples_); }

real_t Field::dw() const { return 2.0 * math_pi / (time_step_ * samples_); }

real_t Field::t(const int& sample) const {
    return static_cast<real_t>(sample) * dt();
}

real_t Field::f(const int& sample) const {
    return center_frequency_ + (sample - samples_ / 2) * df();
}

real_t Field::w(const int& sample) const { return 2 * math_pi * f(sample); }

Domain Field::domain() const {
    return {modes_, samples_, time_step_, center_frequency_};
}

RealVector Field::time_grid() const {
    RealVector time_grid(samples_);
    for (int sample = 0; sample < samples_; ++sample)
        time_grid[sample] = t(sample);

    return time_grid;
}

RealVector Field::frequency_grid() const {
    RealVector frequency_grid(samples_);
    for (int sample = 0; sample < samples_; ++sample)
        frequency_grid[sample] = f(sample);

    return frequency_grid;
}

RealVector Field::temporal_power() const {
    RealVector temporal_power(samples_);
    for (int sample = 0; sample < samples_; ++sample)
        temporal_power[sample] = power(sample);

    return temporal_power;
}

RealVector Field::spectral_power() const {
    RealVector spectral_power(samples_);
    Field spectrum = this->fft().fft_shift();
    for (int sample = 0; sample < samples_; ++sample)
        spectral_power[sample] = spectrum.power(sample);

    return spectral_power;
}

// ----------------------------------------------------------------------
// -------------------- Field resampling methods
// ----------------------------------------------------------------------
Field Field::upsample(const int& factor) const {
    assert(factor > 0);

    Domain upsampled_domain = domain();
    upsampled_domain.samples *= factor;
    upsampled_domain.time_step /= factor;
    Field upsampled(upsampled_domain);
    for (int index = 0; index < size_; ++index)
        upsampled.data_[factor * index] = data_[index];

    return upsampled;
}

Field Field::downsample(const int& factor) const {
    assert(factor > 0);
    assert(samples_ % factor == 0);

    Domain downsampled_domain = domain();
    downsampled_domain.samples /= factor;
    downsampled_domain.time_step *= factor;
    Field downsampled(downsampled_domain);
    for (int index = 0; index < downsampled.size(); ++index)
        downsampled.data_[index] = data_[index * factor];

    return downsampled;
}

Field Field::decimate(const int& factor) const {
    assert(factor > 0);
    assert(samples_ % factor == 0);

    Domain decimated_domain = domain();
    decimated_domain.samples /= factor;
    decimated_domain.time_step *= factor;
    Field decimated(decimated_domain);

    Field spectrum = this->fft();
    int cutoff = spectrum.samples_ / factor / 2;
    int shift = spectrum.samples_ - decimated.samples_;

    for (int mode = 0; mode < modes_; ++mode) {
        for (int sample = 0; sample < cutoff; ++sample)
            decimated(mode, sample) = spectrum(mode, sample);
        for (int sample = cutoff; sample < decimated.samples_; ++sample)
            decimated(mode, sample) = spectrum(mode, sample + shift);
    }

    return decimated.ifft_inplace();
}

// ----------------------------------------------------------------------
// -------------------- Field fourier methods (fft)
// ----------------------------------------------------------------------
Field& Field::fft_shift() {
    int left, right;
    int half_samples = samples_ / 2;
    complex_t swap_buffer;
    for (int mode = 0; mode < modes_; ++mode) {
        for (int sample = 0; sample < half_samples; ++sample) {
            left = sample + mode * samples_;
            right = sample + mode * samples_ + half_samples;

            swap_buffer = data_[left];
            data_[left] = data_[right];
            data_[right] = swap_buffer;
        }
    }

    return *this;
}

Field& Field::fft_inplace() {
    fftw_execute(forward_inplace_);
    (*this) *= (1.0 / static_cast<real_t>(samples_));
    return *this;
}

Field& Field::ifft_inplace() {
    fftw_execute(backward_inplace_);
    return *this;
}

Field Field::fft() const { return Field(*this).fft_inplace(); }

Field Field::ifft() const { return Field(*this).ifft_inplace(); }

Field& Field::filter_inplace(const ComplexVector& impulse_response) {
    ComplexVector zero_padded(samples_, 0);
    int half_size = impulse_response.size() / 2;
    int shift_index = samples_ - half_size;
    for (int i = 0; i < impulse_response.size() / 2; ++i) {
        zero_padded[i] = impulse_response[i + half_size];
        zero_padded[shift_index + i] = impulse_response[i];
    }

    fft_inplace();
    (*this) *= focss::fft(zero_padded);
    ifft_inplace();

    return *this;
}

Field Field::filter(const ComplexVector& impulse_response) const {
    return Field(*this).filter_inplace(impulse_response);
}

// ----------------------------------------------------------------------
// -------------------- Field multiplication
// ----------------------------------------------------------------------
Field& Field::operator*=(const real_t& multiplier) {
    for (int index = 0; index < size_; ++index)
        data_[index] *= multiplier;
    return *this;
}

Field& Field::operator*=(const complex_t& multiplier) {
    for (int index = 0; index < size_; ++index)
        data_[index] *= multiplier;
    return *this;
}

Field& Field::operator*=(const RealVector& multipliers) {
    assert(multipliers.size() == samples_);
    for (int mode = 0; mode < modes_; ++mode)
        for (int sample = 0; sample < samples_; ++sample)
            data_[sample + mode * samples_] *= multipliers[sample];
    return *this;
}

Field& Field::operator*=(const ComplexVector& multipliers) {
    assert(multipliers.size() == samples_);
    for (int mode = 0; mode < modes_; ++mode)
        for (int sample = 0; sample < samples_; ++sample)
            data_[sample + mode * samples_] *= multipliers[sample];
    return *this;
}

Field operator*(const Field& field, const real_t& operand) {
    return Field(field) *= operand;
}

Field operator*(const Field& field, const complex_t& operand) {
    return Field(field) *= operand;
}

Field operator*(const Field& field, const RealVector& operand) {
    return Field(field) *= operand;
}

Field operator*(const Field& field, const ComplexVector& operand) {
    return Field(field) *= operand;
}

Field operator*(const real_t& operand, const Field& field) {
    return Field(field) *= operand;
}

Field operator*(const complex_t& operand, const Field& field) {
    return Field(field) *= operand;
}

Field operator*(const RealVector& operand, const Field& field) {
    return Field(field) *= operand;
}

Field operator*(const ComplexVector& operand, const Field& field) {
    return Field(field) *= operand;
}

// ----------------------------------------------------------------------
// -------------------- Field division
// ----------------------------------------------------------------------
Field& Field::operator/=(const real_t& multiplier) {
    for (int index = 0; index < size_; ++index)
        data_[index] /= multiplier;
    return *this;
}

Field& Field::operator/=(const complex_t& multiplier) {
    for (int index = 0; index < size_; ++index)
        data_[index] /= multiplier;
    return *this;
}

Field& Field::operator/=(const RealVector& multipliers) {
    assert(multipliers.size() == samples_);
    for (int mode = 0; mode < modes_; ++mode)
        for (int sample = 0; sample < samples_; ++sample)
            data_[sample + mode * samples_] /= multipliers[sample];
    return *this;
}

Field& Field::operator/=(const ComplexVector& multipliers) {
    assert(multipliers.size() == samples_);
    for (int mode = 0; mode < modes_; ++mode)
        for (int sample = 0; sample < samples_; ++sample)
            data_[sample + mode * samples_] /= multipliers[sample];
    return *this;
}

Field operator/(const Field& field, const real_t& operand) {
    return Field(field) /= operand;
}

Field operator/(const Field& field, const complex_t& operand) {
    return Field(field) /= operand;
}

Field operator/(const Field& field, const RealVector& operand) {
    return Field(field) /= operand;
}

Field operator/(const Field& field, const ComplexVector& operand) {
    return Field(field) /= operand;
}

Field operator/(const real_t& operand, const Field& field) {
    return Field(field) /= operand;
}

Field operator/(const complex_t& operand, const Field& field) {
    return Field(field) /= operand;
}

Field operator/(const RealVector& operand, const Field& field) {
    return Field(field) /= operand;
}

Field operator/(const ComplexVector& operand, const Field& field) {
    return Field(field) /= operand;
}

// ----------------------------------------------------------------------
// -------------------- Field addition/substraction
// ----------------------------------------------------------------------
Field& Field::operator+=(const Field& other) {
    for (int index = 0; index < size_; ++index)
        data_[index] += other.data_[index];
    return *this;
}

Field& Field::operator-=(const Field& other) {
    for (int index = 0; index < size_; ++index)
        data_[index] -= other.data_[index];
    return *this;
}

Field operator+(const Field& lhs, const Field& rhs) {
    return Field(lhs) += rhs;
}

Field operator-(const Field& lhs, const Field& rhs) {
    return Field(lhs) -= rhs;
}

// ----------------------------------------------------------------------
// -------------------- Field private methods
// ----------------------------------------------------------------------
void Field::prepare_fft_plans() {
    int dimensions[] = {samples_};
    fftw_complex* fftw_data = reinterpret_cast<fftw_complex*>(data_);

    forward_inplace_ = fftw_plan_many_dft(
        1, dimensions, modes_, fftw_data, nullptr, 1, samples_, fftw_data,
        nullptr, 1, samples_, FFTW_FORWARD, FFTW_MEASURE);

    backward_inplace_ = fftw_plan_many_dft(
        1, dimensions, modes_, fftw_data, nullptr, 1, samples_, fftw_data,
        nullptr, 1, samples_, FFTW_BACKWARD, FFTW_MEASURE);
}

void Field::release_resources() {
    fftw_free(reinterpret_cast<fftw_complex*>(data_));
    fftw_destroy_plan(forward_inplace_);
    fftw_destroy_plan(backward_inplace_);
}
}  // namespace focss
