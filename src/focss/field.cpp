#include "field.h"
#include <cassert>
#include <cmath>
#include "focss/core.h"
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
    circular_frequencies_ = nullptr;

    size_ = 0;
    data_ = nullptr;

    forward_inplace_ = nullptr;
    backward_inplace_ = nullptr;
}

Field::Field(const int& samples) : Field(1, samples) {}

Field::Field(const int& modes, const int& samples) {
    assert(samples > 0 && modes > 0);

    modes_ = modes;
    samples_ = samples;
    time_step_ = 0;
    center_frequency_ = 0;
    circular_frequencies_ = nullptr;

    size_ = modes_ * samples_;
    data_ = reinterpret_cast<Complex*>(fftw_alloc_complex(size_));
    for (int i = 0; i < size_; ++i)
        data_[i] = 0;

    calculate_fftw_plans();
}

Field::Field(const ComplexVector& data) {
    modes_ = 1;
    samples_ = data.size();
    time_step_ = 0;
    center_frequency_ = 0;
    circular_frequencies_ = nullptr;

    size_ = modes_ * samples_;
    data_ = reinterpret_cast<Complex*>(fftw_alloc_complex(size_));
    for (int i = 0; i < size_; ++i)
        data_[i] = data[i];

    calculate_fftw_plans();
}

Field::Field(const ComplexVector& x_data, const ComplexVector& y_data) {
    assert(x_data.size() == y_data.size());

    modes_ = 2;
    samples_ = x_data.size();
    time_step_ = 0;
    center_frequency_ = 0;
    circular_frequencies_ = nullptr;

    size_ = modes_ * samples_;
    data_ = reinterpret_cast<Complex*>(fftw_alloc_complex(size_));
    for (int i = 0; i < samples_; ++i)
        data_[i] = x_data[i];
    for (int i = 0; i < samples_; ++i)
        data_[i + samples_] = y_data[i];

    calculate_fftw_plans();
}

Field::Field(const Field& other) {
    modes_ = other.modes_;
    samples_ = other.samples_;
    time_step_ = other.time_step_;
    center_frequency_ = other.center_frequency_;
    circular_frequencies_ = nullptr;
    calculate_frequency_grid();

    size_ = modes_ * samples_;
    data_ = reinterpret_cast<Complex*>(fftw_alloc_complex(size_));
    for (int i = 0; i < size_; ++i)
        data_[i] = other.data_[i];

    calculate_fftw_plans();
}

Field::Field(Field&& other) {
    modes_ = other.modes_;
    samples_ = other.samples_;
    time_step_ = other.time_step_;
    center_frequency_ = other.center_frequency_;
    circular_frequencies_ = other.circular_frequencies_;

    size_ = other.size_;
    data_ = other.data_;

    forward_inplace_ = other.forward_inplace_;
    backward_inplace_ = other.backward_inplace_;

    other.circular_frequencies_ = nullptr;
    other.data_ = nullptr;
    other.forward_inplace_ = nullptr;
    other.backward_inplace_ = nullptr;
}

Field& Field::operator=(const Field& other) {
    if (this != &other) {
        free_resources();

        modes_ = other.modes_;
        samples_ = other.samples_;
        time_step_ = other.time_step_;
        center_frequency_ = other.center_frequency_;
        circular_frequencies_ = nullptr;
        calculate_frequency_grid();

        size_ = modes_ * samples_;
        data_ = reinterpret_cast<Complex*>(fftw_alloc_complex(size_));
        for (int i = 0; i < size_; ++i)
            data_[i] = other.data_[i];

        calculate_fftw_plans();
    }

    return *this;
}

Field& Field::operator=(Field&& other) {
    if (this != &other) {
        free_resources();

        modes_ = other.modes_;
        samples_ = other.samples_;
        time_step_ = other.time_step_;
        center_frequency_ = other.center_frequency_;
        circular_frequencies_ = other.circular_frequencies_;

        size_ = other.size_;
        data_ = other.data_;

        forward_inplace_ = other.forward_inplace_;
        backward_inplace_ = other.backward_inplace_;

        other.circular_frequencies_ = nullptr;
        other.data_ = nullptr;
        other.forward_inplace_ = nullptr;
        other.backward_inplace_ = nullptr;
    }

    return *this;
}

Field::~Field() { free_resources(); }

// ----------------------------------------------------------------------
// -------------------- Field accessors
// ----------------------------------------------------------------------
ComplexVector Field::vector() { return ComplexVector::proxy(size_, data_); }

Complex& Field::operator()(const int& index) {
    assert(0 <= index && index < size_);
    return data_[index];
}

Complex Field::operator()(const int& index) const {
    assert(0 <= index && index < size_);
    return data_[index];
}

Complex Field::operator()(const int& mode, const int& sample) const {
    assert(0 <= mode && mode < modes_);
    assert(0 <= sample && sample < samples_);
    return data_[sample + mode * samples_];
}

Complex& Field::operator()(const int& mode, const int& sample) {
    assert(0 <= mode && mode < modes_);
    assert(0 <= sample && sample < samples_);
    return data_[sample + mode * samples_];
}

Complex& Field::x(const int& sample) {
    assert(0 <= sample && sample < samples_);
    return data_[sample];
}

Complex& Field::y(const int& sample) {
    assert(modes_ > 1);
    assert(0 <= sample && sample < samples_);
    return data_[sample + samples_];
}

Complex Field::x(const int& sample) const {
    assert(0 <= sample && sample < samples_);
    return data_[sample];
}

Complex Field::y(const int& sample) const {
    assert(modes_ > 1);
    assert(0 <= sample && sample < samples_);
    return data_[sample + samples_];
}

ComplexVector Field::x() { return ComplexVector::proxy(samples_, data_); }

ComplexVector Field::y() {
    assert(modes_ > 1);
    return ComplexVector::proxy(samples_, data_ + samples_);
}

ComplexVector Field::x() const { return ComplexVector::proxy(samples_, data_); }

ComplexVector Field::y() const {
    assert(modes_ > 1);
    return ComplexVector::proxy(samples_, data_ + samples_);
}

ComplexVector Field::x(const int& at_begin, const int& at_end) {
    return ComplexVector::proxy(samples_ - at_begin - at_end, data_ + at_begin);
}

ComplexVector Field::y(const int& at_begin, const int& at_end) {
    assert(modes_ > 1);
    return ComplexVector::proxy(samples_ - at_begin - at_end,
                                data_ + samples_ + at_begin);
}

ComplexVector Field::x(const int& at_begin, const int& at_end) const {
    return ComplexVector::proxy(samples_ - at_begin - at_end, data_ + at_begin);
}

ComplexVector Field::y(const int& at_begin, const int& at_end) const {
    assert(modes_ > 1);
    return ComplexVector::proxy(samples_ - at_begin - at_end,
                                data_ + samples_ + at_begin);
}

// ----------------------------------------------------------------------
// -------------------- Field multiplication
// ----------------------------------------------------------------------
Field& Field::operator*=(const Complex& multiplier) {
    ComplexVector::proxy(size_, data_) *= multiplier;
    return *this;
}

Field& Field::operator*=(const ComplexVector& multipliers) {
    for (int i = 0; i < modes_; ++i)
        ComplexVector::proxy(samples_, data_ + i * samples_) *= multipliers;

    return *this;
}

Field Field::operator*(const Complex& multiplier) const {
    return Field(*this) *= multiplier;
}

Field Field::operator*(const ComplexVector& multipliers) const {
    return Field(*this) *= multipliers;
}

// ----------------------------------------------------------------------
// -------------------- Field addition
// ----------------------------------------------------------------------
Field& Field::operator+=(const Complex& summand) {
    ComplexVector::proxy(size_, data_) += summand;
    return *this;
}

Field& Field::operator+=(const ComplexVector& summands) {
    for (int i = 0; i < modes_; ++i)
        ComplexVector::proxy(samples_, data_ + i * samples_) += summands;

    return *this;
}

Field Field::operator+(const Complex& summand) const {
    return Field(*this) += summand;
}

Field Field::operator+(const ComplexVector& summands) const {
    return Field(*this) += summands;
}

// ----------------------------------------------------------------------
// -------------------- Field substraction
// ----------------------------------------------------------------------
Field& Field::operator-=(const Complex& subtrahend) {
    ComplexVector::proxy(size_, data_) -= subtrahend;
    return *this;
}

Field& Field::operator-=(const ComplexVector& subtrahends) {
    for (int i = 0; i < modes_; ++i)
        ComplexVector::proxy(samples_, data_ + i * samples_) -= subtrahends;

    return *this;
}

Field Field::operator-(const Complex& subtrahend) const {
    return Field(*this) -= subtrahend;
}

Field Field::operator-(const ComplexVector& subtrahends) const {
    return Field(*this) -= subtrahends;
}

// ----------------------------------------------------------------------
// -------------------- Field power methods
// ----------------------------------------------------------------------
double Field::power(const int& sample) const {
    assert(0 <= sample && sample < samples_);
    double power = 0;
    for (int index = sample; index < size_; index += samples_)
        power += norm(data_[index]);

    return power;
}

double Field::power(const int& mode, const int& sample) const {
    assert(0 <= mode && mode < modes_);
    assert(0 <= sample && sample < samples_);
    return norm(data_[sample + mode * samples_]);
}

double Field::peak_power() const {
    double peak_power = 0;
    for (int i = 0; i < samples_; ++i)
        if (peak_power < power(i)) peak_power = power(i);

    return peak_power;
}

double Field::average_power() const {
    double average_power = 0;
    for (int i = 0; i < samples_; ++i)
        average_power += power(i);

    return average_power / samples_;
}

Field& Field::peak_normalize(const double& power) {
    return (*this) *= std::sqrt(power / peak_power());
}

Field& Field::average_normalize(const double& power) {
    return (*this) *= std::sqrt(power / average_power());
}

// ----------------------------------------------------------------------
// -------------------- Field grid and domain methods
// ----------------------------------------------------------------------
int Field::size() const { return size_; }

int Field::modes() const { return modes_; }

int Field::samples() const { return samples_; }

double Field::get_time_step() const { return time_step_; }

double Field::get_sampling_rate() const { return 1 / time_step_; }

double Field::get_center_frequency() const { return center_frequency_; }

void Field::set_time_step(const double& time_step) {
    assert(time_step >= 0);
    time_step_ = time_step;
    calculate_frequency_grid();
}

void Field::set_sampling_rate(const double& sampling_rate) {
    set_time_step(1 / sampling_rate);
}

void Field::set_center_frequency(const double& center_frequency) {
    center_frequency_ = center_frequency;
}

double Field::dt() const { return time_step_; }

double Field::df() const { return 1 / (time_step_ * samples_); }

double Field::dw() const { return 2 * math_pi / (time_step_ * samples_); }

double Field::t(const int& sample) const {
    assert(0 <= sample && sample < samples_);
    assert(circular_frequencies_ != nullptr);
    return sample * dt();
}

double Field::f(const int& sample) const {
    assert(0 <= sample && sample < samples_);
    assert(circular_frequencies_ != nullptr);
    return circular_frequencies_[sample] / (2 * math_pi);
}

double Field::w(const int& sample) const {
    assert(0 <= sample && sample < samples_);
    assert(circular_frequencies_ != nullptr);
    return circular_frequencies_[sample];
}

RealVector Field::temporal_power() const {
    RealVector temporal_power(samples_);
    for (int i = 0; i < samples_; ++i)
        temporal_power[i] = power(i);

    return temporal_power;
}

RealVector Field::spectral_power() const {
    RealVector spectral_power(samples_);
    Field spectrum = this->fft();
    for (int i = 0; i < samples_; ++i)
        spectral_power[i] = spectrum.power(i);

    return spectral_power;
}

// ----------------------------------------------------------------------
// -------------------- Field resampling methods
// ----------------------------------------------------------------------
Field Field::upsample(const int& factor) const {
    assert(factor > 0);

    Field upsampled(modes_, factor * samples_);
    upsampled.set_time_step(time_step_ / factor);
    upsampled.set_center_frequency(center_frequency_);
    for (int i = 0; i < size_; ++i)
        upsampled.data_[factor * i] = data_[i];

    return upsampled;
}

Field Field::downsample(const int& factor) const {
    assert(factor > 0);
    assert(samples_ % factor == 0);

    Field downsampled(modes_, samples_ / factor);
    downsampled.set_time_step(time_step_ * factor);
    downsampled.set_center_frequency(center_frequency_);
    for (int i = 0; i < downsampled.size_; ++i)
        downsampled.data_[i] = data_[i * factor];

    return downsampled;
}

Field Field::decimate(const int& factor) const {
    assert(factor > 0);
    assert(samples_ % factor == 0);

    Field decimated(modes_, samples_ / factor);
    decimated.set_time_step(time_step_ * factor);
    decimated.set_center_frequency(center_frequency_);

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
    int i, j;
    int half_samples = samples_ / 2;
    Complex swap_buffer;
    for (int mode = 0; mode < modes_; ++mode) {
        for (int sample = 0; sample < half_samples; ++sample) {
            i = sample + mode * samples_;
            j = sample + mode * samples_ + half_samples;

            swap_buffer = data_[i];
            data_[j] = data_[i];
            data_[i] = swap_buffer;
        }
    }

    return *this;
}

Field& Field::fft_inplace() {
    fftw_execute(forward_inplace_);
    (*this) *= (1 / double(samples_));
    return *this;
}

Field& Field::ifft_inplace() {
    fftw_execute(backward_inplace_);
    return *this;
}

Field Field::fft() const { return Field(*this).fft_inplace(); }

Field Field::ifft() const { return Field(*this).ifft_inplace(); }

Field& Field::apply_filter(const ComplexVector& filter) {
    ComplexVector impulse_response(samples_, 0);
    int half_size = int(filter.size()) / 2;
    int shift_index = samples_ - half_size;
    for (int i = 0; i < filter.size() / 2; ++i) {
        impulse_response[i] = filter[i + half_size];
        impulse_response[shift_index + i] = filter[i];
    }

    fft_inplace() *= focss::fft(impulse_response);
    return ifft_inplace();
}

// ----------------------------------------------------------------------
// -------------------- Field private methods
// ----------------------------------------------------------------------
void Field::calculate_fftw_plans() {
    int sizes_1d[] = {samples_};
    forward_inplace_ = fftw_plan_many_dft(
        1, sizes_1d, modes_, reinterpret_cast<fftw_complex*>(data_), sizes_1d,
        1, samples_, reinterpret_cast<fftw_complex*>(data_), sizes_1d, 1,
        samples_, FFTW_FORWARD, FFTW_ESTIMATE);

    backward_inplace_ = fftw_plan_many_dft(
        1, sizes_1d, modes_, reinterpret_cast<fftw_complex*>(data_), sizes_1d,
        1, samples_, reinterpret_cast<fftw_complex*>(data_), sizes_1d, 1,
        samples_, FFTW_BACKWARD, FFTW_ESTIMATE);
}

void Field::calculate_frequency_grid() {
    if (circular_frequencies_ != nullptr) {
        delete[] circular_frequencies_;
        circular_frequencies_ = nullptr;
    }

    if (time_step_ != 0) {
        circular_frequencies_ = new double[samples_];

        double circ_freq_step = dw();
        for (int i = 0; i <= samples_ / 2; ++i)
            circular_frequencies_[i] = circ_freq_step * i;
        for (int i = samples_ / 2 + 1; i < samples_; ++i)
            circular_frequencies_[i] = circ_freq_step * (i - samples_);
    }
}

void Field::free_resources() {
    delete[] circular_frequencies_;
    fftw_free(reinterpret_cast<fftw_complex*>(data_));
    fftw_destroy_plan(forward_inplace_);
    fftw_destroy_plan(backward_inplace_);
}
}  // namespace focss
