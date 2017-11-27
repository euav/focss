#ifndef FOCSS_FIELD_H_
#define FOCSS_FIELD_H_

#include <fftw3.h>
#include "focss/core.h"

namespace focss {
class Field {
  public:
    Field();
    explicit Field(const int& samples);
    explicit Field(const int& modes, const int& samples);
    Field(const ComplexVector& data);
    Field(const ComplexVector& x_data, const ComplexVector& y_data);

    Field(const Field& other);
    Field(Field&& other);

    Field& operator=(const Field& other);
    Field& operator=(Field&& other);
    ~Field();

  public:
    ComplexVector vector();
    Complex& operator()(const int& index);
    Complex operator()(const int& index) const;
    Complex& operator()(const int& mode, const int& sample);
    Complex operator()(const int& mode, const int& sample) const;

    Complex& x(const int& sample);
    Complex& y(const int& sample);
    Complex x(const int& sample) const;
    Complex y(const int& sample) const;

    ComplexVector x();
    ComplexVector y();
    ComplexVector x() const;
    ComplexVector y() const;

    ComplexVector x(const int& at_begin, const int& at_end);
    ComplexVector y(const int& at_begin, const int& at_end);
    ComplexVector x(const int& at_begin, const int& at_end) const;
    ComplexVector y(const int& at_begin, const int& at_end) const;

  public:
    Field& operator*=(const Complex& multiplier);
    Field& operator*=(const ComplexVector& multipliers);
    Field operator*(const Complex& multiplier) const;
    Field operator*(const ComplexVector& multipliers) const;

    Field& operator+=(const Complex& summand);
    Field& operator+=(const ComplexVector& summands);
    Field operator+(const Complex& summand) const;
    Field operator+(const ComplexVector& summands) const;

    Field& operator-=(const Complex& subtrahend);
    Field& operator-=(const ComplexVector& subtrahends);
    Field operator-(const Complex& subtrahend) const;
    Field operator-(const ComplexVector& subtrahends) const;

  public:
    double power(const int& sample) const;
    double power(const int& mode, const int& sample) const;
    double peak_power() const;
    double average_power() const;
    Field& peak_normalize(const double& power = 1);
    Field& average_normalize(const double& power = 1);

  public:
    int size() const;
    int modes() const;
    int samples() const;
    double get_time_step() const;
    double get_sampling_rate() const;
    double get_center_frequency() const;
    void set_time_step(const double& time_step);
    void set_sampling_rate(const double& rate);
    void set_center_frequency(const double& center_frequency);

    double dt() const;
    double df() const;
    double dw() const;
    double t(const int& sample) const;
    double f(const int& sample) const;
    double w(const int& sample) const;

    RealVector temporal_power() const;
    RealVector spectral_power() const;
    RealVector time_grid() const;
    RealVector frequecny_grid() const;

  public:
    Field upsample(const int& factor) const;
    Field downsample(const int& factor) const;
    Field decimate(const int& factor) const;

  public:
    Field& fft_shift();
    Field& fft_inplace();
    Field& ifft_inplace();
    Field fft() const;
    Field ifft() const;
    Field& apply_filter(const ComplexVector& filter);

  private:
    int modes_;
    int samples_;
    double time_step_;
    double center_frequency_;
    double* circular_frequencies_;

    int size_;
    Complex* data_;

    fftw_plan forward_inplace_;
    fftw_plan backward_inplace_;

    void calculate_fftw_plans();
    void calculate_frequency_grid();
    void free_resources();
};
}  // namespace focss

#endif  // FOCSS_FIELD_H_
