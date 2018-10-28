 /*!
 * @file
 * @brief Class @ref Domain, @ref Field
 */
#ifndef FOCSS_FIELD_H_
#define FOCSS_FIELD_H_

#include <fftw3.h>
#include "focss/definitions.h"

namespace focss {

/*!
 *  @brief Domain structure represents size and frequency of Field.
 *
 *  Every Field instance constructed by Domain instance which responsible
 *  for maintaing of number of samples, modes, time step between samples
 *  and center frequency of field.
 */
struct Domain {
    int modes; //!< number of modes, polarizations also considered as modes
    int samples; //!< number of samples per mode
    real_t time_step; //!< time step between consequent samples
    real_t center_frequency; //!< center frequency of field

    /*!
     * @brief total number of samples in the field within all modes
     * @return number of samples multiplied by number of modes
     */
    inline int size() const { return modes * samples; }

    /*!
     * @brief checks if two domains are identical
     * @param other domain to compare with
     * @return logical value if two domains are the same
     */
    inline bool operator==(const Domain& other) const {
        return modes == other.modes && samples == other.samples &&
               time_step == other.time_step &&
               center_frequency == other.center_frequency;
    }
};

class Field {
  public:
    Field();
    explicit Field(const Domain& domain);
    explicit Field(const Domain& domain, const ComplexVector& data);

    Field(const Field& other);
    Field(Field&& other);

    Field& operator=(const Field& other);
    Field& operator=(Field&& other);
    ~Field();

  public:
    complex_t operator[](const int& index) const;
    complex_t operator()(const int& mode, const int& sample) const;

    complex_t& operator[](const int& index);
    complex_t& operator()(const int& mode, const int& sample);

    ComplexVector operator()(const int& mode);

  public:
    real_t power(const int& sample) const;
    real_t peak_power() const;
    real_t average_power() const;
    Field& peak_normalize(const real_t& power = 1);
    Field& average_normalize(const real_t& power = 1);

  public:
    int size() const;
    int modes() const;
    int samples() const;

    real_t duration() const;
    real_t bandwidth() const;
    real_t time_step() const;
    real_t sampling_rate() const;
    real_t center_frequency() const;

    real_t dt() const;
    real_t df() const;
    real_t dw() const;
    real_t t(const int& sample) const;
    real_t f(const int& sample) const;
    real_t w(const int& sample) const;

    Domain domain() const;

    RealVector time_grid() const;
    RealVector frequency_grid() const;
    RealVector temporal_power() const;
    RealVector spectral_power() const;

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
    Field& filter_inplace(const ComplexVector& impulse_response);
    Field filter(const ComplexVector& impulse_response) const;

  public:
    Field& operator*=(const real_t&);
    Field& operator*=(const complex_t&);
    Field& operator*=(const RealVector&);
    Field& operator*=(const ComplexVector&);

    Field& operator/=(const real_t&);
    Field& operator/=(const complex_t&);
    Field& operator/=(const RealVector&);
    Field& operator/=(const ComplexVector&);

    Field& operator+=(const Field&);
    Field& operator-=(const Field&);

    friend Field operator*(const Field&, const real_t&);
    friend Field operator*(const Field&, const complex_t&);
    friend Field operator*(const Field&, const RealVector&);
    friend Field operator*(const Field&, const ComplexVector&);
    friend Field operator*(const real_t&, const Field&);
    friend Field operator*(const complex_t&, const Field&);
    friend Field operator*(const RealVector&, const Field&);
    friend Field operator*(const ComplexVector&, const Field&);

    friend Field operator/(const Field&, const real_t&);
    friend Field operator/(const Field&, const complex_t&);
    friend Field operator/(const Field&, const RealVector&);
    friend Field operator/(const Field&, const ComplexVector&);
    friend Field operator/(const real_t&, const Field&);
    friend Field operator/(const complex_t&, const Field&);
    friend Field operator/(const RealVector&, const Field&);
    friend Field operator/(const ComplexVector&, const Field&);

    friend Field operator+(const Field&, const Field&);
    friend Field operator-(const Field&, const Field&);

  private:
    void prepare_fft_plans();
    void release_resources();

  private:
    int modes_;
    int samples_;
    real_t time_step_;
    real_t center_frequency_;

    int size_;
    complex_t* data_;

    fftw_plan forward_inplace_;
    fftw_plan backward_inplace_;
};
}  // namespace focss

#endif  // FOCSS_FIELD_H_
