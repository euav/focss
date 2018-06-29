#ifndef FOCSS_EQUALIZER_LINEAR_H_
#define FOCSS_EQUALIZER_LINEAR_H_

/*!
 * @file linear.h
 * @brief Class @ref PhaseShiftEqualizer, @ref HomothetyEqualizer,
 * @ref ScalarEqualizer
 */

#include "focss/field.h"
#include "focss/functions.h"

namespace focss {

/*!
 * @brief rotates signal constellation by constant angle
 */
class PhaseShiftEqualizer {
  public:
    //! Constructor
    PhaseShiftEqualizer();
    explicit PhaseShiftEqualizer(const int& angle_steps);

    void train(const Field& desired, const Field& actual);
    Field equalize(const Field& original) const;
    Field train_equalize(const Field& desired, const Field& actual) const;

    double get_angle() const;
    void set_angle_steps(const int& angle_steps);

  private:
    int angle_steps_;
    double estimated_angle_;

    double estimate(const Field& desired, const Field& actual) const;
};

class HomothetyEqualizer {
  public:
    HomothetyEqualizer();

    void train(const Field& desired, const Field& actual);
    Field equalize(const Field& original) const;
    Field train_equalize(const Field& desired, const Field& actual) const;

    double get_scalar() const;

  private:
    double scalar_;
    double estimate(const Field& desired, const Field& actual) const;
};

class ScalarEqualizer {
  public:
    ScalarEqualizer();

    void train(const Field& desired, const Field& actual);
    Field equalize(const Field& original) const;
    Field train_equalize(const Field& desired, const Field& actual) const;

    complex_t get_scalar() const;

  private:
    complex_t scalar_;
    complex_t estimate(const Field& desired, const Field& actual) const;
};
}  // namespace focss

#endif  // FOCSS_EQUALIZER_LINEAR_H_
