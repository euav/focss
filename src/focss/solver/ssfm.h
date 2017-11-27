#ifndef FOCSS_SOLVER_SSFM_H_
#define FOCSS_SOLVER_SSFM_H_

#include "focss/field.h"
#include "focss/module/fiber.h"

namespace focss {
class SSFM {
    Fiber fiber_;

    bool uniform_grid_;
    int total_steps_;
    double maximum_phase_shift_;
    double maximum_step_;

  public:
    typedef const double& Direction;
    static const double FORWARD;
    static const double BACKWARD;

    SSFM();
    SSFM(const Fiber& fiber);
    SSFM(const Fiber& fiber, const int& total_steps);
    void set_fiber(const Fiber& fiber);
    void set_total_steps(const int& total_steps);
    void set_maximum_shift_and_step(const double& maximum_phase_shift,
                                    const double& maximum_step);

    void run(Field& field, Direction direction = FORWARD) const;

  private:
    void uniform_solve(Field& field, Direction direction) const;
    void adaptive_solve(Field& field, Direction direction) const;
    double estimate_next_step(const Field& field) const;

    void linear_step(Field& field, const double& step) const;
    void nonlinear_step(Field& field, const double& step) const;
    void linear_step_without_fft(Field& field, const double& step) const;
    void nonlinear_step_with_fft(Field& field, const double& step) const;
};
}  // namespace focss

#endif  // FOCSS_SOLVER_SSFM_H_
