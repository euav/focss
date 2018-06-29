#ifndef FOCSS_SOLVER_SSFM_H_
#define FOCSS_SOLVER_SSFM_H_

#include "focss/field.h"
#include "focss/module/fiber.h"

namespace focss {
class SSFM {
  public:
    typedef const int Grid;
    static const int UNIFORM;
    static const int LOGARITHMIC;
    static const int ADAPTIVE;

    typedef const double Direction;
    static const double FORWARD;
    static const double BACKWARD;

  public:
    SSFM();
    SSFM(const Fiber& fiber);
    SSFM(const Fiber& fiber, const int& total_steps);

  public:
    void set_fiber(const Fiber& fiber);
    void set_grid(const Grid& grid);
    void set_total_steps(const int& total_steps);
    void set_maximum_phase_shift(const double& maximum_phase_shift);
    void set_maximum_step_size(const double& maximum_step_size);

    void run(Field& field, const Direction& direction = FORWARD);

  private:
    void uniform_solve(Field& field, const Direction& direction) const;
    void logarithmic_solve(Field& field, const Direction& direction) const;
    void adaptive_solve(Field& field, const Direction& direction) const;

    int estimate_logarithmic_steps(const Field& field) const;
    double estimate_next_adaptive_step(const Field& field) const;

    inline void linear_step(Field& field, const double& step) const;
    inline void nonlinear_step(Field& field, const double& step) const;
    inline void nofft_linear_step(Field& field, const double& step) const;
    inline void fft_nonlinear_step(Field& field, const double& step) const;

  private:
    Fiber fiber_;

    int grid_;
    int total_steps_;
    double maximum_phase_shift_;
    double maximum_step_size_;
    double* disp_cache_;
};
}  // namespace focss

#endif  // FOCSS_SOLVER_SSFM_H_
