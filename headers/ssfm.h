#ifndef SSFM_H_
#define SSFM_H_

#include "fiber.h"
#include "field.h"
#include "utility.h"

class SSFM {
    Fiber fiber;

    bool uniform_grid;
    int total_steps;
    double max_ps;
    double max_step;

    enum Mode { FORWARD, BACKWARD, CD_COMPENSATION };

  public:
    SSFM();
    SSFM(const Fiber& fiber);
    void setFiber(const Fiber& fiber);
    void setTotalSteps(const int& steps);
    void setMaximumShiftAndStep(const double& shift, const double& step);

    void run(Field& field) const;
    void run(Field& field, Mode mode) const;

  private:
    void linearStep(Field& field, const double& step) const;
    void nonlinearStep(Field& field, const double& step) const;
};

#endif  // SSFM_H_
