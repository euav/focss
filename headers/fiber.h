#ifndef FIBER_H_
#define FIBER_H_

#include "signal.h"
#include "utility.h"
const double planck = 6.62607004081e-34;
const double light_speed = 299792458;

class Fiber {
    double alpha;
    double beta2;
    double gamma;
    double length;
    double wavelength;

    bool uniform_grid;
    int total_steps;
    double max_ps;
    double max_step;

  public:
    Fiber(const double& center_wavelength);
    void setAttenuation(const double& alpha);
    void setAttenuationDB(const double& alpha_dB);
    void setDispersionPhysical(const double& dispersion);
    void setDispersionEngineering(const double& beta2);
    void setNonlinearity(const double& gamma);
    void setFiberLength(const double& length);
    void setTotalSteps(const int& steps);
    void setMaximumShiftAndStep(const double& shift, const double& step);

    void propagate(Field& field) const;
    void compensateCD(Field& field, const double& times = 1) const;
    void amplify(Field& field) const;

  private:
    void linearStep(Field& field, const double& step) const;
    void nonlinearStep(Field& field, const double& step) const;
};

#endif  // FIBER_H_
