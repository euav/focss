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
    double fiber_length;
    double wavelength;
    double sampling_rate;
    int total_steps;

  public:
    Fiber(const double& center_wavelength);
    void setAttenuation(const double& alpha);
    void setAttenuationDB(const double& alpha_dB);
    void setDispersionPhysical(const double& dispersion);
    void setDispersionEngineering(const double& beta2);
    void setNonlinearity(const double& gamma);
    void setFiberLength(const double& length);
    void setTotalSteps(const double& steps);
    void setSamplingRate(const double& rate);

    void propagate(Field& field) const;
    void compensateCD(Field& field, const double& times = 1) const;
    void amplify(Field& field) const;

  private:
    Field estimateLinearity(const int& samples) const;
};

#endif  // FIBER_H_
