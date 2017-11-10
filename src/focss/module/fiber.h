#ifndef FOCSS_MODULE_FIBER_H_
#define FOCSS_MODULE_FIBER_H_

namespace focss {
struct Fiber {
    double alpha;
    double beta2;
    double gamma;

    double length;
    double wavelength;
    double kappa;

  public:
    Fiber linear() const;
    Fiber nonlinear() const;
    Fiber dispersive() const;

    Fiber operator-() const;
    Fiber operator*(const double& length_factor) const;
};
}  // namespace focss

#endif  // FOCSS_MODULE_FIBER_H_
