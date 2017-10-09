#ifndef FIBER_H_
#define FIBER_H_

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

#endif  // FIBER_H_
