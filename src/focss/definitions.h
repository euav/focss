#ifndef FOCSS_CORE_H_
#define FOCSS_CORE_H_

#include <complex>
#include "focss/vector.h"

namespace focss {
const double math_pi = 3.14159265358979323846264338;  // [1]
const double planck = 6.62607004081818181818181e-34;  // [m^2 kg / s]
const double light_speed = 299792458;                 // [m/s]

typedef std::complex<double> Complex;
typedef Vector<Complex> ComplexVector;
typedef Vector<double> RealVector;
const Complex i_unit = Complex(0, 1);
}  // namespace focss

#endif  // FOCSS_CORE_H_
