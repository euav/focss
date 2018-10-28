#ifndef FOCSS_DEFINITIONS_H_
#define FOCSS_DEFINITIONS_H_

#include <complex>
#include "focss/vector.h"

namespace focss {
typedef double real_t;
typedef std::complex<real_t> complex_t;
typedef Vector<real_t> RealVector;
typedef Vector<complex_t> ComplexVector;

const real_t math_pi = 3.14159265358979323846264338;  // [1]
const real_t planck = 6.62607004081818181818181e-34;  // [m^2 kg / s]
const real_t light_speed = 299792458;                 // [m / s]
const complex_t i_unit = complex_t(0, 1);
}  // namespace focss

#endif  // FOCSS_DEFINITIONS_H_
