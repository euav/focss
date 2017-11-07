#ifndef FOCSS_FUNCTIONS_H_
#define FOCSS_FUNCTIONS_H_

#include <cmath>
#include <complex>
#include "focss/vector.h"

namespace focss {
typedef std::complex<double> Complex;
typedef Vector<Complex> ComplexVector;
typedef Vector<double> RealVector;

double sinc(const double& x);

double evm2_factor(const ComplexVector& tx, const ComplexVector& rx);
double q2_factor(const ComplexVector& tx, const ComplexVector& rx);

double db_to_linear(const double& db_value);
double dbm_to_watts(const double& dbm_power);
double db_to_natural(const double& db_value);
double disp_to_beta2(const double& dispersion, const double& wavelength);

Complex i_exp(const double& x);

void fft_inplace(ComplexVector* data);
void ifft_inplace(ComplexVector* data);
ComplexVector fft(const ComplexVector& data);
ComplexVector ifft(const ComplexVector& data);
}  // namespace focss

#endif  // FOCSS_FUNCTIONS_H_
