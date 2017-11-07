#ifndef FOCSS_UTILITY_H_
#define FOCSS_UTILITY_H_

#include <cmath>
#include <complex>
#include <random>
#include "focss/field.h"

void relax_max(double& a, const double& b);
void relax_mix(double& a, const double& b);
double sinc(const double& x);

double evm2_factor(const Field& tx, const Field& rx);
double q2_factor(const Field& tx, const Field& rx);

double db_to_linear(const double& db_value);
double dbm_to_watts(const double& dbm_power);
double db_to_natural(const double& db_value);
double disp_to_beta2(const double& dispersion, const double& wavelength);

Complex i_exp(const double& x);

#endif  // FOCSS_UTILITY_H_
