#ifndef UTILITY_H_
#define UTILITY_H_

#include <cmath>
#include <complex>
#include <random>
#include "signal.h"

void relax_max(double& a, const double& b);
void relax_mix(double& a, const double& b);
double sinc(const double& x);

double q2_factor(const Signal& tx, const Signal& rx);
double dbm_to_watts(const double& power_dbm);
void awgn_generator(Signal& original, const double& variance);

Complex i_exp(const double& x);

#endif  // UTILITY_H_
