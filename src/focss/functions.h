#ifndef FOCSS_FUNCTIONS_H_
#define FOCSS_FUNCTIONS_H_

#include "focss/core.h"
#include "focss/field.h"
#include <random>

namespace focss {
std::default_random_engine& global_urng();
void reseed_global_urng();

double evm2_factor(const Field& tx, const Field& rx);
double evm2_factor(const Field& tx, const Field& rx, const int& cut);
double evm2_factor(const ComplexVector& tx, const ComplexVector& rx);

double q2_factor(const Field& tx, const Field& rx);
double q2_factor(const Field& tx, const Field& rx, const int& cut);
double q2_factor(const ComplexVector& tx, const ComplexVector& rx);

double db_to_linear(const double& db_value);
double dbm_to_watts(const double& dbm_power);
double db_to_natural(const double& db_value);
double disp_to_beta2(const double& dispersion, const double& wavelength);

double sinc(const double& x);
Complex i_exp(const double& x);

void fft_inplace(const int& size, Complex* input_output);
void ifft_inplace(const int& size, Complex* input_output);
void fft(const int& size, Complex* input, Complex* output);
void ifft(const int& size, Complex* input, Complex* output);

void fft_inplace(ComplexVector* data);
void ifft_inplace(ComplexVector* data);
ComplexVector fft(const ComplexVector& data);
ComplexVector ifft(const ComplexVector& data);

void save_transmission(const char* filename, const Field& tx, const Field& rx);
void load_transmission(const char* filename, Field* tx, Field* rx);
}  // namespace focss

#endif  // FOCSS_FUNCTIONS_H_
