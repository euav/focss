#ifndef FOCSS_FUNCTIONS_H_
#define FOCSS_FUNCTIONS_H_

#include "focss/definitions.h"
#include "focss/field.h"
#include <random>
#include <string>
#include <iostream>

namespace focss {
void setup();
void reseed_global_urng();
std::default_random_engine& global_urng();

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
complex_t i_exp(const double& x);

void fft_inplace(const int& size, complex_t* input_output);
void ifft_inplace(const int& size, complex_t* input_output);
void fft(const int& size, complex_t* input, complex_t* output);
void ifft(const int& size, complex_t* input, complex_t* output);

void fft_inplace(ComplexVector* data);
void ifft_inplace(ComplexVector* data);
ComplexVector fft(const ComplexVector& data);
ComplexVector ifft(const ComplexVector& data);

void save_transmission(std::ostream& output, const Field& tx, const Field& rx);
void load_transmission(std::istream& input, Field* tx, Field* rx);
void save_transmission(const std::string& filename, const Field& tx, const Field& rx);
void load_transmission(const std::string& filename, Field* tx, Field* rx);
}  // namespace focss

#endif  // FOCSS_FUNCTIONS_H_
