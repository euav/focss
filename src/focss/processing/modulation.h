#ifndef FOCSS_PROCESSING_MODULATION_H_
#define FOCSS_PROCESSING_MODULATION_H_

#include "focss/field.h"
#include "focss/functions.h"

namespace focss {
namespace pulse {
ComplexVector sech(const int& n_samples, const double& width);

ComplexVector rrc(const int& n_symbols,
                  const int& oversampling,
                  const double& roll_off);

ComplexVector lorentzian(const int& n_samples,
                         const double& fwhm,
                         const double& grid_step);
}  // namespace pulse

ComplexVector random_16qam_symbols(const int& n_symbols);

Field rrc_shaping(const Field& symbols, const int& osf, const double& gain = 1);
Field rrc_sampling(const Field& signal, const int& osf, const double& gain = 1);

unsigned int demodulate_16qam(const complex_t& symbol);
Field hard_decision_16qam(const Field& signal);
Field soft_decision_16qam(const Field& signal, const double& uncertainty);
Vector<int> hd_sequence_16qam(const Field& signal);

double bit_error_rate_16qam(const Field& tx, const Field& rx);
double bit_error_rate_16qam(const Field& tx, const Field& rx, const int& cut);
double symbol_error_rate_16qam(const Field& tx, const Field& rx, const int& cut);

Field wdm_mux(const Vector<Field>& channels, const double& channel_bandwidth);

Vector<Field> wdm_demux(const Field& field,
                        const int& n_channels,
                        const double& channel_bandwidth);

Field wdm_select(const Field& field,
                 const int& channel_index,
                 const double& channel_bandwidth);

Field frequency_shift(const Field& field, const double& frequency);
Field bandpass_filter(const Field& field, const double& frequency);
}  // namespace focss

#endif  // FOCSS_PROCESSING_MODULATION_H_
