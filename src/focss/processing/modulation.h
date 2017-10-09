#ifndef MODULATION_H_
#define MODULATION_H_

#include <random>
#include "focss/field.h"

const double rrc_roll_off = 0.01;
const Complex gray_symbols_16qam[16] = {Complex(-3, -3) / sqrt(10),  // 0
                                        Complex(-1, -3) / sqrt(10),  // 1
                                        Complex(3, -3) / sqrt(10),   // 2
                                        Complex(1, -3) / sqrt(10),   // 3
                                        Complex(-3, -1) / sqrt(10),  // 4
                                        Complex(-1, -1) / sqrt(10),  // 5
                                        Complex(3, -1) / sqrt(10),   // 6
                                        Complex(1, -1) / sqrt(10),   // 7
                                        Complex(-3, 3) / sqrt(10),   // 8
                                        Complex(-1, 3) / sqrt(10),   // 9
                                        Complex(3, 3) / sqrt(10),    // A
                                        Complex(1, 3) / sqrt(10),    // B
                                        Complex(-3, 1) / sqrt(10),   // C
                                        Complex(-1, 1) / sqrt(10),   // D
                                        Complex(3, 1) / sqrt(10),    // E
                                        Complex(1, 1) / sqrt(10)};   // F

Field random_16qam_symbols(const unsigned long& length,
                           const double& baudrate = 1);

Field sech_pulse(const unsigned long& nodes_quantity, const double& width);

Field rrc_filter(const double& roll_off,
                 const unsigned long& width,
                 const unsigned long& osf);

Field rrc_modulate(const Field& symbols,
                   const unsigned long& osf,
                   const double& power);
                   
Field rrc_demodulate(const Field& signal, const unsigned long& osf);

void lowpass_inplace(Field& field, const double& cutoff_frequency);

#endif  // MODULATION_H_
