#ifndef MODULATION_H_
#define MODULATION_H_

#include <vector>
#include "utility.h"

typedef int InformationType;
typedef std::vector<InformationType> Information;

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

const Information constellation_16qam = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

Field modulate_16qam(const Information& information);
Information demodulate_16qam(const Field& signal);

Field sech_pulse(const int& nodes_quantity, const double& width);
Field rrc_filter(const int& samples, const double& roll_off, const int& N);
double bit_error_rate(const Information& data_tx, const Information& data_rx);

#endif  // MODULATION_H_
