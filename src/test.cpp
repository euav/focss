#include <fstream>
#include <iostream>
#include "focss/utility.h"

// reference units [m], [s], [W]
const double attenuation = 2e-4;        // [dB/m]
const double dispersion = 1.7e-5;       // [s/m/m]
const double nonlinearity = 1.4e-3;     // [1/W/m]
const double wavelength = 1.55e-6;      // [m]
const double fiber_length = 1e5;        // [m]
const double pol_coupling = 8.0 / 9.0;  // [1]
const double bandwidth = 32e9;          // [baud]
const double noise_figure = 4.5;        // [dB]
const double spans = 10;                // [1]
const double max_phase_shift = 1e-2;    // [1]
const double max_forward_step = 1e3;    // [m]
const int training_length = 1536;        // [1]
const int IDEAL = 0;



int main() {
    std::cout << disp_to_beta2(dispersion, wavelength) << std::endl;

    return 0;
}
