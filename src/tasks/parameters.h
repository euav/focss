#ifndef TASKS_PARAMETERS_H_
#define TASKS_PARAMETERS_H_

// reference units [m], [s], [W]
const double attenuation = 2e-4;                 // [dB/m]
const double dispersion = 1.7e-5;                // [s/m/m]
const double polarization_coupling = 8.0 / 9.0;  // [1]
const double nonlinearity = 1.4e-3;              // [1/W/m]
const double wavelength = 1.55e-6;               // [m]
const double fiber_length = 1e5;                 // [m]
const double bandwidth = 32e9;                   // [baud]
const double noise_figure = 4.5;                 // [dB]
const int spans = 20;                            // [1]
const double max_phase_shift = 1e-3;             // [1]
const double max_forward_step = 1e3;             // [m]
const double channel_spacing = 37.5e9;           // [Hz]
const int n_channels = 1;                        // [1]
const int osf = 16;                              // [1]
const int dbp_osf = 2;                           // [1]

#endif  // TASKS_PARAMETERS_H_
