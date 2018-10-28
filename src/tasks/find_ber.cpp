#include "focss/equalizer/linear.h"
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/processing/modulation.h"
#include "parameters.h"
#include "systems.h"
#include "tasks.h"
using namespace focss;

void find_ber(const int& n_symbols,
              const int& power_index,
              std::ostream& data) {
    double average_power_dbm = 0.5 * static_cast<double>(power_index) - 6.0;
    std::cout << "Experiment at launch power ";
    std::cout << average_power_dbm << " dBm" << std::endl;

    // 1) modulation
    double launch_power = dbm_to_watts(average_power_dbm);
    Domain domain{2, n_symbols, 1.0 / ::bandwidth, light_speed / ::wavelength};
    Vector<Field> tx(::n_channels);
    Vector<Field> shaped(::n_channels);
    for (int ch = 0; ch < ::n_channels; ++ch) {
        tx[ch] = Field(domain, random_16qam_symbols(2 * n_symbols));
        tx[ch] *= std::sqrt(launch_power / 2);
        shaped[ch] = rrc_shaping(tx[ch], ::osf);
    }
    Field signal = wdm_mux(shaped, ::channel_spacing);
    Field tx_coi = tx[::n_channels / 2];

    // 2) propagation
    PhaseShiftEqualizer pse(720);
    for (int span = 1; span <= 1; ++span) {
        forward_propagation(signal, 1);
        Field cdc = signal;
        cd_compensation(cdc, span);
        Field signal_coi = wdm_select(cdc, 0, ::channel_spacing);
        Field rx_coi = rrc_sampling(signal_coi, ::osf);
        rx_coi = pse.train_equalize(tx_coi,rx_coi);

        data << span << ',';
        data << bit_error_rate_16qam(tx_coi, rx_coi) << std::endl;
    }
}
