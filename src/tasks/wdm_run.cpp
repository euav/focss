#include <fstream>
#include <iostream>
#include "focss/equalizer/linear.h"
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/module/amplifier.h"
#include "focss/module/fiber.h"
#include "focss/processing/modulation.h"
#include "focss/solver/ssfm.h"
#include "parameters.h"
#include "systems.h"
#include "tasks.h"
using namespace focss;

const focss::Fiber fiber = {focss::db_to_natural(::attenuation),
                            focss::disp_to_beta2(::dispersion, ::wavelength),
                            ::nonlinearity,
                            ::fiber_length,
                            ::wavelength,
                            ::polarization_coupling};

double power_node(const int& power_index) {
    return 0.5 * static_cast<double>(power_index) - 6.0;
}

Field cdc(const Field& signal, const int& n_spans) {
    Field compensated = signal;
    focss::SSFM solver(::fiber.dispersive() * n_spans, 1);
    solver.run(compensated, focss::SSFM::BACKWARD);
    return compensated;
}

Field dbp(const Field& signal, const int& n_spans, const int& sps) {
    Field compensated = wdm_select(signal, 0, ::channel_spacing);
    compensated = compensated.decimate(osf / dbp_osf);

    focss::SSFM solver(::fiber, sps);
    focss::Amplifier amplifier(
        focss::db_to_linear(::attenuation * ::fiber_length), 0);

    for (int i = 0; i < n_spans; ++i) {
        amplifier.drop_power(compensated);
        solver.run(compensated, focss::SSFM::BACKWARD);
    }

    return compensated;
}

Field pse(const Field& tx, const Field& rx) {
    PhaseShiftEqualizer pse(720);
    pse.train(tx, rx);
    return pse.equalize(rx);
}

void save_field(const Vector<Field>& data,
                const std::string& set,
                const std::string& tag,
                const int& index) {
    std::string filename = set + std::string("/") + tag + std::string("/") +
                           std::to_string(index) + std::string(".bin");
    std::fstream file(filename, std::ios::out | std::ios_base::binary);

    int n_ch = data.size();
    int n_samples = data[0].samples();
    int n_modes = data[0].modes();
    int flags = 0;

    file.write((char*)&n_ch, sizeof(int));
    file.write((char*)&n_samples, sizeof(int));
    file.write((char*)&n_modes, sizeof(int));
    file.write((char*)&flags, sizeof(int));

    for (int sample = 0; sample < n_samples; ++sample)
        for (int ch = 0; ch < n_ch; ++ch)
            for (int mode = 0; mode < n_modes; ++mode)
                file.write((char*)&data[ch](mode, sample), sizeof(complex_t));
}

void wdm_run(const int& n_symbols,
             const int& power_index,
             std::string& set,
             std::string& data) {
    setup();

    // 1) modulation -- ok
    std::cout << "> modulation" << std::endl;
    double channel_power_dbm = power_node(power_index);
    double channel_power = dbm_to_watts(channel_power_dbm);
    Domain domain{2, n_symbols, 1.0 / ::bandwidth, light_speed / ::wavelength};
    Vector<Field> tx(::n_channels);
    Vector<Field> shaped(::n_channels);
    for (int ch = 0; ch < ::n_channels; ++ch) {
        tx[ch] = Field(domain, random_16qam_symbols(2 * n_symbols));
        tx[ch] *= std::sqrt(channel_power / 2);
        shaped[ch] = rrc_shaping(tx[ch], osf);
    }

    // 2) multiplexor -- ok
    std::cout << "> multiplexor" << std::endl;
    Field signal = wdm_mux(shaped, ::channel_spacing);

    // 2) propagation -- ok
    focss::SSFM solver(::fiber);
    solver.set_grid(focss::SSFM::ADAPTIVE);
    solver.set_maximum_phase_shift(::max_phase_shift);
    solver.set_maximum_step_size(::max_forward_step);
    focss::Amplifier amplifier(
        focss::db_to_linear(::attenuation * ::fiber_length),
        focss::db_to_linear(::noise_figure));

    save_field(tx, set, "tx", power_index);
    for (int span = 1; span <= ::spans; ++span) {
        std::cout << "> propagation (span " << span << ")" << std::endl;
        solver.run(signal, focss::SSFM::FORWARD);
        amplifier.amplify(signal);

        // 2.1) compensation
        std::cout << "  * compensation" << std::endl;
        Field cdc0 = cdc(signal, span);
        Vector<Field> rx = wdm_demux(cdc0, ::n_channels, ::channel_spacing);
        Field dbp1 = dbp(signal, span, 1);
        Field dbp2 = dbp(signal, span, 2);

        // 2.2) sampling -- ok
        std::cout << "  * sampling" << std::endl;
        int coi = ::n_channels / 2;
        for (int ch = 0; ch < ::n_channels; ++ch)
            rx[ch] = rrc_sampling(rx[ch], osf);
        dbp1 = rrc_sampling(dbp1, dbp_osf);
        dbp2 = rrc_sampling(dbp2, dbp_osf);

        // 2.3) phase shift equalization
        std::cout << "  * phase-shift eq" << std::endl;
        for (int ch = 0; ch < ::n_channels; ++ch)
            rx[ch] = pse(tx[ch], rx[ch]);
        dbp1 = pse(tx[coi], dbp1);
        dbp2 = pse(tx[coi], dbp2);

        // 2.4) save constellations -- ok
        std::string tag = std::string("span-") + std::to_string(span);
        save_field(rx, set, tag, power_index);

        // 2.5) save measurement -- ok
        std::cout << "  * measurement" << std::endl;
        std::string name = data + std::string("/") + tag + std::string("/") +
                           std::to_string(power_index) + std::string(".csv");
        std::ofstream file(name);
        file.precision(15);
        file << power_node(power_index);
        file << ',' << bit_error_rate_16qam(tx[coi], rx[coi]);
        file << ',' << bit_error_rate_16qam(tx[coi], dbp1);
        file << ',' << bit_error_rate_16qam(tx[coi], dbp2);
        file << std::endl;
    }
}
