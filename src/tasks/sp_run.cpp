#include "focss/equalizer/linear.h"
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/processing/modulation.h"
#include "parameters.h"
#include "systems.h"
#include "tasks.h"
using namespace focss;

void sp_run(const int& n_symbols,
            const double& average_power_dbm,
            std::ostream& set,
            std::ostream& data) {
    std::cout << "Experiment at launch power ";
    std::cout << average_power_dbm << " dBm" << std::endl;

    // 1) modulation
    double launch_power = dbm_to_watts(average_power_dbm);
    Domain domain{1, n_symbols, 1.0 / ::bandwidth, light_speed / ::wavelength};
    Field tx(domain, random_16qam_symbols(n_symbols));
    tx *= std::sqrt(launch_power);
    Field signal = rrc_shaping(tx, osf);
    Field btb = signal;

    // 2) propagation
    forward_propagation(signal);
    back_to_back(btb);

    // 3) save transmission for offline dsp
    Field cdc = signal;
    cd_compensation(cdc);
    Field rx = rrc_sampling(cdc, osf);
    save_transmission(set, tx, rx);

    // 4) btb scoring
    btb = rrc_sampling(btb, osf);
    data << average_power_dbm;
    data << ',' << q2_factor(tx, btb);

    // 5) dbp scoring
    std::cout << ">> dbp scoring" << std::endl;
    for (int sps = 1; sps <= 10; ++sps) {
        Field dbp = signal.decimate(osf / dbp_osf);
        backward_propagation(dbp, sps);
        dbp = rrc_sampling(dbp, dbp_osf);

        PhaseShiftEqualizer pse(3600);
        dbp = pse.train_equalize(tx, dbp);

        data << ',' << q2_factor(tx, dbp);
    }
    data << std::endl;
}
