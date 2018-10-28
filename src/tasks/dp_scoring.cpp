#include "focss/equalizer/linear.h"
#include "focss/equalizer/sino.h"
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/processing/modulation.h"
#include "parameters.h"
#include "tasks.h"
using namespace focss;

void dp_scoring(const int& n_symbols,
                const double& average_power_dbm,
                std::istream& train,
                std::istream& test,
                std::ostream& data) {
    std::cout << "Experiment at launch power ";
    std::cout << average_power_dbm << " dBm" << std::endl;

    Domain domain{2, n_symbols, 1.0 / ::bandwidth, light_speed / ::wavelength};
    Field tx_train(domain);
    Field rx_train(domain);
    Field tx_test(domain);
    Field rx_test(domain);
    load_transmission(train, &tx_train, &rx_train);
    load_transmission(test, &tx_test, &rx_test);
    load_transmission(train, &tx_train, &rx_train);
    load_transmission(test, &tx_test, &rx_test);

    ScalarEqualizer seq;
    PhaseShiftEqualizer pse(3600);
    SINO sino(36 , 1e-4, SinoType::ORDINARY);
    // sino.set_linear_weight(1.0);

    pse.train(tx_train, rx_train);
    sino.train(tx_train, rx_train);
    std::cout << "> |a| = " << std::abs(sino.get_linear_weight()) << std::endl;

    Field pse_test = pse.equalize(rx_test);
    Field pse_ppe_test = sino.equalize(rx_test, pse_test);
    pse_ppe_test = seq.train_equalize(tx_test, pse_ppe_test);

    int cut = 100;
    data << average_power_dbm << ',';
    data << q2_factor(tx_test, pse_test, cut) << ',';
    data << q2_factor(tx_test, pse_ppe_test, cut) << ',';
    data << bit_error_rate_16qam(tx_test, pse_test, cut) << ',';
    data << bit_error_rate_16qam(tx_test, pse_ppe_test, cut) << std::endl;
}
