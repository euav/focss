#include "focss/equalizer/linear.h"
#include "focss/equalizer/sino.h"
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/processing/modulation.h"
#include "parameters.h"
#include "tasks.h"
using namespace focss;

void sp_scoring(const int& n_symbols,
                const double& average_power_dbm,
                std::istream& train,
                std::istream& test,
                std::ostream& data) {
    std::cout << "Experiment at launch power ";
    std::cout << average_power_dbm << " dBm" << std::endl;

    Domain dtrain{1, n_symbols, 1.0 / ::bandwidth, light_speed / ::wavelength};
    Domain dtest{1, 65536, 1.0 / ::bandwidth, light_speed / ::wavelength};
    Field tx_train(dtrain);
    Field rx_train(dtrain);
    Field tx_test(dtest);
    Field rx_test(dtest);
    load_transmission(train, &tx_train, &rx_train);
    load_transmission(test, &tx_test, &rx_test);

    ScalarEqualizer seq;
    PhaseShiftEqualizer pse(3600);
    SINO sino(62, 0, SinoType::FIXED);
    sino.set_linear_weight(1.0);
    // 2.456176373e-2 -- weight for l2
    // 8e-11

    pse.train(tx_train, rx_train);
    sino.train(tx_train, rx_train);

    Field pse_test = pse.equalize(rx_test);
    Field hd_test = hard_decision_16qam(pse_test);

    Field rx_ppe_test = sino.equalize(rx_test, rx_test);
    rx_ppe_test = seq.train_equalize(tx_test, rx_ppe_test);

    Field pse_ppe_test = sino.equalize(rx_test, pse_test);
    pse_ppe_test = seq.train_equalize(tx_test, pse_ppe_test);

    Field hd_ppe_test = sino.equalize(rx_test, hd_test);
    hd_ppe_test = seq.train_equalize(tx_test, hd_ppe_test);

    Field tx_ppe_test = sino.equalize(rx_test, tx_test);
    tx_ppe_test = seq.train_equalize(tx_test, tx_ppe_test);

    int cut = 100;
    data << average_power_dbm << ',';
    data << q2_factor(tx_test, pse_test, cut) << ',';
    data << q2_factor(tx_test, rx_ppe_test, cut) << ',';
    data << q2_factor(tx_test, pse_ppe_test, cut) << ',';
    data << q2_factor(tx_test, hd_ppe_test, cut) << ',';
    data << q2_factor(tx_test, tx_ppe_test, cut) << ',';
    data << bit_error_rate_16qam(tx_test, pse_test, cut) << ',';
    data << bit_error_rate_16qam(tx_test, rx_ppe_test, cut) << ',';
    data << bit_error_rate_16qam(tx_test, pse_ppe_test, cut) << ',';
    data << bit_error_rate_16qam(tx_test, hd_ppe_test, cut) << ',';
    data << bit_error_rate_16qam(tx_test, tx_ppe_test, cut) << std::endl;
}
