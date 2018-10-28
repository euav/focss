#include "focss/equalizer/linear.h"
#include "focss/equalizer/sino.h"
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/processing/modulation.h"
#include "parameters.h"
#include "tasks.h"
using namespace focss;

void redyuk_scoring() {
    Domain d_train{2, 106496, 1.0 / ::bandwidth, light_speed / ::wavelength};
    Field tx_train(d_train), rx_train(d_train);
    load_transmission("data/redyuk/train/full.csv", &tx_train, &rx_train);

    std::ofstream head("data/redyuk/data/dp-full-head.csv");
    std::ofstream tail("data/redyuk/data/dp-full-tail.csv");
    std::ofstream full("data/redyuk/data/dp-full-full.csv");

    ScalarEqualizer seq;
    for (int complexity = 0; complexity <= 50; ++complexity) {
        SINO sino(complexity, 1e-4, SinoType::ORDINARY);
        sino.train(tx_train, rx_train);

        { // head test
            Domain d{2, 53248, 1.0 / ::bandwidth, light_speed / ::wavelength};
            Field tx(d), rx(d);
            load_transmission("data/redyuk/test/head.csv", &tx, &rx);

            Field ppe_pse = sino.equalize(rx, rx);
            ppe_pse = seq.train_equalize(tx, ppe_pse);

            Field ppe_hd = sino.equalize(rx, hard_decision_16qam(rx));
            ppe_hd = seq.train_equalize(tx, ppe_hd);

            int cut = 100;
            head << complexity << ',';
            head << bit_error_rate_16qam(tx, ppe_pse, cut) << ',';
            head << bit_error_rate_16qam(tx, ppe_hd, cut) << std::endl;
        }

        { // tail test
            Domain d{2, 53248, 1.0 / ::bandwidth, light_speed / ::wavelength};
            Field tx(d), rx(d);
            load_transmission("data/redyuk/test/tail.csv", &tx, &rx);

            Field ppe_pse = sino.equalize(rx, rx);
            ppe_pse = seq.train_equalize(tx, ppe_pse);

            Field ppe_hd = sino.equalize(rx, hard_decision_16qam(rx));
            ppe_hd = seq.train_equalize(tx, ppe_hd);

            int cut = 100;
            tail << complexity << ',';
            tail << bit_error_rate_16qam(tx, ppe_pse, cut) << ',';
            tail << bit_error_rate_16qam(tx, ppe_hd, cut) << std::endl;
        }

        { // full test
            Domain d{2, 106496, 1.0 / ::bandwidth, light_speed / ::wavelength};
            Field tx(d), rx(d);
            load_transmission("data/redyuk/test/full.csv", &tx, &rx);

            Field ppe_pse = sino.equalize(rx, rx);
            ppe_pse = seq.train_equalize(tx, ppe_pse);

            Field ppe_hd = sino.equalize(rx, hard_decision_16qam(rx));
            ppe_hd = seq.train_equalize(tx, ppe_hd);

            int cut = 100;
            full << complexity << ',';
            full << bit_error_rate_16qam(tx, ppe_pse, cut) << ',';
            full << bit_error_rate_16qam(tx, ppe_hd, cut) << std::endl;
        }
    }
}
