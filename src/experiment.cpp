#include <fstream>
#include <iostream>
#include "focss/focss.h"
#include "focss/equalizer/nonlinear_phase_recovery.h"
#include "focss/equalizer/sino.h"

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
const Fiber fiber = {db_to_natural(attenuation),
                     disp_to_beta2(dispersion, wavelength),
                     nonlinearity,
                     fiber_length,
                     wavelength,
                     pol_coupling};

void forward_propagation(Field& signal) {
    std::cout << ">> forward propagation" << std::endl;

    SSFM ssfm(fiber);
    ssfm.setMaximumShiftAndStep(max_phase_shift, max_forward_step);
    EDFA edfa(attenuation * fiber_length, noise_figure, EDFA::DECIBELS);
    for (int i = 0; i < spans; ++i) {
        ssfm.run(signal);
        edfa.amplify(signal);
        std::cout << "    * span " << i + 1 << " has done" << std::endl;
    }
    std::cout << std::endl;
}

void cd_compensation(Field& signal) {
    SSFM ssfm(-fiber.dispersive() * spans);
    ssfm.setTotalSteps(1);
    ssfm.run(signal);
}

void backward_propagation(Field& signal, const int& steps = 0) {
    SSFM ssfm(-fiber);
    EDFA edfa(attenuation * fiber_length, noise_figure, EDFA::DECIBELS);

    if (steps > 0)
        ssfm.setTotalSteps(steps);
    else
        ssfm.setMaximumShiftAndStep(max_phase_shift, max_forward_step);

    for (int i = 0; i < spans; ++i) {
        edfa.drop_power(signal);
        ssfm.run(signal);
    }
}

void final_experiment(const int& symbols_length,
                      const double& average_power_dbm,
                      std::ostream& output_stream) {
    output_stream << average_power_dbm << ',';

    // propagation
    Field symbols = random_16qam_symbols(symbols_length, bandwidth);
    Field signal = rrc_modulate(symbols, 16, dbm_to_watts(average_power_dbm));
    forward_propagation(signal);

    std::cout << ">> experimenting with dbp" << std::endl;
    for (int i = 1; i <= 11; ++i) {
        Field dbp_signal = signal.decimate(8);

        // sps parameter
        if (i <= 10)
            backward_propagation(dbp_signal, i);
        else
            backward_propagation(dbp_signal, IDEAL);

        dbp_signal = rrc_demodulate(dbp_signal, 2);

        Field train_x = symbols.chomp(0, symbols_length - training_length);
        Field train_y = dbp_signal.chomp(0, symbols_length - training_length);
        IdealPhaseRecovery pse(360);
        pse.train(train_x, train_y);
        Field valid_x = symbols.chomp(training_length, 0);
        Field valid_y = pse.equalize(dbp_signal).chomp(training_length, 0);

        output_stream << q2_factor(valid_x, valid_y) << ',';
    }

    std::cout << ">> experimenting with lms and pse" << std::endl;
    {
        Field cd_signal = signal;

        SSFM ssfm(-fiber.dispersive() * spans);
        ssfm.setTotalSteps(1);
        ssfm.run(cd_signal);

        cd_signal = rrc_demodulate(cd_signal, 16);

        Field train_x = symbols.chomp(0, symbols_length - training_length);
        Field train_y = cd_signal.chomp(0, symbols_length - training_length);
        IdealPhaseRecovery pse(360);
        LeastMeanSquare lms(3);

        pse.train(train_x, train_y);
        lms.train(train_x, train_y);

        Field valid_x = symbols.chomp(training_length, 0);
        Field valid_pse = pse.equalize(cd_signal).chomp(training_length, 0);
        Field valid_lms = lms.equalize(cd_signal).chomp(training_length, 0);

        output_stream << q2_factor(valid_x, valid_pse) << ',';
        output_stream << q2_factor(valid_x, valid_lms) << std::endl;
    }
    std::cout << ">> experiment has done" << std::endl;
}

void nlpr_experiment(const int& symbols_length,
                      const double& average_power_dbm,
                      std::ostream& output_stream) {
    output_stream << average_power_dbm << ',';

    // propagation
    Field symbols = random_16qam_symbols(symbols_length, bandwidth);
    Field signal = rrc_modulate(symbols, 16, dbm_to_watts(average_power_dbm));
    forward_propagation(signal);

    // cd compensation
    Field cd_signal = signal;
    SSFM ssfm(-fiber.dispersive() * spans);
    ssfm.setTotalSteps(1);
    ssfm.run(cd_signal);
    cd_signal = rrc_demodulate(cd_signal, 16);

    std::cout << ">> experimenting with lms, ipr and nlpr" << std::endl;
    {
        Field train_x = symbols.chomp(0, symbols_length - training_length);
        Field train_y = cd_signal.chomp(0, symbols_length - training_length);

        LeastMeanSquare lms(3);
        IdealPhaseRecovery ipr(360);
        NonlinearPhaseRecovery nlpr3(1);
        NonlinearPhaseRecovery nlpr7(3);
        NonlinearPhaseRecovery nlpr11(5);
        NonlinearPhaseRecovery nlpr15(7);
        NonlinearPhaseRecovery nlpr27(13);
        NonlinearPhaseRecovery nlpr41(20);
        SINO sino1(1);
        SINO sino2(2);
        SINO sino3(3);

        lms.train(train_x, train_y);
        ipr.train(train_x, train_y);
        nlpr3.train(train_x, train_y);
        nlpr7.train(train_x, train_y);
        nlpr11.train(train_x, train_y);
        nlpr15.train(train_x, train_y);
        nlpr27.train(train_x, train_y);
        nlpr41.train(train_x, train_y);
        sino1.train(train_x, train_y);
        sino2.train(train_x, train_y);
        sino3.train(train_x, train_y);

        Field valid_lms = lms.equalize(cd_signal).chomp(training_length, 20);
        Field valid_ipr = ipr.equalize(cd_signal).chomp(training_length, 20);
        Field valid_nlpr7 = nlpr7.equalize(cd_signal).chomp(training_length, 20);
        Field valid_nlpr11 = nlpr11.equalize(cd_signal).chomp(training_length, 20);
        Field valid_nlpr15 = nlpr15.equalize(cd_signal).chomp(training_length, 20);
        Field valid_nlpr27 = nlpr27.equalize(cd_signal).chomp(training_length, 20);
        Field valid_nlpr41 = nlpr41.equalize(cd_signal).chomp(training_length, 20);
        Field valid_sino1 = sino1.equalize(cd_signal).chomp(training_length, 20);
        Field valid_sino2 = sino2.equalize(cd_signal).chomp(training_length, 20);
        Field valid_sino3 = sino3.equalize(cd_signal).chomp(training_length, 20);
     
        Field valid_nlpr_lms;
        {
            Field train_y_nlpr = nlpr27.equalize(train_y).chomp(13, 13);
            lms.train(train_x.chomp(13, 13), train_y_nlpr);
            valid_nlpr_lms = lms.equalize(cd_signal).chomp(training_length, 20);
        }

        Field valid_lms_nlpr;
        {
            lms.train(train_x, train_y);
            Field train_y_lms = lms.equalize(train_y).chomp(2, 0);
            nlpr27.train(train_x.chomp(2, 0), train_y_lms);
            valid_lms_nlpr = nlpr27.equalize(cd_signal).chomp(training_length, 20);
        }
        
        Field valid_x = symbols.chomp(training_length, 20);

        output_stream << q2_factor(valid_x, valid_lms) << ',';
        output_stream << q2_factor(valid_x, valid_ipr) << ',';
        output_stream << q2_factor(valid_x, valid_nlpr7) << ',';
        output_stream << q2_factor(valid_x, valid_nlpr11) << ',';
        output_stream << q2_factor(valid_x, valid_nlpr15) << ',';
        output_stream << q2_factor(valid_x, valid_nlpr27) << ',';
        output_stream << q2_factor(valid_x, valid_nlpr41) << ',';
        output_stream << q2_factor(valid_x, valid_lms_nlpr) << ',';
        output_stream << q2_factor(valid_x, valid_nlpr_lms) << ',';
        output_stream << q2_factor(valid_x, valid_sino1) << ',';
        output_stream << q2_factor(valid_x, valid_sino2) << ',';
        output_stream << q2_factor(valid_x, valid_sino3) << std::endl;
    }
    std::cout << ">> experiment has done" << std::endl;
}

int main() {
    std::ofstream fout("nlpr.txt");
    fout << "#AVG(dBm),LMS(3),IPR,NLPR(7,11,15,27,41),LMS-NLPR,NLPR-LMS,SINO(1),SINO(2),SINO(3)" << std::endl;
    fout.precision(15);

    for (int i = 0; i <= 40; ++i) {
        double average_power_dbm = double(i) / 2.0 - 10.0;
        std::cout << "\nstarting experiment at power ";
        std::cout << average_power_dbm << " dBm" << std::endl;

        nlpr_experiment(2048, average_power_dbm, fout);
    }
    fout.close();

    return 0;
}
