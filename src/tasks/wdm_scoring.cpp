#include "focss/equalizer/wdm_sino.h"
#include "focss/field.h"
#include "focss/functions.h"
#include "focss/processing/modulation.h"
#include "parameters.h"
#include "tasks.h"
using namespace focss;

const complex_t gray_symbols_16qam[16] = {
    complex_t(-3.0, -3.0) / sqrt(10.0), complex_t(-1.0, -3.0) / sqrt(10.0),
    complex_t(3.0, -3.0) / sqrt(10.0),  complex_t(1.0, -3.0) / sqrt(10.0),
    complex_t(-3.0, -1.0) / sqrt(10.0), complex_t(-1.0, -1.0) / sqrt(10.0),
    complex_t(3.0, -1.0) / sqrt(10.0),  complex_t(1.0, -1.0) / sqrt(10.0),
    complex_t(-3.0, 3.0) / sqrt(10.0),  complex_t(-1.0, 3.0) / sqrt(10.0),
    complex_t(3.0, 3.0) / sqrt(10.0),   complex_t(1.0, 3.0) / sqrt(10.0),
    complex_t(-3.0, 1.0) / sqrt(10.0),  complex_t(-1.0, 1.0) / sqrt(10.0),
    complex_t(3.0, 1.0) / sqrt(10.0),   complex_t(1.0, 1.0) / sqrt(10.0)};

double power_nonode(const int& power_index) {
    return 0.5 * static_cast<double>(power_index) - 6.0;
}

void pse(const arma::cx_mat& tx, arma::cx_mat& rx) {
    double angle, rms;
    double argmin_angle = 0;
    double min_rms = norm(tx - rx, "fro");
    for (int i = 1; i < 720; ++i) {
        angle = 2 * math_pi * static_cast<double>(i) / 720.0;
        rms = norm(tx - i_exp(angle) * rx, "fro");

        if (min_rms > rms) {
            min_rms = rms;
            argmin_angle = angle;
        }
    }

    rx *= i_exp(argmin_angle);
}

double q2(const arma::cx_mat& tx,
          const arma::cx_mat& rx,
          const arma::uword& cut) {
    double numerator = 0;
    for (arma::uword sample = cut; sample < tx.n_cols - cut; ++sample)
        for (arma::uword mode = 0; mode < tx.n_rows; ++mode)
            numerator += std::norm(tx(mode, sample) - rx(mode, sample));

    double denominator = 0;
    for (arma::uword sample = cut; sample < tx.n_cols - cut; ++sample)
        for (arma::uword mode = 0; mode < tx.n_rows; ++mode)
            denominator += std::norm(tx(mode, sample));

    return 10 * std::log10(denominator / numerator);
}

double ber(const arma::cx_mat& tx,
           const arma::cx_mat& rx,
           const arma::uword& cut) {
    double scale = std::sqrt(static_cast<double>(tx.n_elem)) / norm(tx, "fro");

    arma::uword errors = 0;
    for (arma::uword sample = cut; sample < tx.n_cols - cut; ++sample) {
        for (arma::uword mode = 0; mode < tx.n_rows; ++mode) {
            arma::uword tx_code = demodulate_16qam(scale * tx(mode, sample));
            arma::uword rx_code = demodulate_16qam(scale * rx(mode, sample));
            arma::uword diff = tx_code ^ rx_code;

            while (diff) {
                diff &= diff - 1;
                errors++;
            }
        }
    }

    return static_cast<double>(errors) /
           static_cast<double>(4 * tx.n_elem - 8 * cut);
}

arma::cx_mat load_field(const std::string& set,
                        const std::string& tag,
                        const int& index) {
    std::string filename = set + std::string("/") + tag + std::string("/") +
                           std::to_string(index) + std::string(".bin");
    std::fstream file(filename, std::ios::in | std::ios_base::binary);

    int n_ch;
    int n_samples;
    int n_modes;
    int flags;

    file.read((char*)&n_ch, sizeof(int));
    file.read((char*)&n_samples, sizeof(int));
    file.read((char*)&n_modes, sizeof(int));
    file.read((char*)&flags, sizeof(int));

    int n_envelopes = n_ch * n_modes;
    arma::cx_mat data(n_envelopes, n_samples, arma::fill::zeros);
    for (int sample = 0; sample < n_samples; ++sample)
        for (int env = 0; env < n_envelopes; ++env)
            file.read((char*)&data(env, sample), sizeof(complex_t));

    return data;
}

double find_certainty(const arma::cx_mat& tx, const arma::cx_mat& rx) {
    double scale = std::sqrt(static_cast<double>(tx.n_elem)) / norm(tx, "fro");
    double threshold = 1.0 / std::sqrt(5.0) / scale;

    complex_t guess;
    arma::uword tx_code, rx_code;
    for (arma::uword sample = 0; sample < tx.n_cols; ++sample) {
        for (arma::uword mode = 0; mode < tx.n_rows; ++mode) {
            tx_code = demodulate_16qam(scale * tx(mode, sample));
            rx_code = demodulate_16qam(scale * rx(mode, sample));

            if (tx_code ^ rx_code) {
                guess = gray_symbols_16qam[rx_code] / scale;
                if (threshold > std::abs(guess - rx(mode, sample))) {
                    threshold = std::abs(guess - rx(mode, sample));
                }
            }
        }
    }

    return std::sqrt(5.0) * scale * threshold;
}

arma::cx_mat sd(const arma::cx_mat& field, const double& certainty) {
    double amp =
        norm(field, "fro") / std::sqrt(static_cast<double>(field.n_elem));

    complex_t guess;
    double threshold = certainty / std::sqrt(5.0);
    arma::cx_mat result = field / amp;
    for (arma::uword sample = 0; sample < field.n_cols; ++sample) {
        for (arma::uword mode = 0; mode < field.n_rows; ++mode) {
            guess = gray_symbols_16qam[demodulate_16qam(result(mode, sample))];
            if (std::abs(guess - result(mode, sample)) < threshold)
                result(mode, sample) = guess;
        }
    }

    return amp * result;
}

arma::cx_mat hd(const arma::cx_mat& field) {
    double amp =
        norm(field, "fro") / std::sqrt(static_cast<double>(field.n_elem));

    arma::cx_mat result = field / amp;
    for (arma::uword sample = 0; sample < field.n_cols; ++sample)
        for (arma::uword mode = 0; mode < field.n_rows; ++mode)
            result(mode, sample) =
                gray_symbols_16qam[demodulate_16qam(result(mode, sample))];

    return amp * result;
}

void wdm_scoring(const int& power_index,
                 const int& span,
                 const std::string& train_set,
                 const std::string& test_set,
                 const std::string& data) {
    std::string tag = std::string("span-") + std::to_string(span);
    arma::cx_mat tx_train = load_field(train_set, "tx", power_index);
    arma::cx_mat rx_train = load_field(train_set, tag, power_index);
    arma::cx_mat tx_test = load_field(test_set, "tx", power_index);
    arma::cx_mat rx_test = load_field(test_set, tag, power_index);
    std::cout << "> data has been loaded" << std::endl;

    // learning
    WdmPPE ppe(30, 0);
    ppe.train(tx_train, rx_train);
    std::cout << "> |a| = " << std::abs(ppe.get_linear_weight()) << std::endl;

    // applying
    arma::cx_mat hd_test = hd(rx_test);
    arma::cx_mat ppe_pse = ppe.equalize(rx_test, rx_test);
    arma::cx_mat ppe_hd = ppe.equalize(hd_test, rx_test);

    // measure
    std::string name = data + std::string("/") + tag + std::string("/") +
                       std::to_string(power_index) + std::string(".csv");
    std::ofstream file(name);
    file.precision(15);

    arma::uword cut = 36;
    file << power_nonode(power_index) << ',';
    file << q2(tx_test.rows(2, 3), ppe_pse, cut) << ',';
    file << q2(tx_test.rows(2, 3), ppe_hd, cut) << ',';
    file << ber(tx_test.rows(2, 3), ppe_pse, cut) << ',';
    file << ber(tx_test.rows(2, 3), ppe_hd, cut) << ',';
    file << std::abs(ppe.get_linear_weight()) << ',';
    file << std::arg(ppe.get_linear_weight()) << std::endl;
}

void check_model(const int& power_index,
                 const int& span,
                 const std::string& train_set,
                 const std::string& test_set,
                 const std::string& data) {
    std::string tag = std::string("span-") + std::to_string(span);
    arma::cx_mat tx_train = load_field(train_set, "tx", power_index);
    arma::cx_mat rx_train = load_field(train_set, tag, power_index);
    arma::cx_mat tx_test = load_field(test_set, "tx", power_index);
    arma::cx_mat rx_test = load_field(test_set, tag, power_index);
    std::cout << "> data has been loaded" << std::endl;

    // learning
    WdmPPE ppe(35, 0);
    ppe.train(tx_train, rx_train);
    std::cout << "> |a| = " << std::abs(ppe.get_linear_weight()) << std::endl;

    // applying
    arma::cx_mat ppe_pse = ppe.equalize(rx_test, rx_test);
    arma::cx_mat ppe_rx = ppe.equalize(tx_test, rx_test);

    arma::cx_mat ppe_model = ppe.model(tx_test);

    // measure
    std::string name = data + std::string("/") + tag + std::string("/") +
                       std::to_string(power_index) + std::string(".csv");
    std::ofstream file(name);
    file.precision(15);

    arma::uword cut = 36;
    file << power_nonode(power_index) << ',';
    file << q2(tx_test.rows(2, 3), ppe_pse, cut) << ',';
    file << q2(tx_test.rows(2, 3), ppe_rx, cut) << ',';
    file << ber(tx_test.rows(2, 3), ppe_pse, cut) << ',';
    file << ber(tx_test.rows(2, 3), ppe_rx, cut) << ',';
    file << std::abs(ppe.get_linear_weight()) << ',';
    file << std::arg(ppe.get_linear_weight()) << std::endl;
}
