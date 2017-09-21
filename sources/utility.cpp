#include "../headers/utility.h"

void relax_max(double& a, const double& b) {
    if (a < b) a = b;
}

void relax_mix(double& a, const double& b) {
    if (a > b) a = b;
}

double sinc(const double& x) {
    if (x == 0) return 1;
    return std::sin(M_PI * x) / M_PI / x;
}

/*
void save_complex(const Complex &number, std::ostream &os) {
    os << number.real();
    if (number.imag() > 0) os << "+";
    os << number.imag();
    os << "j";
}

void save_transmission(const Field &tx, const Field &rx, std::ostream &os) {
    int transmission_size = std::min(tx.size(), rx.size());
    if (transmission_size > 0) {
        os << "id,";
        os << "transmitted_real,transmitted_imag,";
        os << "received_real,received_imag\n";
        os << "1,";
        save_complex(tx[0], os);
        os << ",";
        save_complex(rx[0], os);
    }

    for (int i = 1; i < transmission_size; ++i) {
        os << "\n" << i + 1 << ",";
        save_complex(tx[i], os);
        os << ",";
        save_complex(rx[i], os);
    }
}
*/
double evm2_factor(const Field& tx, const Field& rx) {
    double numerator = 0;
    double denominator = 0;

    int transmission_size = std::min(tx.size(), rx.size());
    for (int i = 0; i < transmission_size; ++i) {
        numerator += norm(tx[i] - rx[i]);
        denominator += norm(tx[i]);
    }

    return numerator / denominator;
}

double q2_factor(const Field& tx, const Field& rx) {
    return 10 * std::log10(1 / evm2_factor(tx, rx));
}

double dbm_to_watts(const double& power_dbm) {
    return 1e-3 * std::pow(10, power_dbm / 10);
}

void awgn_generator(Field& original, const double& variance) {
    std::default_random_engine generator(time(0));
    std::normal_distribution<double> noise(0, std::sqrt(variance / 2));
    for (int i = 0; i < original.size(); ++i)
        original[i] += Complex(noise(generator), noise(generator));
}

Complex i_exp(const double& x) {
    return Complex(std::cos(x), std::sin(x));
}
