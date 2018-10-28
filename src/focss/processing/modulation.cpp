#include "modulation.h"
#include <algorithm>
#include <cassert>
#include <random>
#include "focss/field.h"
#include "focss/functions.h"

namespace focss {
const complex_t gray_symbols_16qam[16] = {
    complex_t(-3.0, -3.0) / sqrt(10.0), complex_t(-1.0, -3.0) / sqrt(10.0),
    complex_t(3.0, -3.0) / sqrt(10.0),  complex_t(1.0, -3.0) / sqrt(10.0),
    complex_t(-3.0, -1.0) / sqrt(10.0), complex_t(-1.0, -1.0) / sqrt(10.0),
    complex_t(3.0, -1.0) / sqrt(10.0),  complex_t(1.0, -1.0) / sqrt(10.0),
    complex_t(-3.0, 3.0) / sqrt(10.0),  complex_t(-1.0, 3.0) / sqrt(10.0),
    complex_t(3.0, 3.0) / sqrt(10.0),   complex_t(1.0, 3.0) / sqrt(10.0),
    complex_t(-3.0, 1.0) / sqrt(10.0),  complex_t(-1.0, 1.0) / sqrt(10.0),
    complex_t(3.0, 1.0) / sqrt(10.0),   complex_t(1.0, 1.0) / sqrt(10.0)};

ComplexVector pulse::sech(const int& n_samples, const double& width) {
    ComplexVector pulse(n_samples);

    double argument;
    for (int i = 0; i < n_samples; ++i) {
        argument = double(i) / (n_samples - 1.0);
        argument = 20 * argument - 20 / 2.0;
        pulse[i] = 1.0 / std::cosh(argument / width);
        pulse[i] /= width;
    }

    return pulse;
}

ComplexVector pulse::lorentzian(const int& n_samples,
                                const double& fwhm,
                                const double& grid_step) {
    ComplexVector pulse(n_samples);

    double argument;
    for (int i = 0; i < n_samples; ++i) {
        argument = grid_step * double(i - n_samples / 2);
        pulse[i] = 1 / (1 + 4 * argument * argument / fwhm / fwhm);
    }

    return pulse;
}

ComplexVector pulse::rrc(const int& n_symbols,
                         const int& oversampling,
                         const double& roll_off) {
    ComplexVector filter(n_symbols * oversampling);

    for (int i = 0; i < filter.size(); ++i) {
        double t = double(i) / oversampling - n_symbols / 2.0;

        if (t == 0) {
            filter[i] = 1 - roll_off + 4 * roll_off / math_pi;
        } else if (std::abs(t) == 1 / roll_off / 4) {
            filter[i] += (1 + 2 / math_pi) * std::sin(math_pi / roll_off / 4);
            filter[i] += (1 - 2 / math_pi) * std::cos(math_pi / roll_off / 4);
            filter[i] *= roll_off / std::sqrt(2.0);
        } else {
            filter[i] += std::cos(math_pi * t * (1 + roll_off));
            filter[i] *= 4 * roll_off * t;
            filter[i] += std::sin(math_pi * t * (1 - roll_off));
            filter[i] /= 1 - 16 * roll_off * roll_off * t * t;
            filter[i] /= math_pi * t;
        }
    }

    double power = 0;
    for (int i = 0; i < filter.size(); ++i)
        power += norm(filter[i]);
    return filter *= filter.size() / power;
}

ComplexVector random_16qam_symbols(const int& n_symbols) {
    ComplexVector symbols(n_symbols);
    for (int i = 0; i < n_symbols; ++i)
        symbols[i] = gray_symbols_16qam[i % 16];

    std::shuffle(symbols.raw(), symbols.raw() + n_symbols, global_urng());
    return symbols;
}

const double rrc_roll_off = 0.01;
Field rrc_shaping(const Field& symbols, const int& osf, const double& gain) {
    auto rrc = pulse::rrc(symbols.samples(), osf, rrc_roll_off);
    auto signal = symbols.upsample(osf).filter_inplace(rrc);
    return (std::sqrt(gain) * osf) * signal;
}

Field rrc_sampling(const Field& signal, const int& osf, const double& gain) {
    auto rrc = pulse::rrc(signal.samples() / osf, osf, rrc_roll_off);
    auto symbols = signal.filter(rrc).downsample(osf);
    return std::sqrt(gain) * symbols;
}

unsigned int demodulate_16qam(const complex_t& symbol) {
    unsigned int argmin = 0;
    double min = norm(symbol - gray_symbols_16qam[0]);

    for (unsigned int code = 1; code < 16; ++code) {
        if (min > norm(symbol - gray_symbols_16qam[code])) {
            min = norm(symbol - gray_symbols_16qam[code]);
            argmin = code;
        }
    }

    return argmin;
}

Field hard_decision_16qam(const Field& signal) {
    double power = signal.average_power() / signal.modes();

    auto hd = signal / std::sqrt(power);
    for (int i = 0; i < signal.size(); ++i)
        hd[i] = gray_symbols_16qam[demodulate_16qam(hd[i])];

    return std::sqrt(power) * hd;
}

Field soft_decision_16qam(const Field& signal, const double& uncertainty) {
    double power = signal.average_power() / signal.modes();

    complex_t guess;
    double pivot = 2.0 / sqrt(10);
    double threshold = sqrt(10) * uncertainty;
    auto sd = signal / std::sqrt(power);
    for (int i = 0; i < signal.size(); ++i) {
        guess = gray_symbols_16qam[demodulate_16qam(sd[i])];

        if (std::abs(sd[i].real()) < threshold ||
            std::abs(sd[i].real() - pivot) < threshold ||
            std::abs(sd[i].real() + pivot) < threshold)
            guess.real(sd[i].real());

        if (std::abs(sd[i].imag()) < threshold ||
            std::abs(sd[i].imag() - pivot) < threshold ||
            std::abs(sd[i].imag() + pivot) < threshold)
            guess.imag(sd[i].imag());

        sd[i] = guess;
    }

    return std::sqrt(power) * sd;
}

Vector<int> hd_sequence_16qam(const Field& signal) {
    double scale = std::sqrt(signal.modes() / signal.average_power());

    Vector<int> sequence(signal.size());
    for (int i = 0; i < signal.size(); ++i)
        sequence[i] = demodulate_16qam(scale * signal[i]);

    return sequence;
}

unsigned int bit_dist(const unsigned int& a, const unsigned int& b) {
    // Brian Kernighan's algorithm
    unsigned int diff = a ^ b;
    unsigned int distance = 0;
    while (diff) {
        diff &= diff - 1;
        distance++;
    }

    return distance;
}

double bit_error_rate_16qam(const Field& tx, const Field& rx) {
    assert(tx.size() == rx.size());

    double scale = std::sqrt(tx.modes() / tx.average_power());

    unsigned int errors = 0;
    for (int i = 0; i < tx.size(); ++i) {
        unsigned int tx_code = demodulate_16qam(scale * tx[i]);
        unsigned int rx_code = demodulate_16qam(scale * rx[i]);
        errors += bit_dist(tx_code, rx_code);
    }

    return static_cast<double>(errors) / static_cast<double>(4 * tx.size());
}

double bit_error_rate_16qam(const Field& tx, const Field& rx, const int& cut) {
    assert(tx.modes() == rx.modes());
    assert(tx.samples() == rx.samples());

    double scale = std::sqrt(tx.modes() / tx.average_power());

    unsigned int errors = 0;
    for (int mode = 0; mode < tx.modes(); ++mode) {
        for (int sample = cut; sample < tx.samples() - cut; ++sample) {
            unsigned int tx_code = demodulate_16qam(scale * tx(mode, sample));
            unsigned int rx_code = demodulate_16qam(scale * rx(mode, sample));
            errors += bit_dist(tx_code, rx_code);
        }
    }

    return static_cast<double>(errors) /
           static_cast<double>(4 * tx.size() - 8 * cut);
}

double symbol_error_rate_16qam(const Field& tx,
                               const Field& rx,
                               const int& cut) {
    assert(tx.modes() == rx.modes());
    assert(tx.samples() == rx.samples());

    double scale = std::sqrt(tx.modes() / tx.average_power());

    unsigned int errors = 0;
    for (int mode = 0; mode < tx.modes(); ++mode) {
        for (int sample = cut; sample < tx.samples() - cut; ++sample) {
            unsigned int tx_code = demodulate_16qam(scale * tx(mode, sample));
            unsigned int rx_code = demodulate_16qam(scale * rx(mode, sample));
            if (bit_dist(tx_code, rx_code) > 0) errors++;
        }
    }

    return static_cast<double>(errors) /
           static_cast<double>(tx.size() - tx.modes() * cut);
}

Field wdm_mux(const Vector<Field>& channels, const double& channel_bandwidth) {
    assert(channels.size() > 0);
    Field result(channels[0].domain());

    for (int index = 0; index < channels.size(); ++index) {
        double frequency = (index - channels.size() / 2) * channel_bandwidth;
        result += frequency_shift(channels[index], frequency);
    }

    return result;
}

Vector<Field> wdm_demux(const Field& field,
                        const int& n_channels,
                        const double& channel_bandwidth) {
    Vector<Field> channels(n_channels, Field(field.domain()));
    for (int ch = 0; ch < n_channels; ++ch)
        channels[ch] =
            wdm_select(field, ch - n_channels / 2, channel_bandwidth);

    return channels;
}

Field wdm_select(const Field& field,
                 const int& channel_index,
                 const double& channel_bandwidth) {
    Field channel = frequency_shift(field, -channel_bandwidth * channel_index);
    return bandpass_filter(channel, 0.5 * channel_bandwidth);
}

Field frequency_shift(const Field& field, const double& frequency) {
    Field shifted = field;

    double wdt = 2 * math_pi * frequency * field.dt();
    for (int mode = 0; mode < shifted.modes(); ++mode)
        for (int sample = 0; sample < shifted.samples(); ++sample)
            shifted(mode, sample) *= i_exp(wdt * sample);

    return shifted;
}

Field bandpass_filter(const Field& field, const double& frequency) {
    Field filtered = field;
    int cut = static_cast<int>(frequency * field.duration()) + 1;

    filtered.fft_inplace();
    for (int mode = 0; mode < filtered.modes(); ++mode)
        for (int sample = cut; sample <= filtered.samples() - cut; ++sample)
            filtered(mode, sample) = 0;
    filtered.ifft_inplace();

    return filtered;
}
}  // namespace focss
