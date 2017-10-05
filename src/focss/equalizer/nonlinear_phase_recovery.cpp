#include "focss/equalizer/nonlinear_phase_recovery.h"

double median(std::vector<double>& numbers) {
    std::sort(numbers.begin(), numbers.end());

    unsigned long size = numbers.size();
    if (size % 2 == 0) 
        return 0.5 * (numbers[size / 2] + numbers[size / 2 - 1]);
    else
        return numbers[size / 2];
}

NonlinearPhaseRecovery::NonlinearPhaseRecovery() : trained(false) {}

void NonlinearPhaseRecovery::train(const Field& desired, const Field& actual) {
    if (desired.size() == actual.size()) {
        unsigned long size = actual.size();
        std::vector<double> x(size);
        std::vector<double> y(size);
        for (unsigned long i = 0; i < size; ++i) {
            x[i] = norm(actual[i]);
            y[i] = arg(desired[i] / actual[i]);
        }

        std::vector<double> slopes;
        for (unsigned long i = 0; i < size; ++i)
            for (unsigned long j = i + 1; j < size; ++j)
                slopes.push_back((y[i] - y[j]) / (x[i] - x[j]));
        
        b = median(slopes);

        std::vector<double> intercepts(size);
        for (unsigned long i = 0; i < size; ++i)
            intercepts[i] = y[i] - b * x[i];

        a = median(intercepts);

        trained = true;
    } else {
        a = 0;
        b = 0;

        trained = false;
    }
}

RealVector NonlinearPhaseRecovery::getWeights() const {
    return {a, b};
}

Field NonlinearPhaseRecovery::equalize(const Field& original) const {
    if (!trained) return original;

    Field equalized = original;
    for (unsigned long i = 0; i < original.size(); ++i)
        equalized[i] *= i_exp(a + b * norm(original[i]));

    return equalized;
}
