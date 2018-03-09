#include "focss/equalizer/nonlinear_phase_recovery.h"
#include <algorithm>
#include "focss/functions.h"

namespace focss {
double median(RealVector numbers) {
    std::sort(numbers.raw(), numbers.raw() + numbers.size());

    int size = numbers.size();
    if (size % 2 == 0)
        return 0.5 * (numbers[size / 2] + numbers[size / 2 - 1]);
    else
        return numbers[size / 2];
}

NonlinearPhaseRecovery::NonlinearPhaseRecovery() : trained_(false) {}

void NonlinearPhaseRecovery::train(const Field& desired, const Field& actual) {
    if (desired.size() == actual.size()) {
        int size = actual.size();
        RealVector x(size);
        RealVector y(size);
        for (int i = 0; i < size; ++i) {
            x[i] = norm(actual[i]);
            y[i] = arg(desired[i] / actual[i]);
        }

        int index = 0;
        RealVector slopes((size * (size - 1)) / 2);
        for (int j = 0; j < size; ++j)
            for (int i = 0; i < j; ++i)
                slopes[(j * (j - 1)) / 2 + i] = (y[i] - y[j]) / (x[i] - x[j]);

        b_ = median(slopes);

        RealVector intercepts(size);
        for (int i = 0; i < size; ++i)
            intercepts[i] = y[i] - b_ * x[i];

        a_ = median(intercepts);

        trained_ = true;
    } else {
        a_ = 0;
        b_ = 0;

        trained_ = false;
    }
}

Field NonlinearPhaseRecovery::equalize(const Field& original) const {
    if (!trained_) return original;

    Field equalized = original;
    for (int index = 0; index < original.size(); ++index)
        equalized[index] *= i_exp(a_ + b_ * norm(original[index]));

    return equalized;
}

RealVector NonlinearPhaseRecovery::get_weights() const { return {a_, b_}; }
}  // namespace focss
