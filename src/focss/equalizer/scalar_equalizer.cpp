#include <algorithm>
#include "linear.h"

namespace focss {
ScalarEqualizer::ScalarEqualizer() : scalar_(1.0) {}

void ScalarEqualizer::train(const Field& desired, const Field& actual) {
    scalar_ = estimate(desired, actual);
}

Field ScalarEqualizer::equalize(const Field& original) const {
    return scalar_ * original;
}

Field ScalarEqualizer::train_equalize(const Field& desired,
                                      const Field& actual) const {
    return estimate(desired, actual) * actual;
}

complex_t ScalarEqualizer::get_scalar() const { return scalar_; }

complex_t ScalarEqualizer::estimate(const Field& desired,
                                    const Field& actual) const {
    complex_t xx_var = 0.0;
    complex_t xy_var = 0.0;

    for (int i = 0; i < actual.size(); ++i)
        xx_var += norm(actual[i]);
    for (int i = 0; i < actual.size(); ++i)
        xy_var += conj(actual[i]) * desired[i];

    return xy_var / xx_var;
}
}  // namespace focss