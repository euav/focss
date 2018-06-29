#include <algorithm>
#include "linear.h"

namespace focss {
HomothetyEqualizer::HomothetyEqualizer() : scalar_(1.0) {}

void HomothetyEqualizer::train(const Field& desired, const Field& actual) {
    scalar_ = estimate(desired, actual);
}

Field HomothetyEqualizer::equalize(const Field& original) const {
    return scalar_ * original;
}

Field HomothetyEqualizer::train_equalize(const Field& desired,
                                         const Field& actual) const {
    return estimate(desired, actual) * actual;
}

double HomothetyEqualizer::get_scalar() const { return scalar_; }

double HomothetyEqualizer::estimate(const Field& desired,
                                    const Field& actual) const {
    complex_t xx_var = 0.0;
    complex_t xy_var = 0.0;

    for (int i = 0; i < actual.size(); ++i)
        xx_var += norm(actual[i]);
    for (int i = 0; i < actual.size(); ++i)
        xy_var += conj(actual[i]) * desired[i];

    return std::abs(xy_var / xx_var);
}
}  // namespace focss