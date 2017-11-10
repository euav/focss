#include "fiber.h"

namespace focss {
Fiber Fiber::linear() const {
    Fiber new_fiber = *this;
    new_fiber.gamma = 0;
    return new_fiber;
}

Fiber Fiber::nonlinear() const {
    Fiber new_fiber = *this;
    new_fiber.alpha = 0;
    new_fiber.beta2 = 0;
    return new_fiber;
}

Fiber Fiber::dispersive() const {
    Fiber new_fiber = *this;
    new_fiber.alpha = 0;
    new_fiber.gamma = 0;
    return new_fiber;
}

Fiber Fiber::operator-() const {
    Fiber new_fiber = *this;
    new_fiber.alpha = -alpha;
    new_fiber.beta2 = -beta2;
    new_fiber.gamma = -gamma;
    return new_fiber;
}

Fiber Fiber::operator*(const double& length_factor) const {
    Fiber new_fiber = *this;

    if (length_factor >= 0) {
        new_fiber.length *= length_factor;
    } else {
        new_fiber.alpha = -alpha;
        new_fiber.beta2 = -beta2;
        new_fiber.gamma = -gamma;
        new_fiber.length *= -length_factor;
    }

    return new_fiber;
}
}  // namespace focss
