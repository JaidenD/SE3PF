// Mock of <franka/rate_limiting.h> (libfranka 0.9.2): the torque-rate limiter of the torque interface. See
// mock/franka/README.md.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace franka {

    constexpr double kDeltaT = 1e-3;
    constexpr double kLimitEps = 1e-3;
    constexpr std::array<double, 7> kMaxTorqueRate{{1000 - kLimitEps, 1000 - kLimitEps, 1000 - kLimitEps,
                                                    1000 - kLimitEps, 1000 - kLimitEps, 1000 - kLimitEps,
                                                    1000 - kLimitEps}};

    inline std::array<double, 7> limitRate(const std::array<double, 7>& max_derivatives,
                                           const std::array<double, 7>& commanded_values,
                                           const std::array<double, 7>& last_commanded_values)
    {
        for (double d : commanded_values)
            if (!std::isfinite(d))
                throw std::invalid_argument("Commanding value is infinite or NaN.");
        std::array<double, 7> limited{};
        for (std::size_t i = 0; i < 7; ++i) {
            const double derivative = (commanded_values[i] - last_commanded_values[i]) / kDeltaT;
            limited[i] = last_commanded_values[i]
                         + std::max(std::min(derivative, max_derivatives[i]), -max_derivatives[i]) * kDeltaT;
        }
        return limited;
    }

}  // namespace franka
