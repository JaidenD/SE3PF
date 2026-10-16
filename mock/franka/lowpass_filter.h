// Mock of <franka/lowpass_filter.h> (libfranka 0.9.2). See mock/franka/README.md.
#pragma once
#include <cmath>
#include <stdexcept>

namespace franka {

    constexpr double kMaxCutoffFrequency = 1000.0;
    constexpr double kDefaultCutoffFrequency = 100.0;

    // As libfranka: first-order low-pass filter with sample time and cutoff frequency; throws std::invalid_argument
    // for a non-finite input.
    inline double lowpassFilter(double sample_time, double y, double y_last, double cutoff_frequency)
    {
        if (sample_time < 0 || !std::isfinite(sample_time))
            throw std::invalid_argument("lowpass-filter: sample_time is negative, infinite or NaN.");
        if (cutoff_frequency <= 0 || !std::isfinite(cutoff_frequency))
            throw std::invalid_argument("lowpass-filter: cutoff_frequency is zero, negative, infinite or NaN.");
        if (!std::isfinite(y) || !std::isfinite(y_last))
            throw std::invalid_argument("lowpass-filter: current or past input value of the signal to be "
                                        "filtered is infinite or NaN.");
        const double gain = sample_time / (sample_time + (1.0 / (2.0 * M_PI * cutoff_frequency)));
        return gain * y + (1 - gain) * y_last;
    }

}  // namespace franka
