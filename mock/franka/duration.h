// Mock of <franka/duration.h> (libfranka 0.9.2): the subset the driver uses. See mock/franka/README.md.
#pragma once
#include <chrono>
#include <cstdint>

namespace franka {

    class Duration {
        public:
            Duration() noexcept = default;
            explicit Duration(uint64_t milliseconds) noexcept : ms_(milliseconds) {}
            Duration(std::chrono::duration<uint64_t, std::milli> duration) noexcept : ms_(duration.count()) {}
            double toSec() const noexcept { return static_cast<double>(ms_) / 1000.0; }
            uint64_t toMSec() const noexcept { return ms_; }
            Duration operator+(const Duration& rhs) const noexcept { return Duration(ms_ + rhs.ms_); }
            Duration& operator+=(const Duration& rhs) noexcept
            {
                ms_ += rhs.ms_;
                return *this;
            }
            Duration operator-(const Duration& rhs) const noexcept { return Duration(ms_ - rhs.ms_); }
            bool operator==(const Duration& rhs) const noexcept { return ms_ == rhs.ms_; }
            bool operator!=(const Duration& rhs) const noexcept { return ms_ != rhs.ms_; }
            bool operator<(const Duration& rhs) const noexcept { return ms_ < rhs.ms_; }

        private:
            uint64_t ms_ = 0;
    };

}  // namespace franka
