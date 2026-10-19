#pragma once
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace se3pf {

    // Periodic cubic Hermite interpolation over uniform grid s_k = k L / n, k = 0, ..., n.
    // Given y_k and derivatives at the grid points, eval returns the interpolated value and derivative at s.
    // T is a double or any fixed eigen type (like Mat65 for N, Vec7 for q*).
    template <class T>
    class PeriodicHermite {
        public:
            // Compute the quantities required for interpolation, given length L and samples y, dy.
            void build(double L, std::vector<T> y, std::vector<T> dy)
            {
                if (!(L > 0.0) || y.size() < 2 || dy.size() != y.size())
                    throw std::invalid_argument("PeriodicHermite: need L > 0 and n + 1 >= 2 values and derivatives");

                n_ = static_cast<int>(y.size()) - 1;
                L_ = L;
                h_ = L / n_;
                y_ = y;
                dy_ = dy;
                y_[n_] = y_[0];
                dy_[n_] = dy_[0];
            }

            // y = interpolant at s; if dy is not null, *dy = its derivative with respect to s. Online: no allocation,
            // no exceptions. A non-finite s gives NaN in every entry of y and *dy (the weights below are NaN), never
            // an out-of-range node index. Call build first (checked by an assert in Debug builds).
            void eval(double s, T &y, T *dy = nullptr) const
            {
                assert(n_ > 0 && "PeriodicHermite::eval before build");
                const bool finite = std::isfinite(s);
                double x = finite ? std::fmod(s, L_) : 0.0;  // in (-L, L)
                if (x < 0.0)
                    x += L_;
                int k = static_cast<int>(x / h_);
                if (k >= n_)  // rounding can land exactly on L
                    k = n_ - 1;
                const double t = finite ? x / h_ - k : std::numeric_limits<double>::quiet_NaN();
                const double t2 = t * t, t3 = t2 * t;

                y = (2 * t3 - 3 * t2 + 1) * y_[k] + ((t3 - 2 * t2 + t) * h_) * dy_[k]
                  + (-2 * t3 + 3 * t2) * y_[k + 1] + ((t3 - t2) * h_) * dy_[k + 1];
                if (dy != nullptr)
                    *dy = ((6 * t2 - 6 * t) / h_) * y_[k] + (3 * t2 - 4 * t + 1) * dy_[k]
                        + ((-6 * t2 + 6 * t) / h_) * y_[k + 1] + (3 * t2 - 2 * t) * dy_[k + 1];
            }

            double length() const { return L_; }
            int intervals() const { return n_; }

            // Fourth-order periodic central differences of node values y_0..y_n (y_n = y_0) on the grid of length L:
            // y'_k ~ (y_{k-2} - 8 y_{k-1} + 8 y_{k+1} - y_{k+2}) / (12 h). For q*' and X' at the lift's nodes.
            static std::vector<T> centralDifferences(const std::vector<T> &y, double L)
            {
                const int n = static_cast<int>(y.size()) - 1;
                if (n < 4)
                    throw std::invalid_argument("PeriodicHermite::centralDifferences: need at least 5 samples");
                const double h = L / n;
                auto at = [&](int k) -> const T & { return y[((k % n) + n) % n]; };
                std::vector<T> d(y.size());
                for (int k = 0; k < n; ++k)
                    d[k] = (1.0 / (12.0 * h)) * (at(k - 2) - 8.0 * at(k - 1) + 8.0 * at(k + 1) - at(k + 2));
                d[n] = d[0];
                return d;
            }

        private:
            double L_ = 1.0, h_ = 1.0;
            int n_ = 0;
            std::vector<T> y_, dy_;
    };

}  // namespace se3pf
