#pragma once
#include "se3pf/types.hpp"

namespace se3pf {
    // Left-invariant metric on SE(3) of form W = diag(wR I, wp I) (eq:metric, Example 1). Exp, Log and
    // Gamma are the closed forms of Example 1, which hold for this type only.
    class Metric {
        public:
            // Throws std::invalid_argument unless W = diag(wR, wR, wR, wp, wp, wp) with wR > 0 and wp > 0 (to a
            // relative 1e-12).
            explicit Metric(const Mat6& W);
            static Metric product(double wR, double wp);

            double inner(const Vec6& v, const Vec6& w) const;       // <v, w>_W = v^T W w
            Vec6 flat(const Vec6& v) const;                         // flat map: W v, the covector <v, .>_W
            double norm(const Vec6& v) const;                       // |v|_W
            const Mat6& matrix() const { return W_; }

            double operator()(const Vec6& v, const Vec6& w) const { return inner(v, w); }
            Vec6 operator()(const Vec6& v) const { return flat(v); }

            Pose Exp(const Vec6& z) const;                          // Exp_I
            bool Log(const Pose& g, Vec6& z) const;                 // Log_I; false if the rotation angle is pi
            Mat6 Gamma(const Vec6& z) const;                        // derivative of Exp_I (eq:Gamma)
            Mat6 GammaInv(const Vec6& z) const;
            double distance(const Pose& g1, const Pose& g2) const;  // |Log_I(g1^-1 g2)|_W (eq:dist); NaN if Log fails

        private:
            Mat6 W_;
    };
}
