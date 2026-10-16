#include "se3pf/metric.hpp"
#include "se3pf/se3.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace se3pf {
    // Only the product type has the closed forms below. A W that is not, e.g. diag(1, 1, 2, 1, 1, 1), would be
    // accepted by everything else and silently give the wrong geodesics, so refuse it here.
    Metric::Metric(const Mat6& W) : W_(W)
    {
        const double wR = W(0, 0), wp = W(3, 3);
        Mat6 P = Mat6::Zero();
        P.diagonal() << wR, wR, wR, wp, wp, wp;
        if (!(wR > 0.0) || !(wp > 0.0) || !W.allFinite() || (W - P).cwiseAbs().maxCoeff() > 1e-12 * std::max(wR, wp))
            throw std::invalid_argument("Metric: only the product type W = diag(wR I, wp I), wR, wp > 0, is "
                                        "implemented (Example 1)");
    }

    Metric Metric::product(double wR, double wp)
    {
        Mat6 W = Mat6::Identity();
        W.diagonal().head<3>() *= wR;
        W.diagonal().tail<3>() *= wp;
        return Metric(W);
    }

    double Metric::inner(const Vec6 &v, const Vec6 &w) const
    {
        return v.dot(W_ * w);
    }

    Vec6 Metric::flat(const Vec6 &v) const
    {
        return W_ * v;
    }

    double Metric::norm(const Vec6 &v) const
    {
        return std::sqrt(inner(v, v));
    }

    // Example 1: Exp_I(z) = (exp(hat(omega)), v). The translation is v itself, not the screw motion's J_l(omega) v.
    Pose Metric::Exp(const Vec6 &z) const
    {
        const Vector3d w = z.head<3>();
        const Vector3d v = z.tail<3>();
        return makePose(expSO3(w), v);
    }

    // Example 1: Log_I(g) = (log(R), p) when tr(R) != -1, i.e. rotation angle below pi.
    bool Metric::Log(const Pose &g, Vec6 &z) const
    {
        const Matrix3d R = g.block<3,3>(0,0);
        if (R.trace() < -1.0 + 1e-12)  // angle pi: two minimizing geodesics, no unique logarithm
            return false;
        z.head<3>() = logSO3(R);
        z.tail<3>() = g.block<3,1>(0,3);
        return true;
    }

    // Example 1: Gamma(z) = diag(Jr(omega), exp(hat(omega))^T).
    Mat6 Metric::Gamma(const Vec6 &z) const
    {
        const Vector3d w = z.head<3>();
        Mat6 G = Mat6::Zero();
        G.block<3,3>(0,0) = Jr(w);
        G.block<3,3>(3,3) = expSO3(w).transpose();
        return G;
    }

    Mat6 Metric::GammaInv(const Vec6 &z) const
    {
        const Vector3d w = z.head<3>();
        Mat6 G = Mat6::Zero();
        G.block<3,3>(0,0) = JrInv(w);
        G.block<3,3>(3,3) = expSO3(w);
        return G;
    }

    double Metric::distance(const Pose &g1, const Pose &g2) const
    {
        Vec6 z;
        if (!Log(inverse(g1) * g2, z))
            return std::numeric_limits<double>::quiet_NaN();
        return norm(z);
    }
}
