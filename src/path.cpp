#include "se3pf/path.hpp"
#include "se3pf/se3.hpp"

#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>

using namespace Eigen;

namespace se3pf {
    namespace {
        // Rotation by th about the z-axis.
        Matrix3d Rz(double th)
        {
            return expSO3(Vector3d(0.0, 0.0, th));
        }

        // Stack (omega, v) into a twist.
        Vec6 twist(const Vector3d &w, const Vector3d &v)
        {
            Vec6 z;
            z << w, v;
            return z;
        }

        // Speed sigma = |zeta~(0)|_W, checked to be constant along the curve: the s <-> theta conversion below
        // assumes it. Both curves have constant speed in every product-type W (the only kind Metric accepts), so
        // this is a safeguard for a new curve or a future general W.
        template <class BodyVelocity>
        double constantSpeed(const BodyVelocity &zetaTilde, const Metric &W)
        {
            const double sigma = W.norm(zetaTilde(0.0));
            for (int k = 1; k < 16; ++k) {
                const double th = 2.0 * EIGEN_PI * k / 16.0;
                if (std::abs(W.norm(zetaTilde(th)) - sigma) > 1e-12 * sigma)
                    throw std::invalid_argument("Path: the speed is not constant in this metric");
            }
            return sigma;
        }
    }  // namespace

    // ---------------------------------------------------------------------------------------------------------
    // Coning circle: R = R0 Rz(th) A Rz(-th), A = exp(cone e1^), p = c + r (cos th, sin th, 0).
    // Body velocity in theta: omega~ = Rz(th) w0 with w0 = A^T e3 - e3, v~ = R^T p'.

    ConingCircle::ConingCircle(const Matrix3d &R0, const Vector3d &center, double radius, double cone,
                               const Metric &W)
        : R0_(R0), A_(expSO3(Vector3d(cone, 0.0, 0.0))), c_(center), r_(radius)
    {
        w0_ = A_.transpose() * Vector3d::UnitZ() - Vector3d::UnitZ();
        auto zetaTilde = [this](double th) {
            const Vector3d pd = r_ * Vector3d(-std::sin(th), std::cos(th), 0.0);
            return twist(Rz(th) * w0_, attitude(th).transpose() * pd);
        };
        sigma_ = constantSpeed(zetaTilde, W);
    }

    double ConingCircle::length() const
    {
        return 2.0 * EIGEN_PI * sigma_;
    }

    Matrix3d ConingCircle::attitude(double th) const
    {
        return R0_ * Rz(th) * A_ * Rz(-th);
    }

    Pose ConingCircle::pose(double s) const
    {
        const double th = s / sigma_;
        return makePose(attitude(th), c_ + r_ * Vector3d(std::cos(th), std::sin(th), 0.0));
    }

    Vec6 ConingCircle::velocity(double s) const
    {
        const double th = s / sigma_;
        const Vector3d pd = r_ * Vector3d(-std::sin(th), std::cos(th), 0.0);  // dp/dtheta
        return twist(Rz(th) * w0_, attitude(th).transpose() * pd) / sigma_;
    }

    Vec6 ConingCircle::acceleration(double s) const
    {
        const double th = s / sigma_;
        const Matrix3d R = attitude(th);
        const Vector3d w = Rz(th) * w0_;                                 // omega~
        const Vector3d dw = Rz(th) * Vector3d::UnitZ().cross(w0_);       // d omega~ / d theta
        const Vector3d pd = r_ * Vector3d(-std::sin(th), std::cos(th), 0.0);
        const Vector3d pdd = -r_ * Vector3d(std::cos(th), std::sin(th), 0.0);
        const Vector3d v = R.transpose() * pd;                           // v~
        const Vector3d dv = R.transpose() * pdd - w.cross(v);            // d v~ / d theta
        return twist(dw, dv) / (sigma_ * sigma_);
    }

    // ---------------------------------------------------------------------------------------------------------
    // Tumbling circle: R = exp(th e1^) exp(th e3^), p = r (cos th, sin th, 0).
    // Body velocity in theta: omega~ = (cos th, -sin th, 1), v~ = R^T p'.

    TumblingCircle::TumblingCircle(double radius, const Metric &W) : r_(radius)
    {
        auto zetaTilde = [this](double th) {
            const Vector3d pd = r_ * Vector3d(-std::sin(th), std::cos(th), 0.0);
            return twist(Vector3d(std::cos(th), -std::sin(th), 1.0), attitude(th).transpose() * pd);
        };
        sigma_ = constantSpeed(zetaTilde, W);
    }

    double TumblingCircle::length() const
    {
        return 2.0 * EIGEN_PI * sigma_;
    }

    Matrix3d TumblingCircle::attitude(double th) const
    {
        return expSO3(Vector3d(th, 0.0, 0.0)) * expSO3(Vector3d(0.0, 0.0, th));
    }

    Pose TumblingCircle::pose(double s) const
    {
        const double th = s / sigma_;
        return makePose(attitude(th), r_ * Vector3d(std::cos(th), std::sin(th), 0.0));
    }

    Vec6 TumblingCircle::velocity(double s) const
    {
        const double th = s / sigma_;
        const Vector3d pd = r_ * Vector3d(-std::sin(th), std::cos(th), 0.0);
        return twist(Vector3d(std::cos(th), -std::sin(th), 1.0), attitude(th).transpose() * pd) / sigma_;
    }

    Vec6 TumblingCircle::acceleration(double s) const
    {
        const double th = s / sigma_;
        const Matrix3d R = attitude(th);
        const Vector3d w(std::cos(th), -std::sin(th), 1.0);
        const Vector3d dw(-std::sin(th), -std::cos(th), 0.0);
        const Vector3d pd = r_ * Vector3d(-std::sin(th), std::cos(th), 0.0);
        const Vector3d pdd = -r_ * Vector3d(std::cos(th), std::sin(th), 0.0);
        const Vector3d v = R.transpose() * pd;
        const Vector3d dv = R.transpose() * pdd - w.cross(v);
        return twist(dw, dv) / (sigma_ * sigma_);
    }

}  // namespace se3pf

namespace se3pf {

    // ---------------------------------------------------------------------------------------------------------
    // SampledPath

    namespace {
        constexpr double kTwoPi = 2.0 * EIGEN_PI;
        constexpr double kH = 1e-3;  // step of the central differences in theta

        // Fourth-order central difference of f at x with step h
        template <class F>
        auto centralDiff(const F& f, double x, double h)
        {
            return ((f(x - 2.0 * h) - f(x + 2.0 * h)) + 8.0 * (f(x + h) - f(x - h))) / (12.0 * h);
        }
    }  // namespace

    SampledPath::SampledPath(Curve curve, const Metric& W, int nodes) : curve_(std::move(curve)), W_(W), n_(nodes)
    {
        if (n_ < 16 || !curve_)
            throw std::invalid_argument("SampledPath: need a curve and at least 16 nodes");
        auto speed = [this](double th) { return W_.norm(curveVelocity(th)); };
        // Gauss-Legendre with 4 points on each interval
        const double gx[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
        const double gw[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
        const double h = kTwoPi / n_;
        s_.assign(n_ + 1, 0.0);
        ds_.assign(n_ + 1, 0.0);
        for (int k = 0; k <= n_; ++k) {
            ds_[k] = speed(k * h);
            if (!std::isfinite(ds_[k]) || !(ds_[k] > 0.0))
                throw std::invalid_argument("SampledPath: the curve's speed must be positive and finite");
            if (k == n_)
                break;
            double acc = 0.0;
            for (int i = 0; i < 4; ++i)
                acc += gw[i] * speed((k + 0.5 + 0.5 * gx[i]) * h);
            s_[k + 1] = s_[k] + 0.5 * h * acc;
        }
        L_ = s_[n_];
    }

    double SampledPath::theta(double s) const
    {
        s = s - L_ * std::floor(s / L_);
        if (!std::isfinite(s))
            return std::numeric_limits<double>::quiet_NaN();
        const double h = kTwoPi / n_;
        int k = static_cast<int>(std::upper_bound(s_.begin(), s_.end(), s) - s_.begin()) - 1;
        k = std::clamp(k, 0, n_ - 1);
        // Cubic Hermite s(tau) on [theta_k, theta_k+1], tau in [0, 1], solved for s by Newton's method
        const double s0 = s_[k], s1 = s_[k + 1], m0 = ds_[k] * h, m1 = ds_[k + 1] * h;
        double tau = (s - s0) / (s1 - s0);
        for (int it = 0; it < 8; ++it) {
            const double t2 = tau * tau, t3 = t2 * tau;
            const double f = (2 * t3 - 3 * t2 + 1) * s0 + (t3 - 2 * t2 + tau) * m0 + (-2 * t3 + 3 * t2) * s1
                             + (t3 - t2) * m1 - s;
            const double df = (6 * t2 - 6 * tau) * s0 + (3 * t2 - 4 * tau + 1) * m0 + (-6 * t2 + 6 * tau) * s1
                              + (3 * t2 - 2 * tau) * m1;
            const double step = f / df;
            tau -= step;
            if (std::abs(step) < 1e-15)
                break;
        }
        return (k + tau) * h;
    }

    Vec6 SampledPath::curveVelocity(double th) const
    {
        const Pose g = curve_(th);
        const Matrix4d dg = centralDiff([this](double x) { return Matrix4d(curve_(x)); }, th, kH);
        const Matrix3d R = g.block<3, 3>(0, 0);
        const Matrix3d Om = R.transpose() * dg.block<3, 3>(0, 0);
        Vec6 z;
        z << vee(Matrix3d(0.5 * (Om - Om.transpose()))), R.transpose() * dg.block<3, 1>(0, 3);
        return z;
    }

    Vec6 SampledPath::unitTangent(double th) const
    {
        const Vec6 z = curveVelocity(th);
        return z / W_.norm(z);
    }

    Pose SampledPath::pose(double s) const
    {
        return curve_(theta(s));
    }

    Vec6 SampledPath::velocity(double s) const
    {
        return unitTangent(theta(s));
    }

    Vec6 SampledPath::acceleration(double s) const
    {
        const double th = theta(s);
        const Vec6 du = centralDiff([this](double x) { return Vec6(unitTangent(x)); }, th, kH);
        return du / W_.norm(curveVelocity(th));
    }

    // ---------------------------------------------------------------------------------------------------------
    // Paper paths

    Matrix3d attitudeAlong(const Vector3d& z, const Vector3d& a)
    {
        const Vector3d zu = z.normalized();
        const Vector3d x = (a - a.dot(zu) * zu).normalized();
        Matrix3d R;
        R << x, zu.cross(x), zu;
        return R;
    }

    SampledPath::Curve inspectionOrbit(const Matrix3d& R0, const InspectionOrbit& o)
    {
        const Vector3d a = R0.col(0);
        return [a, o](double th) {
            const Vector3d p = o.center + Vector3d(o.rx * std::cos(th), o.ry * std::sin(th), o.wave * std::sin(2.0 * th));
            const Vector3d target = o.center - Vector3d(0.0, 0.0, o.depth);
            return makePose(attitudeAlong(target - p, a), p);
        };
    }

    SampledPath::Curve straightRail(const Matrix3d& R0, const StraightRail& r)
    {
        const double phase = std::asin(std::min(1.0, 0.1 / r.halfLength));
        return [R0, r, phase](double th) {
            const double u = th - phase;
            const double x = -r.halfWidth * std::tanh(r.sharpness * std::cos(u)) / std::tanh(r.sharpness);
            return makePose(R0, Vector3d(r.center + Vector3d(x, r.halfLength * std::sin(u), 0.0)));
        };
    }

    SampledPath::Curve domeFigureEight(const Matrix3d& R0, const DomeFigureEight& d)
    {
        const Vector3d a = R0.col(0);
        return [a, d](double th) {
            const Vector3d u = Vector3d(d.tx * std::sin(th), d.ty * std::sin(th) * std::cos(th), 1.0).normalized();
            const Matrix3d R = attitudeAlong(-u, a) * expSO3(Vector3d(0.0, 0.0, d.twist * std::cos(th)));
            return makePose(R, Vector3d(d.center + d.radius * u));
        };
    }

}  // namespace se3pf
