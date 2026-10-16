#include "se3pf/path.hpp"
#include "se3pf/se3.hpp"

#include <cmath>
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
