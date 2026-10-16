#pragma once
#include "se3pf/types.hpp"
#include "se3pf/metric.hpp"

namespace se3pf {

    // gamma: R/LZ -> C, parameterized by arc length s wrt metric W.

    class Path {
        public:
            virtual ~Path() = default;
            virtual double length() const = 0;              // L
            virtual Pose pose(double s) const = 0;          // gamma(s)
            virtual Vec6 velocity(double s) const = 0;      // zeta_gamma(s), unit W-norm
            virtual Vec6 acceleration(double s) const = 0;  // zeta_gamma'(s)
    };

    // Example: R = R0 Rz(th) A Rz(-th), A = exp(cone * e1^), p = center + r (cos th, sin th, 0)
    class ConingCircle final : public Path {
        public:
            ConingCircle(const Matrix3d& R0, const Vector3d& center, double radius, double cone, const Metric& W);
            double length() const override;
            Pose pose(double s) const override;
            Vec6 velocity(double s) const override;
            Vec6 acceleration(double s) const override;

        private:
            Matrix3d attitude(double th) const;   // R(theta), shared by the three functions above
            Matrix3d R0_, A_;
            Vector3d c_, w0_;                     // w0 = A^T e3 - e3
            double r_;
            double sigma_;                        // constant speed in the metric: s = sigma * theta
    };

    // Example: R = exp(th e1^) exp(th e3^), p = r (cos th, sin th, 0)
    class TumblingCircle final : public Path {
        public:
            TumblingCircle(double radius, const Metric& W);
            double length() const override;
            Pose pose(double s) const override;
            Vec6 velocity(double s) const override;
            Vec6 acceleration(double s) const override;

        private:
            Matrix3d attitude(double th) const;
            double r_;
            double sigma_;
    };

}  // namespace se3pf
