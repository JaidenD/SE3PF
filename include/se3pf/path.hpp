#pragma once
#include "se3pf/types.hpp"
#include "se3pf/metric.hpp"

#include <functional>
#include <vector>

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

    // Any closed path: a smooth 2 pi-periodic pose curve g(theta) (at least three times differentiable, with
    // g'(theta) != 0), reparameterized by arc length in W numerically. The arc length s(theta) is tabulated on
    // `nodes` uniform values of theta (Gauss-Legendre on each interval) and inverted by Newton's method on its
    // cubic Hermite interpolant (relative error near 1e-12). The body velocity in theta, (g^-1 dg/dtheta)^vee, and
    // its derivative come from fourth-order central differences of the curve, so velocity(s) has unit W-norm to
    // rounding and acceleration(s) is accurate to about 1e-8. The curve must not cross itself in SE(3) (the tube
    // needs an embedded path). Construction is offline (allocates, throws std::invalid_argument for nodes < 16 or a
    // curve whose speed is not positive and finite); pose, velocity and acceleration are online: they call the
    // curve (a std::function, which does not allocate) and allocate nothing themselves.
    class SampledPath final : public Path {
        public:
            using Curve = std::function<Pose(double)>;
            SampledPath(Curve curve, const Metric& W, int nodes = 2048);
            double length() const override { return L_; }
            Pose pose(double s) const override;
            Vec6 velocity(double s) const override;
            Vec6 acceleration(double s) const override;

            double theta(double s) const;              // theta(s), the curve parameter at arc length s
            Vec6 curveVelocity(double theta) const;    // (g^-1 dg/dtheta)^vee

        private:
            Vec6 unitTangent(double theta) const;      // curveVelocity / |curveVelocity|_W
            Curve curve_;
            Metric W_;
            int n_;
            double L_ = 0.0;
            std::vector<double> s_;      // s(theta_k), theta_k = 2 pi k / n, k = 0..n
            std::vector<double> ds_;     // ds/dtheta at theta_k
    };

    // The attitude whose z-axis is the unit vector z and whose x-axis is the reference direction a projected onto
    // the plane normal to z ("no twist" about z relative to a). a must not be parallel to z.
    Matrix3d attitudeAlong(const Vector3d& z, const Vector3d& a);

    // Paper paths for the Panda (2026-11-02), as curves for SampledPath. R0 is the rest attitude (z down); its
    // x-axis is the reference direction of attitudeAlong.
    //
    // Inspection orbit: the position runs around an ellipse with a vertical wave,
    //     p(theta) = c + (rx cos theta, ry sin theta, h sin 2 theta),
    // and the tool z-axis always points at the target c - (0, 0, d), as a camera inspecting a part.
    // The defaults (sized 2026-11-02 for the Panda with the rest-posture lift): 52 deg of rotation over the loop,
    // joint margin 0.41 rad (joint 4), sigma_min(J_b) 0.114, tube radius 0.093, joint speeds at most 12.5 % of
    // their limits at 0.05 per second. Nearer the base (centre x 0.45) joint 4 comes within 0.03 rad of its limit.
    struct InspectionOrbit {
        Vector3d center = Vector3d(0.50, 0.0, 0.50);
        double rx = 0.12, ry = 0.10;   // m
        double wave = 0.03;            // h, m
        double depth = 0.25;           // d, m: the target below the centre
    };
    SampledPath::Curve inspectionOrbit(const Matrix3d& R0, const InspectionOrbit& o);
    //
    // Dome figure-eight: the tool stays normal to a dome (a sphere of radius rho around c), as in polishing, and its
    // contact point traces a figure-eight over the top: the normal is u = (X, Y, 1) / |(X, Y, 1)| with the
    // lemniscate X = tx sin theta, Y = ty sin theta cos theta, p = c + rho u, tool z = -u. The tool also twists
    // about its axis by psi0 cos theta, so the two passes through the top (theta = 0, pi) differ by 2 psi0 and the
    // path does not cross itself in SE(3).
    // The defaults (sized 2026-11-02 like the orbit): 92 deg of rotation over the loop, joint margin 0.34 rad
    // (joint 4), sigma_min(J_b) 0.101, tube radius 0.053 (set by the ends of the lobes, not by the crossing),
    // joint speeds at most 14.4 % of their limits at 0.05 per second. Not with PandaExample::hardwareLiftWeights()
    // (weight 10 on joint 7): the twist needs joint 7, and that lift fails part way round.
    struct DomeFigureEight {
        Vector3d center = Vector3d(0.50, 0.0, 0.27);
        double radius = 0.25;          // rho, m
        double tx = 0.60, ty = 0.60;   // tan of the tilt: up to atan(tx) = 31 deg about y, atan(ty / 2) = 17 deg about x
        double twist = 0.8;            // psi0, rad
    };
    SampledPath::Curve domeFigureEight(const Matrix3d& R0, const DomeFigureEight& d);
    //
    // Straight rail (2026-11-02): a long, thin, smooth racetrack in a horizontal plane at constant attitude R0, for
    // pushing the end effector along a line by hand. With u = theta - phase,
    //     p = c + (-w tanh(k cos u) / tanh(k), a sin u, 0),
    // the two long sides are nearly straight lines along y at x = c_x -+ w (within w (1 - tanh(k/2)/tanh(k)) where
    // |cos u| > 1/2) and the turns are at y = c_y +- a. phase puts s = PandaExample::s0 = 0.1 near the middle of the
    // near side (x = c_x - w): sin(phase) = 0.1 / a. The tube radius is at most w (the sides are 2w apart).
    // The defaults (sized 2026-11-02): the near side is straight to 1.4 mm over +-10 cm and 4.4 mm over +-15 cm
    // around its middle (where s0 is); tube radius 0.057 (k = 4 and w = 0.06 gave 0.017: the turns limit it); joint
    // margin 0.57 rad (joint 4), sigma_min(J_b) 0.162, joint speeds at most 8.5 % of their limits at 0.05 per second.
    struct StraightRail {
        Vector3d center = Vector3d(0.55, 0.0, 0.35);
        double halfLength = 0.25;   // a, m: the turns are at y = +-a
        double halfWidth = 0.10;    // w, m
        double sharpness = 2.0;     // k: larger gives straighter sides and tighter turns
    };
    SampledPath::Curve straightRail(const Matrix3d& R0, const StraightRail& r);

}  // namespace se3pf
