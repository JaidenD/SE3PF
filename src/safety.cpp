#include "se3pf/safety.hpp"
#include "se3pf/controller.hpp"  // quinticRamp

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace se3pf {

    const char* toString(Violation v) noexcept
    {
        switch (v) {
            case Violation::None: return "none";
            case Violation::Projection: return "controller step failed (left the tube, A singular or non-finite)";
            case Violation::JointLimit: return "joint position margin (0.10 rad from a limit)";
            case Violation::JointVelocity: return "joint velocity margin";
            case Violation::Acceleration: return "commanded joint acceleration margin";
            case Violation::Torque: return "joint torque margin (delivered torque, gravity included)";
            case Violation::TorqueRate: return "commanded torque rate margin";
            case Violation::FlangeSpeed: return "flange speed margin";
            case Violation::Tube: return "transverse error |xi| above its bound";
            case Violation::Redundancy: return "redundancy output |beta| above its bound";
            case Violation::Conditioning: return "A(q) badly conditioned (sigma_min below its bound)";
            case Violation::NonFinite: return "non-finite value";
            case Violation::Timeout: return "not at rest within the duration + 5 s";
        }
        return "unknown";
    }

    namespace {
        // A requested limit must be a finite number >= 0 (offline check; throws).
        void requireValid(double value, const char* name)
        {
            if (!std::isfinite(value) || value < 0.0)
                throw std::invalid_argument(std::string("SafetyMonitor: ") + name + " must be finite and >= 0");
        }

        // "At most" bounds (fractions, speeds, xi, beta): the smaller of the request and the bound.
        double atMost(double requested, double bound, bool& clamped)
        {
            if (requested > bound) {
                clamped = true;
                return bound;
            }
            return requested;
        }

        // "At least" bounds (joint margin, sigma_min): the larger of the request and the bound.
        double atLeast(double requested, double bound, bool& clamped)
        {
            if (requested < bound) {
                clamped = true;
                return bound;
            }
            return requested;
        }
    }  // namespace

    SafetyMonitor::SafetyMonitor() : SafetyMonitor(SafetyLimits()) {}

    SafetyMonitor::SafetyMonitor(const SafetyLimits& r)
    {
        requireValid(r.jointMargin, "jointMargin");
        requireValid(r.dqFrac, "dqFrac");
        requireValid(r.ddqFrac, "ddqFrac");
        requireValid(r.tauFrac, "tauFrac");
        requireValid(r.tauRateFrac, "tauRateFrac");
        requireValid(r.flangeSpeed, "flangeSpeed");
        requireValid(r.xi, "xi");
        requireValid(r.beta, "beta");
        requireValid(r.sigmaMin, "sigmaMin");

        using B = HardwareBounds;
        bool c = false;
        limits_.jointMargin = atLeast(r.jointMargin, B::jointMargin, c);
        limits_.dqFrac = atMost(r.dqFrac, B::dqFrac, c);
        limits_.ddqFrac = atMost(r.ddqFrac, B::ddqFrac, c);
        limits_.tauFrac = atMost(r.tauFrac, B::tauFrac, c);
        limits_.tauRateFrac = atMost(r.tauRateFrac, B::tauRateFrac, c);
        limits_.flangeSpeed = atMost(r.flangeSpeed, B::flangeSpeed, c);
        limits_.xi = atMost(r.xi, B::xi, c);
        limits_.beta = atMost(r.beta, B::beta, c);
        limits_.sigmaMin = atLeast(r.sigmaMin, B::sigmaMin, c);
        clamped_ = c;

        const Vec7 margin = Vec7::Constant(limits_.jointMargin);
        qLow_ = PandaArm::qMin() + margin;
        qHigh_ = PandaArm::qMax() - margin;
        dqLimit_ = limits_.dqFrac * PandaArm::dqMax();
        ddqLimit_ = limits_.ddqFrac * PandaArm::ddqMax();
        tauLimit_ = limits_.tauFrac * PandaArm::tauMax();
        tauRateLimit_ = limits_.tauRateFrac * PandaArm::tauRateMax;
    }

    Violation SafetyMonitor::checkOutputs(const Vec5& xi, double beta, double sigmaMin) const noexcept
    {
        if (!xi.allFinite() || !std::isfinite(beta) || !std::isfinite(sigmaMin))
            return Violation::NonFinite;
        if (!(xi.norm() <= limits_.xi))
            return Violation::Tube;
        if (!(std::abs(beta) <= limits_.beta))
            return Violation::Redundancy;
        if (!(sigmaMin >= limits_.sigmaMin))
            return Violation::Conditioning;
        return Violation::None;
    }

    Violation SafetyMonitor::checkState(const Vec7& q, const Vec7& dq) const noexcept
    {
        if (!q.allFinite() || !dq.allFinite())
            return Violation::NonFinite;
        for (int i = 0; i < 7; ++i)
            if (!(q[i] >= qLow_[i] && q[i] <= qHigh_[i]))
                return Violation::JointLimit;
        for (int i = 0; i < 7; ++i)
            if (!(std::abs(dq[i]) <= dqLimit_[i]))
                return Violation::JointVelocity;
        // Translational speed of the flange: the v part of the body twist J_b(q) qd (its norm is the same in every
        // frame).
        Mat67 Jb = Mat67::Zero();
        PandaArm::output(q, &Jb);
        const Vec6 zeta = Jb * dq;
        if (!(zeta.tail<3>().norm() <= limits_.flangeSpeed))
            return Violation::FlangeSpeed;
        return Violation::None;
    }

    Violation SafetyMonitor::checkAcceleration(const Vec7& qdd) const noexcept
    {
        if (!qdd.allFinite())
            return Violation::NonFinite;
        for (int i = 0; i < 7; ++i)
            if (!(std::abs(qdd[i]) <= ddqLimit_[i]))
                return Violation::Acceleration;
        return Violation::None;
    }

    Violation SafetyMonitor::checkTorque(const Vec7& tauDelivered, const Vec7& tauCmd, const Vec7& tauCmdPrev,
                                         double dt) const noexcept
    {
        if (!tauDelivered.allFinite() || !tauCmd.allFinite() || !tauCmdPrev.allFinite() || !std::isfinite(dt)
            || !(dt > 0.0))
            return Violation::NonFinite;
        for (int i = 0; i < 7; ++i)
            if (!(std::abs(tauDelivered[i]) <= tauLimit_[i]))
                return Violation::Torque;
        for (int i = 0; i < 7; ++i)
            if (!(std::abs(tauCmd[i] - tauCmdPrev[i]) <= tauRateLimit_ * dt))
                return Violation::TorqueRate;
        return Violation::None;
    }

    // ---------------------------------------------------------------------------------------------------------
    // ApproachMotion

    ApproachMotion::ApproachMotion(const Vec7& q0, const Vec7& q1) : q0_(q0), q1_(q1)
    {
        if (!q0.allFinite() || !q1.allFinite())
            throw std::invalid_argument("ApproachMotion: the configurations must be finite");
        // Peaks of the smoothstep's derivatives: max r' = 15/8, max |r''| = 10/sqrt(3), max |r'''| = 60, so joint i
        // moving by d_i peaks at 15/8 d_i / T, 10/sqrt(3) d_i / T^2 and 60 d_i / T^3.
        const Vec7 d = (q1 - q0).cwiseAbs();
        const Vec7 dqMax = PandaArm::dqMax(), ddqMax = PandaArm::ddqMax(), jerkMax = PandaArm::jerkMax();
        double T = kMinTime;
        for (int i = 0; i < 7; ++i) {
            T = std::max(T, 1.875 * d[i] / (kVelocityFraction * dqMax[i]));
            T = std::max(T, std::sqrt(10.0 / std::sqrt(3.0) * d[i] / (kAccelerationFraction * ddqMax[i])));
            T = std::max(T, std::cbrt(60.0 * d[i] / (kJerkFraction * jerkMax[i])));
        }
        T_ = T;
    }

    Vec7 ApproachMotion::position(double t) const
    {
        if (t >= T_)
            return q1_;
        return q0_ + quinticRamp(t / T_) * (q1_ - q0_);
    }

    Vec7 ApproachMotion::velocity(double t) const
    {
        const double x = std::clamp(t / T_, 0.0, 1.0);
        return (30.0 * x * x * (1.0 - x) * (1.0 - x) / T_) * (q1_ - q0_);
    }

    Vec7 ApproachMotion::acceleration(double t) const
    {
        const double x = std::clamp(t / T_, 0.0, 1.0);
        return (60.0 * x * (1.0 - x) * (1.0 - 2.0 * x) / (T_ * T_)) * (q1_ - q0_);
    }

    Vec7 ApproachMotion::jerk(double t) const
    {
        if (t < 0.0 || t > T_)
            return Vec7::Zero();
        const double x = t / T_;
        return ((60.0 - 360.0 * x + 360.0 * x * x) / (T_ * T_ * T_)) * (q1_ - q0_);
    }

}  // namespace se3pf
