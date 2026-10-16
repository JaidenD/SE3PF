#pragma once
#include "se3pf/types.hpp"
#include "se3pf/panda.hpp"

namespace se3pf {

    // ---------------------------------------------------------------------------------------------------------
    // The online safety monitor of the robot driver (design doc, Hardware; docs/HARDWARE.md)
    //
    // Nothing moves on the robot unless the preflight (preflight.hpp) has simulated the motion with every check below
    // at every control period, and the driver runs the same checks at every period of the 1 kHz torque loop: on the
    // first violation it switches to joint damping and brings the arm to rest. The checks are fixed-size, allocation
    // free and noexcept, so the 1 kHz callback may call them.

    // Why a run stops (or a preflight fails). Every violation stops the arm by joint damping on the robot.
    enum class Violation {
        None,
        Projection,     // Controller::step failed: the pose left the tube, A nearly singular, or a non-finite result
        JointLimit,     // a joint closer than jointMargin to one of its position limits
        JointVelocity,  // |qd_i| above dqFrac of joint i's velocity limit
        Acceleration,   // commanded |qdd*_i| above ddqFrac of joint i's acceleration limit
        Torque,         // delivered |tau_i + G_i| above tauFrac of joint i's torque limit
        TorqueRate,     // commanded |tau_i - tau_i,prev| / dt above tauRateFrac of 1000 N m/s
        FlangeSpeed,    // translational speed of the flange above flangeSpeed
        Tube,           // |xi| above xi (the tube of the outputs)
        Redundancy,     // |beta| above beta
        Conditioning,   // sigma_min(A) below sigmaMin (Proposition 2: A invertible near the lift)
        NonFinite,      // a measured or computed value is NaN or infinite (or a period dt <= 0)
        Timeout         // the arm was not at rest within the experiment's duration + 5 s
    };

    // A short English description, for the driver's messages and the preflight report.
    const char* toString(Violation v) noexcept;

    // The values a SafetyMonitor checks against: the HardwareBounds (panda.hpp) as runtime values. An experiment may
    // ask for tighter values; SafetyMonitor clamps every value to its HardwareBounds bound, so a monitor can never be
    // looser than the bounds, whatever an experiment asks for.
    struct SafetyLimits {
        double jointMargin = HardwareBounds::jointMargin;  // rad from every position limit, at least
        double dqFrac = HardwareBounds::dqFrac;            // |qd_i| <= dqFrac * dqMax_i
        double ddqFrac = HardwareBounds::ddqFrac;          // commanded |qdd*_i| <= ddqFrac * ddqMax_i
        double tauFrac = HardwareBounds::tauFrac;          // delivered |tau_i + G_i| <= tauFrac * tauMax_i
        double tauRateFrac = HardwareBounds::tauRateFrac;  // commanded torque rate <= tauRateFrac * tauRateMax
        double flangeSpeed = HardwareBounds::flangeSpeed;  // m/s, translational
        double xi = HardwareBounds::xi;                    // |xi|, in the units of the metric
        double beta = HardwareBounds::beta;                // |beta|, rad
        double sigmaMin = HardwareBounds::sigmaMin;        // sigma_min(A), at least
    };

    // Checks one control period against SafetyLimits. Each check returns the first violation it finds (None if
    // there is none). A bound is kept when the value equals it (|qd_i| = 50 % passes) and violated beyond it; every
    // comparison is written so that a NaN fails. Each joint is checked against its own limit (PandaArm::dqMax, ...).
    // Construction is offline (it may throw); the checks are online: no allocation, no exceptions.
    class SafetyMonitor {
        public:
            // The HardwareBounds themselves.
            SafetyMonitor();

            // Tighter values for one experiment, each clamped to its HardwareBounds bound: a jointMargin or sigmaMin
            // below the bound is raised to it, a fraction, speed, xi or beta above the bound is lowered to it, and
            // tighter values are kept. clamped() says whether any value had to be clamped. Throws
            // std::invalid_argument for a value that is NaN, infinite or negative.
            explicit SafetyMonitor(const SafetyLimits& requested);

            const SafetyLimits& limits() const noexcept { return limits_; }
            bool clamped() const noexcept { return clamped_; }

            // NonFinite, then Tube (|xi| > xi), Redundancy (|beta| > beta), Conditioning (sigmaMin(A) < sigmaMin).
            Violation checkOutputs(const Vec5& xi, double beta, double sigmaMin) const noexcept;

            // NonFinite, then JointLimit (any joint closer than jointMargin to a limit), JointVelocity (any |qd_i| >
            // dqFrac dqMax_i), FlangeSpeed (|v| of the flange twist J_b(q) qd above flangeSpeed).
            Violation checkState(const Vec7& q, const Vec7& dq) const noexcept;

            // NonFinite, then Acceleration (any commanded |qdd*_i| > ddqFrac ddqMax_i).
            Violation checkAcceleration(const Vec7& qdd) const noexcept;

            // tauDelivered = tau_cmd + G is the torque the joints deliver (libfranka is sent tau_cmd and adds the
            // robot's G); tauCmd and tauCmdPrev are the commanded torques of this period and the previous one (zero
            // before the first period: the arm rests under the robot's gravity compensation). NonFinite (also for
            // dt <= 0), then Torque (any |tauDelivered_i| > tauFrac tauMax_i), TorqueRate (any |tauCmd_i -
            // tauCmdPrev_i| / dt > tauRateFrac tauRateMax).
            Violation checkTorque(const Vec7& tauDelivered, const Vec7& tauCmd, const Vec7& tauCmdPrev,
                                  double dt) const noexcept;

        private:
            SafetyLimits limits_;
            bool clamped_ = false;
            // Per-joint limits, computed once
            Vec7 qLow_, qHigh_;    // qMin + jointMargin, qMax - jointMargin
            Vec7 dqLimit_;         // dqFrac * dqMax
            Vec7 ddqLimit_;        // ddqFrac * ddqMax
            Vec7 tauLimit_;        // tauFrac * tauMax
            double tauRateLimit_;  // tauRateFrac * tauRateMax, N m/s
    };

    // ---------------------------------------------------------------------------------------------------------
    // The approach to a start configuration (design doc, Hardware: approach)

    // A minimum-jerk joint motion q(t) = q0 + r(t / T) (q1 - q0), with r the quintic smoothstep (quinticRamp,
    // controller.hpp): velocity, acceleration and jerk start and end at rest except the jerk, and the peaks of r', r''
    // and r''' are 15/8, 10/sqrt(3) and 60. T is the shortest time that keeps every joint within kVelocityFraction,
    // kAccelerationFraction and kJerkFraction of its own limits, and at least kMinTime. The driver commands it with
    // franka::JointPositions from the robot's q_d at the first callback. Construction is offline (throws
    // std::invalid_argument for a non-finite q0 or q1); the evaluation is online.
    class ApproachMotion {
        public:
            static constexpr double kVelocityFraction = 0.20;      // of PandaArm::dqMax
            static constexpr double kAccelerationFraction = 0.10;  // of PandaArm::ddqMax
            static constexpr double kJerkFraction = 0.05;          // of PandaArm::jerkMax
            static constexpr double kMinTime = 3.0;                // s

            ApproachMotion(const Vec7& q0, const Vec7& q1);

            double duration() const { return T_; }
            const Vec7& from() const { return q0_; }
            const Vec7& to() const { return q1_; }
            // max_i |q1_i - q0_i|, rad
            double maxJointDistance() const { return (q1_ - q0_).cwiseAbs().maxCoeff(); }

            // At time t (clamped to [0, T]): q1 exactly for t >= T.
            Vec7 position(double t) const;
            Vec7 velocity(double t) const;
            Vec7 acceleration(double t) const;
            Vec7 jerk(double t) const;  // zero outside [0, T]

        private:
            Vec7 q0_, q1_;
            double T_ = kMinTime;
    };

}  // namespace se3pf
