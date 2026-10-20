#pragma once
#include "se3pf/types.hpp"

#include <vector>

namespace se3pf {

    // ---------------------------------------------------------------------------------------------------------
    // Joint identification (CTRL-1): the breakaway torque and the command-side inertia of each joint, measured on
    // the robot with torque pulses, one joint at a time, while the damped hold (SafetyStop::torque) keeps the
    // others still. Not path following: the robot driver runs it as its own command, with the same checks, the same
    // safety stop and the same logging as an experiment. Why: the first session (2026-10-19) found breakaway
    // torques two to four times Gaz et al. 2019 and a wrist that franka::Model's link-side inertia makes the
    // controller treat as 30 times lighter than it is; both are plant-model terms of eq:fb, and neither was measured.
    //
    // For every joint in the list, in each direction s = +1, -1:
    //   ramp   tau_j = s r (t - t0) from zero at rampRate until the joint moves (|dq_j| >= motionSpeed for
    //          motionHold), or until torqueCap_j: the breakaway torque (command side; the driver logs the joint
    //          torque sensor next to it);
    //   step   tau_j = s (breakaway + stepExtra[k]) for k = 0, 1, held for stepDuration or until the joint has
    //          travelled travelLimit or reached speedLimit; two amplitudes, so that the inertia follows from their
    //          difference without knowing the kinetic friction (tau2 - tau1 = M (a2 - a1));
    // with settle seconds of the hold between pulses and the hold re-latched at the joint's new position after each
    // pulse. Every change of the command is blended over blend seconds (the torque-rate margin). The other joints
    // are held at the start configuration throughout. Online: no allocation, no exceptions.

    struct IdentifySettings {
        double rampRate = 0.5;        // N m/s
        Vec7 torqueCap = (Vec7() << 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5).finished();  // N m, the ramp's ceiling
        // N m above the breakaway, the two steps, per joint: the big joints (1 kg m^2 and more) need 1 to 2 N m for
        // a measurable acceleration, the wrist (0.1 kg m^2 command side, identify_1 of 2026-10-20) 0.1 to 0.2 N m
        // or it hits the speed limit within 0.1 s
        Vec7 stepExtra1 = (Vec7() << 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1).finished();
        Vec7 stepExtra2 = (Vec7() << 2.0, 2.0, 2.0, 2.0, 0.2, 0.2, 0.2).finished();
        double stepDuration = 0.3;    // s
        double travelLimit = 0.08;    // rad from the pulse's start: the pulse ends there
        double speedLimit = 0.5;      // rad/s: the pulse ends there
        double motionSpeed = 0.02;    // rad/s (identify_1: the hold's settling on the gravity-loaded joints 2 to 4
                                      // tripped 0.01 within 0.1 to 0.6 s of the ramp) ...
        double motionHold = 0.03;     // ... for this long: the joint has broken away
        double settle = 2.0;          // s of hold between pulses (the hold must be still before a ramp)
        double blend = 0.05;          // s over which a change of the command is blended
    };

    enum class IdentifyPhase { Hold = 0, Ramp = 1, Step = 2 };

    // One control cycle of the identification run, for the log (the driver's CSV: apps/robot/robot_log.hpp)
    struct IdentifyLogRow {
        double t = 0.0;
        int phase = 0;           // IdentifyPhase, or 3 for the safety stop
        int joint = 0;           // 1-based active joint (0: none)
        int direction = 0;       // +1, -1, 0
        int pulse = 0;           // 0 ramp, 1 step 1, 2 step 2 (-1 none)
        double pulseTorque = 0;  // the active joint's commanded pulse torque, N m
        Vec7 tau = Vec7::Zero();   // commanded torque (the robot adds G)
        Vec7 q = Vec7::Zero(), dq = Vec7::Zero();
        Vec7 tauJ = Vec7::Zero();  // the measured link-side torque (driver only)
        Vec7 G = Vec7::Zero();     // the model's gravity torque
        double computeUs = 0.0, periodMs = 0.0;
    };

    class JointIdentification {
        public:
            // joints: zero-based indices, in the order to identify; q0: the hold configuration. Offline; throws
            // std::invalid_argument for an empty or out-of-range list or invalid settings.
            JointIdentification(const IdentifySettings& settings, const std::vector<int>& joints, const Vec7& q0);

            // One cycle at time t: the torque to command. done is set once the sequence is over (after the last
            // settle). Online.
            Vec7 step(double t, const Vec7& q, const Vec7& dq, bool& done) noexcept;

            IdentifyPhase phase() const noexcept { return phase_; }
            int activeJoint() const noexcept { return joint_; }      // zero-based, -1 if none
            int direction() const noexcept { return dir_; }
            int pulse() const noexcept { return pulseKind_; }        // 0 ramp, 1, 2 steps, -1 none
            double pulseTorque() const noexcept { return pulseTau_; }
            // The breakaway found so far for (joint, direction +1 -> 0, -1 -> 1), NaN if not found (ramp to cap)
            double breakaway(int joint, int dirIndex) const noexcept { return breakaway_[joint][dirIndex]; }
            double duration() const noexcept;  // worst-case length of the whole sequence, s

        private:
            struct Pulse {
                int joint, dir, kind;  // kind 0 ramp, 1, 2 steps
            };
            Vec7 holdTorque(const Vec7& q, const Vec7& dq) const noexcept;
            void beginSettle(double t, const Vec7& q) noexcept;
            void beginPulse(double t, const Vec7& q) noexcept;

            IdentifySettings s_;
            std::vector<Pulse> pulses_;
            std::size_t next_ = 0;
            Vec7 qHold_;
            IdentifyPhase phase_ = IdentifyPhase::Hold;
            int joint_ = -1, dir_ = 0, pulseKind_ = -1;
            double pulseTau_ = 0.0;
            double t0_ = 0.0;             // start of the current phase
            double q0_ = 0.0;             // the active joint's position at the pulse's start
            double movingSince_ = -1.0;   // start of the current stretch with |dq_j| >= motionSpeed
            double breakaway_[7][2];
            double stepTau_ = 0.0;        // the step's amplitude
            Vec7 lastTau_ = Vec7::Zero();  // the command of the last cycle, for the blend
            Vec7 blendFrom_ = Vec7::Zero();
            double blendT0_ = -1.0;
            bool done_ = false;
    };

    // What the identification found, fitted from the log (offline)
    struct JointIdentificationResult {
        bool attempted[2] = {false, false};  // per direction (+, -): a ramp was run
        bool found[2] = {false, false};      // the ramp moved the joint
        double breakawayRaw[2] = {0, 0};     // the commanded torque when the joint reached motionSpeed, N m
        // The breakaway corrected for the ramp's lag: after the breakaway at tau_b the net torque grows at the
        // ramp rate r, so a joint of inertia M reaches motionSpeed v_m only sqrt(2 M v_m / r) later, when the ramp
        // is at tau_b + sqrt(2 M v_m r) (0.1 to 0.2 N m for this arm's joints); M from the steps, else uncorrected
        double breakawayCmd[2] = {0, 0};     // N m
        double breakawaySensor[2] = {0, 0};  // tau_J - G at that moment, N m (0 without sensor data)
        bool inertiaFound[2] = {false, false};
        double inertiaCmd[2] = {0, 0};       // M from the two steps, command side: (tau2 - tau1) / (a2 - a1), kg m^2
        double inertiaSensor[2] = {0, 0};    // the same with tau_J - G: the link-side inertia the model should match
        double kineticFriction[2] = {0, 0};  // tau1 - M a1, command side, N m
        double accel[2][2] = {{0, 0}, {0, 0}};  // the mean accelerations of the two steps
    };

    // Fits the results of every joint from the rows of a run; joints that were not identified stay all-false.
    std::vector<JointIdentificationResult> fitIdentification(const std::vector<IdentifyLogRow>& rows,
                                                             const IdentifySettings& settings);

    // One line per joint and direction, for the driver's summary
    std::string identificationSummary(const std::vector<JointIdentificationResult>& r);

}  // namespace se3pf
