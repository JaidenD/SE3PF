#pragma once
#include "se3pf/types.hpp"
#include "se3pf/controller.hpp"
#include "se3pf/panda.hpp"
#include "se3pf/panda_example.hpp"
#include "se3pf/safety.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace se3pf {

    // ---------------------------------------------------------------------------------------------------------
    // The hardware experiments (docs/HARDWARE.md): what the robot driver runs, as plain data, and the offline setup
    // that turns one into a controller and a start configuration. Every experiment uses the paper's Panda example
    // (PandaExample: l = 0.1 m, the coning circle, the rest-posture lift and beta = X^T (q - q*) of Lemma 3), starts
    // at rest at s0 = 0.1 and ends at rest.

    // A smooth speed profile for a velocity reference (eq:tan-vel): nu_d rises from 0 to nuMax over rampUp
    // (quintic smoothstep), stays at nuMax for cruise, falls back to 0 over rampDown (smoothstep again) and stays 0.
    // nu_d and its first two time derivatives are continuous, so the reference starts and ends at rest.
    struct SpeedProfile {
        double nuMax = 0.0;     // in the units of the metric per second
        double rampUp = 0.0;    // s
        double cruise = 0.0;    // s
        double rampDown = 0.0;  // s

        double end() const { return rampUp + cruise + rampDown; }
        // Arc length travelled by the reference: the smoothstep averages 1/2 over a ramp.
        double distance() const { return nuMax * (0.5 * rampUp + cruise + 0.5 * rampDown); }
        // nu_d(t) and its exact derivative d nu_d / dt; zero before 0 and after end(). Online.
        void eval(double t, double& nu, double& dnuDt) const;
    };

    // One hardware experiment. A run starts at rest at the configuration with pose Psi(s0, xi0) and beta = beta0
    // (s0 = PandaExample::s0 = 0.1), and the reference ends at `duration` (T); the run then lasts until the arm is at
    // rest (|qd| < RunConstants::restSpeed), or fails with a Timeout at T + RunConstants::timeoutExtra.
    struct ExperimentSpec {
        std::string name;     // the command-line name
        std::string title;
        std::string purpose;  // what it shows
        std::string watch;    // what the operator watches, and when the run counts as passed

        Vec5 xi0 = Vec5::Zero();  // start: transverse output (|xi0| = 0.0312 is the paper's start)
        double beta0 = 0.0;       // start: redundancy output, rad

        // Tangential reference: Velocity follows nu_d(t) = speed (eq:tan-vel); Position holds eta_d(t) = s0
        // (eq:tan-pos), the point of the path where the run starts.
        TangentialReference::Mode mode = TangentialReference::Mode::Velocity;
        SpeedProfile speed;

        double duration = 0.0;  // T, s
        Gains gains;
        double rampTime = 1.0;  // T of the mu ramp (Controller::setRampTime), s
        SafetyLimits limits;    // the monitor's values; only tighter than HardwareBounds take effect
    };

    // The experiments in the order the lab runs them: hold, loops, recover, paper, self-motion.
    const std::vector<ExperimentSpec>& experiments();
    // The experiment of that name, or nullptr.
    const ExperimentSpec* findExperiment(const std::string& name);

    // What every run on the robot (and its preflight) shares.
    struct RunConstants {
        static constexpr double dt = 1e-3;            // control period, s (libfranka's 1 kHz)
        static constexpr double restSpeed = 0.02;     // a run is finished once t >= T and |qd| < restSpeed, rad/s
        static constexpr double timeoutExtra = 5.0;   // ... or fails with Timeout at T + timeoutExtra, s
        // The safety stop (SafetyStop) ends once |qd| < stopSpeed has held for stopSettle, or after stopTime.
        static constexpr double stopSpeed = 0.01;     // rad/s
        static constexpr double stopSettle = 0.1;     // s
        static constexpr double stopTime = 2.0;       // s
        // Every joint's stop torque is clamped to stopTorqueFrac of its torque limit (the robot adds G: the delivered
        // torque stays below 60 % + |G_i|, about 85 %, below the robot's own limit).
        static constexpr double stopTorqueFrac = HardwareBounds::tauFrac;
        // The stop's joint stiffness K_s = (300, 300, 300, 300, 50, 50, 20) N m/rad: a weight the robot's G does not
        // know deflects the hold by less than 0.02 rad (the 0.76 kg hand of the MuJoCo Panda missing in G: 0.014 rad
        // on joint 2, 0.018 on joint 6; with 200 on joints 1-4 it was 0.0202).
        static Vec7 stopStiffness();
        // The stop's joint damping D = (30, 30, 30, 30, 3, 3, 1.5) N m s/rad: the earlier implementation's damping
        // stop (10, 10, 10, 10, 3, 3, 1.5) with joints 1-4 raised to damp the hold (damping ratio about 0.8 with
        // MuJoCo's joint inertias; with 10 the hold rang for 1.9 s after a stop at 50 % speed). At |qd_i| = 50 % of
        // the velocity limits it asks for at most 33 N m (joints 1-4) and 3.9 N m (5-7). Every gain of K_s and D is
        // at most 60 % of the same gain of libfranka's joint impedance example, which runs on the Panda (K 600, 600,
        // 600, 600, 250, 150, 50; D 50, 50, 50, 50, 30, 25, 15).
        static Vec7 stopDamping();
    };

    // The safety stop (docs/HARDWARE.md, section 7). On the first violation the driver (and the MuJoCo closed loop)
    // stops following and holds the arm at q_stop, the measured configuration of the period of the violation, with a
    // damped joint hold (the robot adds G):
    //     tau = -K_s (q - q_stop) - D qd,
    // each joint clamped to +-stopTorqueFrac tauMax_i. D dissipates the motion; K_s holds a weight that G does not
    // know, which damping alone cannot (an end effector missing in Desk: with D alone the flange fell 0.3 m in the
    // 2 s of the stop in MuJoCo). On a joint whose q or q_stop is not finite only the damping acts (step latches
    // q_stop for such a joint at its first finite q), on a joint whose qd is not finite only the stiffness (zero
    // torque if neither is finite). The stop is over once |qd| <
    // stopSpeed has held for stopSettle, or after stopTime; then the torque loop ends and the robot's own controller
    // takes the arm. Online: no allocation, no exceptions.
    class SafetyStop {
        public:
            // Starts the stop at time t at the measured configuration q (latched as q_stop).
            void begin(double t, const Vec7& q) noexcept;
            bool active() const noexcept { return active_; }
            double startTime() const noexcept { return t0_; }
            const Vec7& holdConfiguration() const noexcept { return qStop_; }
            // The last step ended the stop because the arm had been at rest for stopSettle (not at stopTime).
            bool settled() const noexcept { return settled_; }

            // One period of the stop at time t: the torque to command; done is set once the stop is over (that
            // period's torque is still sent, with MotionFinished). Zero torque and done if the stop is not active.
            Vec7 step(double t, const Vec7& q, const Vec7& dq, bool& done) noexcept;

            // The stop's torque law at (q, qd) for the hold configuration qStop.
            static Vec7 torque(const Vec7& q, const Vec7& dq, const Vec7& qStop) noexcept;

        private:
            bool active_ = false;
            double t0_ = 0.0;
            double restSince_ = -1.0;  // start of the current stretch with |qd| < stopSpeed; -1: none
            bool settled_ = false;
            Vec7 qStop_ = Vec7::Zero();
    };

    // An experiment made concrete with the Panda example: its start configuration, its tangential reference and its
    // controller. Keeps a reference to the example, which must outlive it (and every controller it builds).
    // Construction is offline: throws std::invalid_argument for an invalid spec (duration not finite and positive,
    // a speed profile longer than the duration or with negative times, a negative ramp time), and std::runtime_error
    // if no start configuration exists.
    class ExperimentSetup {
        public:
            ExperimentSetup(const PandaExample& example, const ExperimentSpec& spec);
            ExperimentSetup(const PandaExample&& example, const ExperimentSpec& spec) = delete;  // no temporaries

            const PandaExample& example() const { return example_; }
            const ExperimentSpec& spec() const { return spec_; }
            static constexpr double s0() { return PandaExample::s0; }
            const Vec7& start() const { return q0_; }   // q0: pose Psi(s0, xi0), beta(q0) = beta0
            double duration() const { return spec_.duration; }
            double timeout() const { return spec_.duration + RunConstants::timeoutExtra; }

            // nu_d(t) of the speed profile (Velocity), or eta_d(t) = s0 (Position)
            TangentialReference reference() const;
            // Controller<PandaArm> with the example's tube and beta, the spec's gains, reference() and the mu ramp,
            // reset to eta-hat = s0. Every user (driver, preflight, MuJoCo) builds its controller here.
            Controller<PandaArm> controller() const;
            // The spec's limits, clamped to HardwareBounds
            SafetyMonitor monitor() const { return SafetyMonitor(spec_.limits); }

        private:
            const PandaExample& example_;
            ExperimentSpec spec_;
            Vec7 q0_;
    };

    // The model terms of eq:fb at one state: M(q), c = C(q, qd) qd and G(q). On the robot they come from
    // franka::Model (with the end-effector load configured in Desk), in the preflight from a DynamicsModel.
    struct ModelTerms {
        Mat<7, 7> M = Mat<7, 7>::Zero();
        Vec7 c = Vec7::Zero();
        Vec7 G = Vec7::Zero();
    };

    // One control cycle of an experiment, shared by the robot driver, the preflight and the MuJoCo closed loop, so
    // that all three run exactly the same controller and the same checks in the same order. Construction is offline;
    // follow and observe are online (no allocation, no exceptions).
    class ExperimentLoop {
        public:
            using Status = Controller<PandaArm>::Status;

            // What one cycle of the following phase decided
            struct Cycle {
                Violation violation = Violation::None;
                bool finished = false;      // t >= T, |qd| < restSpeed, no violation: send tau and end the run
                Vec7 qdd = Vec7::Zero();    // qdd* = A^-1 (mu - Adot qd)
                Vec7 tau = Vec7::Zero();    // commanded torque M qdd* + c (what libfranka is sent; it adds G)
                Vec7 tauG = Vec7::Zero();   // tau + G, the torque the joints deliver
                Status status;
            };

            // The setup must outlive the loop (the controller keeps references into its example).
            ExperimentLoop(const ExperimentSetup& setup, const SafetyMonitor& monitor);

            const ExperimentSetup& setup() const { return setup_; }
            const SafetyMonitor& monitor() const { return monitor_; }
            const Controller<PandaArm>& controller() const { return controller_; }

            // One cycle at time t (t = 0 at the first cycle) from the measured q and qd, and the model terms at
            // (q, qd), or nullptr (kinematic only: tau = tauG = 0 and no torque check). In order: finite q and qd
            // (NonFinite), Controller::step (Projection), checkOutputs, checkState, checkAcceleration, then with a
            // model tau = M qdd* + c and checkTorque(tau + G, tau, tau of the previous cycle (0 before the first),
            // dt), then Timeout (t > T + timeoutExtra), then finished. The first violation is reported; after one,
            // the caller switches to the safety stop (SafetyStop) and calls follow no more.
            void follow(double t, const Vec7& q, const Vec7& dq, const ModelTerms* model, Cycle& out);

            // Outputs at q with eta_dot and beta_dot from qd, for the log of the stopping phase (no reference, no
            // command). False if the projection fails.
            bool observe(const Vec7& q, const Vec7& dq, Status& st) const;

        private:
            const ExperimentSetup& setup_;
            SafetyMonitor monitor_;
            Controller<PandaArm> controller_;
            Vec7 tauPrev_ = Vec7::Zero();
            double eta_ = PandaExample::s0;  // last projected eta, the start of observe's projection
    };

    // ---------------------------------------------------------------------------------------------------------
    // Logs and their peaks

    // One row of a run's log, at one control cycle: the robot driver writes one per callback, the preflight and the
    // MuJoCo closed loop one per period (CSV columns: apps/robot/robot_log.hpp).
    struct ExperimentLogRow {
        double t = 0.0;
        int phase = 0;              // 0 following, 1 safety stop (damped hold, SafetyStop)
        double eta = 0.0, etaDot = 0.0;
        double ref = 0.0;           // nu_d (velocity reference) or eta_d (position reference)
        double eTan = 0.0;          // e_v = eta_dot - nu_d, or e_eta = eta - eta_d wrapped
        Vec5 xi = Vec5::Zero();
        double xiNorm = 0.0;
        double beta = 0.0, betaDot = 0.0;
        double sigmaMin = 0.0;      // sigma_min(A)
        double ramp = 0.0;          // the mu ramp factor r(t / T)
        int newton = 0;             // Newton updates of the projection
        Vec7 q = Vec7::Zero(), dq = Vec7::Zero();
        Vec7 qdd = Vec7::Zero();    // qdd* (zero in the stopping phase)
        Vec7 tau = Vec7::Zero();    // commanded torque (M qdd* + c, or the stop's -K_s (q - q_stop) - D qd)
        Vec7 tauG = Vec7::Zero();   // tau + G
        double computeUs = 0.0;     // driver only: time spent in the callback, microseconds
        double periodMs = 0.0;      // driver only: libfranka's period of this callback, ms
        double successRate = 0.0;   // driver only: RobotState::control_command_success_rate
    };

    // A row from one following cycle (phase 0) or from the outputs observed in the stopping phase (phase 1).
    ExperimentLogRow makeLogRow(double t, int phase, const ExperimentLoop::Status& st, const Vec7& q, const Vec7& dq,
                                const Vec7& qdd, const Vec7& tau, const Vec7& tauG);

    // The peaks of a run over its following rows (phase 0), as the preflight reports them, measured the way the
    // SafetyMonitor checks them. Joints are zero-based. Offline (add evaluates the kinematics for the flange speed).
    struct RunPeaks {
        bool torques = false;            // set if the rows carry torques (tau, tauG); otherwise the torque peaks stay 0
        long rows = 0;                   // following rows seen
        double minJointMargin = 0.0;     // rad
        int minMarginJoint = -1;
        double maxDqFrac = 0.0;          // max_i |qd_i| / dqMax_i
        int maxDqJoint = -1;
        double maxDdqFrac = 0.0;         // max_i |qdd*_i| / ddqMax_i
        int maxDdqJoint = -1;
        double maxTauFrac = 0.0;         // max_i |tau_i + G_i| / tauMax_i
        int maxTauJoint = -1;
        double maxTauRate = 0.0;         // max_i |tau_i - tau_i,prev| / dt of the commanded torque, N m/s
        double maxFlangeSpeed = 0.0;     // m/s
        double maxXi = 0.0;              // max |xi|
        double maxBeta = 0.0;            // max |beta|, rad
        double minSigma = 0.0;           // min sigma_min(A)
        int maxNewton = 0;
        double finalT = 0.0, finalXi = 0.0, finalBeta = 0.0, finalEtaDot = 0.0, finalETan = 0.0;  // the last row

        // Adds one row; rows of phase 1 are ignored. The torque rate compares with the previous row's tau (zero
        // before the first row), over RunConstants::dt.
        void add(const ExperimentLogRow& row);

        // One line: the peaks with the bound of each in brackets
        std::string summary(const SafetyLimits& limits) const;

        private:
            Vec7 tauPrev_ = Vec7::Zero();
    };

}  // namespace se3pf
