#pragma once
// The hardware experiments (se3pf/experiments.hpp, docs/HARDWARE.md) on the MuJoCo Menagerie Panda: the best
// prediction of the robot before it moves. Part of se3pf_mujoco (built with -DSE3PF_WITH_MUJOCO=ON), never of se3pf.
//   1. MujocoDynamics: MuJoCo's M (joint armature included), C(q, qd) qd and G as a DynamicsModel, so the preflight
//      can check every experiment with realistic torques (the robot uses franka::Model instead).
//   2. simulateExperimentMujoco: an experiment closed loop at the torque level, run as the robot driver runs it
//      (ExperimentLoop::follow with the model terms, the safety stop on a violation, the same end conditions), with
//      the torque command held for each 1 ms period, libfranka's low-pass filter and rate limiter on the command,
//      the robot's G added to it, and the plant's joint damping (1 N m s/rad) left uncompensated, as run C of
//      panda_mujoco (docs/MUJOCO.md).
#include "mujoco_panda.hpp"
#include "se3pf/experiments.hpp"
#include "se3pf/preflight.hpp"
#include "se3pf/safety.hpp"

#include <memory>
#include <string>
#include <vector>

namespace se3pf {

    // MuJoCo's dynamics of one model file as a DynamicsModel: M(q) from mj_fullM (the 7 x 7 arm block, armature
    // included), C(q, qd) qd = qfrc_bias - G and G(q) = qfrc_bias at rest. For panda.xml the arm block of M ignores
    // the finger dofs' coupling (the fingers counted as rigid). The passive forces (joint damping) are not part of it.
    // Keeps its own MujocoPanda, which the const functions evaluate (MuJoCo needs a scratch state), so it is not
    // thread safe. Throws std::runtime_error if the model cannot be loaded.
    class MujocoDynamics final : public DynamicsModel {
        public:
            explicit MujocoDynamics(const std::string& xmlPath);
            Mat<7, 7> mass(const Vec7& q) const override;
            Vec7 coriolis(const Vec7& q, const Vec7& dq) const override;
            Vec7 gravity(const Vec7& q) const override;

        private:
            std::unique_ptr<MujocoPanda> arm_;
    };

    struct MujocoExperimentOptions {
        std::string plant;           // the plant's model file (panda_nohand.xml, or panda.xml: the hand as a load)
        std::string model;           // model file of M and C qd in the torque; empty: the plant's own
        bool plantGravity = true;    // G added to the command is the plant's (as the robot adds G with the load
                                     // configured in Desk); false: the model's (a load missing in Desk)
        bool damping = true;         // joint damping 1 N m s/rad in the plant (the robot's friction), uncompensated
        bool filter = true;          // libfranka's low-pass filter (100 Hz) and rate limiter (1000 N m/s) on the command
        int substeps = 1;            // plant steps per 1 ms period (RK4)
        int pushJoint = -1;          // zero-based joint of an external torque push (-1: none) ...
        double pushTorque = 0.0;     // ... its torque, N m ...
        double pushStart = 0.0;      // ... from this time, s ...
        double pushDuration = 0.0;   // ... for this long, s
    };

    struct MujocoExperimentRun {
        std::vector<ExperimentLogRow> log;     // one row per period, as the driver's log
        Violation violation = Violation::None; // the safety stop's reason (None if the run finished)
        double tStop = -1.0;                   // time of the violation (the safety stop starts), s
        double maxStopDrift = 0.0;             // max_i |q_i - q_stop,i| over the safety stop, rad
        bool stopSettled = false;              // the safety stop ended at rest (|qd| < stopSpeed for stopSettle)
        bool finished = false;                 // ended at rest after the reference, no violation
        double tEnd = 0.0;                     // last period, s
        double endSpeed = 0.0;                 // |qd| at the last period, rad/s
        RunPeaks peaks;                        // over the following periods (torques included)
        double steadyXi = 0.0, steadyBeta = 0.0;  // max |xi| and |beta| over the last 5 s of following
        double cruiseEtaDot = 0.0;             // mean eta_dot over the cruise of a velocity reference (0 otherwise)
        int maxContacts = 0;
        long clampedPeriods = 0;               // periods with a torque beyond the actuator range
        int warnings = 0;                      // MuJoCo warnings
    };

    // Runs the experiment on the MuJoCo plant: at rest at setup.start(), every 1 ms the driver's cycle with the
    // model's M, C qd and G at the plant's state, the command (filtered and rate limited as libfranka does) plus G
    // held for the period, until the run is finished, or until the safety stop after a violation ends (SafetyStop:
    // at rest for stopSettle, or stopTime), or a MuJoCo warning. Offline. Throws std::invalid_argument for bad options and
    // std::runtime_error if a model cannot be loaded.
    MujocoExperimentRun simulateExperimentMujoco(const ExperimentSetup& setup, const SafetyMonitor& monitor,
                                                 const MujocoExperimentOptions& options);

    // The safety stop alone (SafetyStop), from any state: what the arm does after a violation at (q0, qd0).
    struct MujocoStopRun {
        double duration = 0.0;       // until the stop ended (settled, or stopTime), s
        bool settled = false;        // at rest (|qd| < stopSpeed) for stopSettle before stopTime
        double endSpeed = 0.0;       // |qd| at the end, rad/s
        Vec7 drift = Vec7::Zero();   // max |q_i - q0_i| over the stop, per joint, rad
        double maxTorque = 0.0;      // max_i |tau_i| commanded by the stop (before libfranka's filter), N m
        double minJointMargin = 0.0; // the closest any joint came to a position limit, rad
        int warnings = 0;            // MuJoCo warnings
    };

    // Runs SafetyStop on the plant from (q0, qd0) as the driver runs it after a violation: the stop's torque through
    // libfranka's filter and rate limiter (the last command before the stop: zero), plus the robot's G, held for each
    // 1 ms period, with the plant's joint damping (1 N m s/rad). gravityModel: the model file whose G the robot adds
    // (empty: the plant's own; panda_nohand.xml under a panda.xml plant: the hand missing in Desk). Offline.
    MujocoStopRun simulateStopMujoco(const std::string& plant, const std::string& gravityModel, const Vec7& q0,
                                     const Vec7& qd0);

}  // namespace se3pf
