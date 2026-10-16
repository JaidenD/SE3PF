#pragma once
#include "se3pf/types.hpp"
#include "se3pf/experiments.hpp"
#include "se3pf/safety.hpp"

#include <string>
#include <vector>

namespace se3pf {

    // ---------------------------------------------------------------------------------------------------------
    // The preflight (design doc, Hardware): before anything moves, simulate exactly the motion the robot will run and
    // refuse to move if any check of the safety monitor fails at any control period.

    // The rigid-body dynamics of the arm for eq:fb: M(q) qdd + C(q, qd) qd + G(q) = tau. On the robot the driver
    // wraps franka::Model with the end-effector load configured in Desk; in tests and in MuJoCo other models stand in.
    class DynamicsModel {
        public:
            virtual ~DynamicsModel() = default;
            virtual Mat<7, 7> mass(const Vec7& q) const = 0;                // M(q), joint inertia included
            virtual Vec7 coriolis(const Vec7& q, const Vec7& dq) const = 0;  // C(q, qd) qd
            virtual Vec7 gravity(const Vec7& q) const = 0;                   // G(q)
    };

    // What the preflight found. ok is true only if every check passed at every period, the run ended at rest (the
    // end condition of the robot run) and the outputs converged (final |xi| and |beta| at most kConvergence).
    struct PreflightReport {
        static constexpr double kConvergence = 1e-3;

        bool ok = false;
        Violation violation = Violation::None;  // the violation that failed it (None if ok, or if it did not converge)
        std::string reason;                     // empty if ok
        double tFail = -1.0;                    // time of the failure, s (-1 if ok)
        bool torqueChecked = false;             // false without a dynamics model: torque and torque-rate unchecked
        bool converged = false;                 // final |xi| and |beta| <= kConvergence
        double simulated = 0.0;                 // s of motion simulated (to the end condition or the failure)
        double startXi = 0.0, startBeta = 0.0;  // the outputs at the start configuration
        SafetyLimits limits;                    // the monitor's limits it was checked against
        RunPeaks peaks;                         // over every simulated period

        // "OK" or "FAILED: <reason> at t = ..." followed by the peaks and their bounds
        std::string summary() const;
    };

    // Simulates the experiment exactly as the robot driver runs it (ExperimentLoop with setup.controller(): the same
    // gains, reference and mu ramp): at rest at setup.start() at t = 0, one cycle every RunConstants::dt with the
    // acceleration qdd* held over the period (q += dt qd + dt^2 / 2 qdd*, qd += dt qdd*: eq:fb with an exact model),
    // through the end condition (t >= T and |qd| < restSpeed) or the first violation (a Timeout at T + 5 s
    // included). With dyn, every cycle also gets tau = M qdd* + C qd and G from dyn and the torque checks; with
    // dyn == nullptr they are skipped and the report says so (torqueChecked = false). If log is not null it receives
    // one row per simulated cycle (the prediction a robot log can be compared with). If start is not null the motion
    // starts at rest at *start instead of setup.start(): the robot driver checks the run again from the configuration
    // the arm actually reached after the approach. Offline.
    PreflightReport preflight(const ExperimentSetup& setup, const DynamicsModel* dyn, const SafetyMonitor& monitor,
                              std::vector<ExperimentLogRow>* log = nullptr, const Vec7* start = nullptr);

}  // namespace se3pf
