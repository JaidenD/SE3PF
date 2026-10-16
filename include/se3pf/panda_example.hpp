#pragma once
#include "se3pf/types.hpp"
#include "se3pf/controller.hpp"
#include "se3pf/lift.hpp"
#include "se3pf/metric.hpp"
#include "se3pf/panda.hpp"
#include "se3pf/path.hpp"
#include "se3pf/redundancy.hpp"
#include "se3pf/tube.hpp"

#include <optional>
#include <vector>

namespace se3pf {

    // The paper's Panda example (design doc, "Panda simulation" and Appendix: Panda), built once and owned in one
    // place: the metric, the coning-circle path, its tube, the rest-posture lift, the kernel frame, the redundancy
    // output beta of Lemma 3 (and the joint coupling q_3 = q_1 for comparison runs), and the two start
    // configurations. The members are declared, and so built, in dependency order: the tube keeps a reference to
    // the path and beta to the lift and the kernel frame, all members of this object. So it cannot be copied or
    // moved (a copy's tube would still point at the original's path). Build it once and keep it for the whole run.
    // Construction is offline: it throws std::invalid_argument for l <= 0 and std::runtime_error if the lift or a
    // start configuration cannot be computed. In Release it takes well under a second.
    class PandaExample {
        public:
            // The two closed loops of the design doc (Panda simulation), both with nu_d = 0.08 and a 1 ms period:
            //   Paper:      at rest at Psi(s0, xi0) with beta = 0, for 20 s (the paper's Panda results);
            //   SelfMotion: on the path at gamma(s0) with beta = 0.2, moving along the path at nu_d, for one loop
            //               (checks the redundancy channel; not in the paper).
            enum class RunKind { Paper, SelfMotion };

            // l = sqrt(w_R / w_p), the metres one radian counts for: W = diag(l^2 I, I). The paper uses l = 0.1 m.
            explicit PandaExample(double l = 0.1);
            PandaExample(const PandaExample&) = delete;
            PandaExample& operator=(const PandaExample&) = delete;

            // The path: coning circle around p0 with radius r and cone angle alpha, R0 = flange attitude at q_r
            static constexpr double radius = 0.12;       // r, m
            static constexpr double cone = 0.25;         // alpha, rad
            static Vector3d center() { return Vector3d(0.40, 0.0, 0.45); }  // p0, m

            // Starts: at rest at s0, at Psi(s0, xi0) with beta = 0 (the paper), or on the path with beta = 0.2
            static constexpr double s0 = 0.1;
            // Simulation only: above the hardware bound |beta| <= 0.05 rad (HardwareBounds::beta, panda.hpp); a
            // self-motion demo on the robot must start at |beta| <= 0.04.
            static constexpr double selfMotionBeta0 = 0.2;
            static Vec5 xi0() { return (Vec5() << 0.02, -0.015, 0.01, -0.005, 0.015).finished(); }

            // Controller and simulation
            static constexpr double nuD = 0.08;          // velocity reference nu_d (about 9.6 s per loop)
            static constexpr double lambda = 4.0;        // transverse poles at -lambda: K_p = lambda^2, K_d = 2 lambda
            static constexpr double kt = 5.0;            // tangential velocity gain
            static constexpr double KpBeta = 16.0;       // redundancy gains (eq:beta-law)
            static constexpr double KdBeta = 8.0;
            static constexpr double dt = 1e-3;           // control period, s
            static constexpr double T = 20.0;            // duration, s
            // The longest run simulate accepts, in control steps (stepCount): 2.8 h at 1 kHz, a log of about 3 GB.
            static constexpr long maxSteps = 10'000'000;

            // The parts. Not on a temporary example (deleted &&-overloads): a controller or a redundancy output
            // keeps a reference to them, which would dangle at the end of the full expression.
            const Metric& metric() const& { return metric_; }
            const Matrix3d& R0() const { return R0_; }
            const ConingCircle& path() const& { return path_; }
            const Tube& tube() const& { return tube_; }
            const Lift& lift() const& { return lift_; }
            const KernelFrame& kernelFrame() const& { return kernel_; }
            const LiftRedundancy& redundancy() const& { return beta_; }      // beta = X^T (q - q*), the paper's
            const CouplingRedundancy& coupling() const& { return coupling_; }  // beta = q_3 - q_1
            const Metric& metric() const&& = delete;
            const ConingCircle& path() const&& = delete;
            const Tube& tube() const&& = delete;
            const Lift& lift() const&& = delete;
            const KernelFrame& kernelFrame() const&& = delete;
            const LiftRedundancy& redundancy() const&& = delete;
            const CouplingRedundancy& coupling() const&& = delete;

            // The paper's start: pose Psi(s0, xi0), beta = 0
            const Vec7& paperStart() const { return paperStart_; }
            // The self-motion run's start: pose gamma(s0) (xi0 = 0), beta = selfMotionBeta0
            const Vec7& selfMotionStart() const { return selfMotionStart_; }

            // Any other start: pose Psi(s0, xi) with beta(q) = beta0 for the given redundancy output, by Newton's
            // method from q*(s0) (startConfiguration). Throws std::runtime_error if it fails.
            Vec7 start(const RedundancyOutput& beta, const Vec5& xi, double beta0) const;

            // Gains of both closed loops: transverse and redundancy poles at -lambda (K_p = K_p,beta = 16,
            // K_d = K_d,beta = 8) and k_t = 5.
            static Gains gains();

            // The velocity reference nu_d(t, eta) = nuD, a constant speed along the path (eq:tan-vel)
            static TangentialReference reference();

            // y0 (1 + lambda t) e^{-lambda t}: the designed response of an output that starts at y0 with zero rate
            // and obeys y'' = -lambda^2 y - 2 lambda y' (eq:transverse, eq:beta-law). Each xi_i follows it, so |xi|
            // does too: |xi(t)| = |xi(0)| (1 + 4t) e^{-4t}; and beta(t) = beta(0) (1 + 4t) e^{-4t}.
            static double designedResponse(double y0, double t);

            // How long a run lasts unless told otherwise: T = 20 s for the paper run, one loop L / nu_d (9.626 s)
            // for the self-motion run.
            double duration(RunKind which) const;

        private:
            Metric metric_;
            Matrix3d R0_;
            ConingCircle path_;
            Tube tube_;
            Lift lift_;
            KernelFrame kernel_;  // X of Lemma 3 (tube_ holds the normal frame N)
            LiftRedundancy beta_;
            CouplingRedundancy coupling_;
            Vec7 paperStart_;
            Vec7 selfMotionStart_;
    };

    // One row of a Panda closed-loop log, at a control instant t_k = k dt (before the arm moves on).
    struct PandaLogRow {
        double t = 0.0;
        double eta = 0.0, etaDot = 0.0;
        double ref = 0.0;            // nu_d
        double eTan = 0.0;           // e_v = eta_dot - nu_d
        Vec5 xi = Vec5::Zero();
        double xiNorm = 0.0;
        double xiDesign = 0.0;       // |xi0| (1 + 4t) e^{-4t}: |xi0| = 0.03122 (paper), 0 (self-motion)
        double beta = 0.0;
        double betaDesign = 0.0;     // beta0 (1 + 4t) e^{-4t}: 0 (paper), beta0 = 0.2 (self-motion)
        Vec7 q = Vec7::Zero();       // joint angles, rad
        Vec7 qd = Vec7::Zero();      // joint velocities, rad/s
        Vec7 qdd = Vec7::Zero();     // qdd* = A^-1 (mu - Adot qdot), held over [t, t + dt), rad/s^2
        double sigmaMin = 0.0;       // smallest singular value of A(q)
        double jointMargin = 0.0;    // distance of q to the nearest joint limit, rad
        int iterations = 0;          // Newton updates of the projection at q
    };

    // The log and the summary numbers of one Panda closed loop (the Panda rows of the Closed loops table).
    // Maxima and minima are over the log rows; "joint" fields are zero-based (joint 4 of the Panda is 3).
    struct PandaRun {
        PandaExample::RunKind kind = PandaExample::RunKind::Paper;
        double T = 0.0;                  // duration asked for, s
        double dt = 0.0;                 // control period, s
        std::vector<PandaLogRow> log;

        // Start: the pose Pi(q0) against gamma(s0), and the outputs at q0
        double startAttitudeDeg = 0.0;   // angle between the two attitudes, degrees
        double startPosition = 0.0;      // distance between the two positions, m
        double startXiNorm = 0.0;        // |xi| at t = 0
        double startBeta = 0.0;          // beta at t = 0

        // Outputs against their designed responses
        double maxXiNorm = 0.0;                                      // max |xi|
        double maxXiDeviation = 0.0, maxXiDeviationTime = 0.0;      // max ||xi| - xiDesign|
        double maxBeta = 0.0, maxBetaTime = 0.0;                     // max |beta|
        double maxBetaDeviation = 0.0, maxBetaDeviationTime = 0.0;   // max |beta - betaDesign|
        double finalXiNorm = 0.0;        // at the last row
        double finalEtaDot = 0.0;
        double finalBeta = 0.0;

        // The arm along the run
        double minJointMargin = 0.0, minMarginTime = 0.0;            // smallest distance to a joint limit, rad
        int minMarginJoint = -1;
        double minSigma = 0.0, minSigmaTime = 0.0;                   // smallest sigma_min(A) (Proposition 2)
        double maxVelocityFraction = 0.0, maxVelocityTime = 0.0;     // largest |qd_i| / qd_max,i
        int maxVelocityJoint = -1;
        double maxAccelerationFraction = 0.0, maxAccelerationTime = 0.0;  // largest |qdd*_i| / qdd_max,i
        int maxAccelerationJoint = -1;
        int maxIterations = 0;           // most Newton updates in one cycle (projection at q only)

        // Distance from the lift at the last row, q - q*(eta): back on the lift after the self-motion run
        double finalLiftGap = 0.0;       // max_i |q_i - q*_i(eta)|
        double finalLiftGapNorm = 0.0;   // |q - q*(eta)|

        // false if a step failed (the pose left the tube, or A became singular); the run stops there and the log
        // holds the rows before it. ok covers the controller only: the model qdd = qdd* has no joint limits, so a
        // run can leave the joint range and still be ok. minJointMargin < 0 says it did; on the robot, the driver's
        // safety checks stop the arm at its margins instead.
        bool ok = false;
        double failureTime = -1.0;       // time of the failed step, or -1
    };

    // The Panda closed loop at the joint-acceleration level (design doc, Panda simulation): with an exact model,
    // eq:fb makes qdd = qdd*, so no dynamic model is needed. Controller<PandaArm> with the example's tube, lift
    // redundancy output, gains() and reference(). Start (RunKind):
    //   Paper:      q0 = paperStart(), qd0 = 0 (at rest);
    //   SelfMotion: q0 = selfMotionStart(), qd0 from A(q0) qd0 = (nu_d, 0, ..., 0), so that eta_dot = nu_d and
    //               xi_dot = beta_dot = 0 at t = 0, as the designed response beta0 (1 + 4t) e^{-4t} assumes.
    // Then reset(s0), and at every t_k = k dt, k = 0, ..., floor(T / dt): step gives qdd*, and the arm moves with
    // qdd* held over the period, q <- q + dt qd + dt^2 / 2 qdd*, qd <- qd + dt qdd* (exact for a held
    // acceleration). Stops at the first failed step (ok = false). T defaults to example.duration(which): T = 20 s
    // gives the 20001 rows t = 0, ..., 20 s (floor allows for rounding, so t = T is included when T is a multiple
    // of dt); one loop, L / nu_d = 9.6261 s, gives the 9627 rows t = 0, ..., 9.626 s.
    // Offline: the log is reserved for the whole run before the loop, which neither allocates nor throws. Throws
    // std::invalid_argument unless T >= 0 and dt > 0 are finite and T / dt <= PandaExample::maxSteps (stepCount),
    // and std::runtime_error if the self-motion run's initial velocity cannot be computed.
    PandaRun simulate(const PandaExample& example, PandaExample::RunKind which, std::optional<double> T = std::nullopt,
                      double dt = PandaExample::dt);

}  // namespace se3pf
