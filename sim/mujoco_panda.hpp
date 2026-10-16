#pragma once
// The Panda in MuJoCo, at the torque level: the last check before hardware (docs/MUJOCO.md).
//
// Not part of the se3pf library, which never includes MuJoCo: this is the separate library se3pf_mujoco, built only
// with -DSE3PF_WITH_MUJOCO=ON. It has two parts:
//   1. MujocoPanda, a small RAII wrapper of one MuJoCo Menagerie Panda model (panda_nohand.xml or panda.xml) with its
//      seven arm actuators turned into torque motors: state, MuJoCo's own kinematics (for a cross-check of
//      PandaArm::output), the arm's dynamics M, C qdot + G and the passive forces, and two ways to step.
//   2. simulateMujoco, the closed loop of the paper's Panda example with the plant in MuJoCo, used by the app
//      panda_mujoco and by tests_mujoco.
#include "se3pf/types.hpp"
#include "se3pf/controller.hpp"
#include "se3pf/panda.hpp"
#include "se3pf/panda_example.hpp"

#include <memory>
#include <string>
#include <vector>

// MuJoCo's two main types, declared as mujoco.h declares them, so that this header does not need mujoco.h.
typedef struct mjModel_ mjModel;
typedef struct mjData_ mjData;

namespace se3pf {

    // ---------------------------------------------------------------------------------------------------------
    // The plant

    // One MuJoCo Menagerie Panda: an mjModel and its mjData, owned (RAII) and not copyable.
    //
    // The Menagerie actuators are position servos (force = kp (ctrl - q) - kd qdot). The constructor turns the seven
    // arm actuators into torque motors, force = ctrl, clamped to the Franka torque limits (+-87 / +-12 N m), so
    // ctrl = tau. The finger actuator of panda.xml is kept, and the constructor opens the hand (keyframe "home":
    // fingers at 0.04 m, ctrl = 255), so a plant with the hand carries it as a payload.
    //
    // Joint i is qpos[i], dof i and ctrl[i] (checked). The flange (Franka frame F) is body link7 moved 0.107 m along
    // its own z-axis, with link7's attitude (the attachment_site of panda_nohand.xml is turned 135 deg about z, so
    // its attitude is not used). MuJoCo's link frames are the modified-DH frames of the design doc.
    //
    // Physics that matter for the controller: joint armature 0.1 kg m^2 (part of M), joint damping 1 N m s/rad
    // (qfrc_passive = -D qdot; switched off unless Options::damping), no friction loss, no contacts along the lift.
    // Both integrators treat the damping explicitly, as a force in qacc = M^-1 (tau - bias + passive): the constructor
    // sets mjDSBL_EULERDAMP, because MuJoCo's Euler would otherwise integrate the damping implicitly, and the qacc it
    // reports would not be the acceleration it integrates. MuJoCo's automatic reset on a bad state is switched off,
    // so a blow-up shows as a warning and NaN, not as a silent restart.
    class MujocoPanda {
        public:
            static constexpr int kArm = 7;  // arm joints; panda.xml adds two finger joints after them

            enum class Integrator {
                Euler,  // MuJoCo's semi-implicit Euler: qdot += h qdd, then q += h qdot (with the new qdot); the
                        // damping is explicit in qdd (mjDSBL_EULERDAMP)
                RK4     // classic Runge-Kutta, with ctrl held over the step
            };

            struct Options {
                double timestep = 1e-3;                  // h, s (the Menagerie default is 2 ms)
                Integrator integrator = Integrator::Euler;
                bool damping = false;                    // joint damping 1 N m s/rad on or off (mjDSBL_DAMPER)
            };

            // Loads the model (MJCF). Throws std::invalid_argument for a timestep <= 0, and std::runtime_error if the
            // file cannot be loaded or is not a Panda whose actuator i drives joint i with unit gear (so that
            // ctrl[i] = tau_i).
            MujocoPanda(const std::string& xmlPath, const Options& options);
            // The default options. (Not a default argument "= Options()": a nested struct with default member
            // initializers cannot be used in its enclosing class before the class is complete; GCC says "default
            // member initializer required before the end of its enclosing class".)
            explicit MujocoPanda(const std::string& xmlPath) : MujocoPanda(xmlPath, Options()) {}
            // Owns its mjModel and mjData: a copy would free them twice. (The unique_ptr members also delete the
            // copy; saying so here is for the reader.)
            MujocoPanda(const MujocoPanda&) = delete;
            MujocoPanda& operator=(const MujocoPanda&) = delete;

            const Options& options() const { return options_; }
            bool hasHand() const;                    // panda.xml: nq = nv = 9 (two finger joints), nu = 8
            Vec7 torqueLimit() const;                // ctrlrange upper bounds: 87, 87, 87, 87, 12, 12, 12 N m

            // State of the arm joints. setState leaves the fingers (if any) as they are and evaluates the model at
            // the new state (mj_forward), so every quantity below belongs to it.
            void setState(const Vec7& q, const Vec7& qd);
            Vec7 q() const;
            Vec7 qd() const;
            Vec7 qacc() const;                       // acceleration of the last evaluation or step
            double time() const;

            // The arm's dynamics at the current state (one mj_forward): M (7 x 7 block, joint armature included),
            // bias = C(q, qdot) qdot + G(q) and passive = -D qdot. With the hand, M and bias include the hand's mass.
            void dynamics(Mat<7, 7>& M, Vec7& bias, Vec7& passive);

            // G(q) at any q with qdot = 0, in a scratch mjData (the state is not touched).
            Vec7 gravity(const Vec7& q);

            // MuJoCo's flange pose and body Jacobian (rows (omega, v) in the flange frame, J_b = [R^T jacr; R^T jacp]
            // of mj_jac at the flange point on link7), in the scratch mjData: to cross-check PandaArm::output.
            Pose flange(const Vec7& q, Mat67* Jb = nullptr);

            // Exact model (eq:fb): mj_step1 evaluates M, C qdot + G and the passive forces at the current state,
            // tau = M qdd + bias - passive goes to ctrl, and mj_step2 integrates one step of the Euler integrator:
            // qacc = qdd and qdot <- qdot + h qdd (to round-off), damping on or off, when no constraint is active.
            // Returns tau. Exact for panda_nohand.xml only: in panda.xml the finger dofs couple to the arm through M,
            // and the product M qdd gives them zero acceleration (simulateMujoco refuses an exact run with the hand).
            Vec7 exactStep(const Vec7& qdd);

            // Zero-order hold: ctrl = tau, held while MuJoCo takes `substeps` steps of the timestep (mj_step).
            void holdStep(const Vec7& tau, int substeps = 1);

            // An external joint torque on the arm (qfrc_applied), kept until changed: a push for the safety tests of
            // the hardware experiments (sim/mujoco_experiments.hpp). Zero by default.
            void setAppliedTorque(const Vec7& tau);

            // Constraints of the last evaluation (mj_forward, mj_step1 or the last stage of mj_step): contacts and
            // rows of the constraint problem (contacts, active joint limits, equalities such as the finger coupling
            // of panda.xml).
            int contacts() const;
            int constraintRows() const;
            // Number of MuJoCo warnings so far (bad qacc, singular inertia, ...): 0 in a healthy run.
            int warnings() const;
            // True if every |tau_i| is within the actuator's range, so MuJoCo applies tau without clamping.
            bool withinTorqueLimits(const Vec7& tau) const;

            // The raw MuJoCo objects, for anything this class does not wrap (read-only).
            const mjModel* model() const { return m_.get(); }
            const mjData* data() const { return d_.get(); }

        private:
            // RAII: mj_deleteModel / mj_deleteData run when the object goes away, also if the constructor throws
            // half-way. The deleters are defined in the .cpp, where mujoco.h is included.
            struct ModelDeleter { void operator()(mjModel* m) const; };
            struct DataDeleter { void operator()(mjData* d) const; };

            Options options_;
            std::unique_ptr<mjModel, ModelDeleter> m_;
            std::unique_ptr<mjData, DataDeleter> d_;        // the plant
            std::unique_ptr<mjData, DataDeleter> scratch_;  // for gravity and flange: does not disturb the plant
            // Buffers of the model's size (nv = 7 or 9), allocated once in the constructor
            std::vector<double> fullM_;                     // nv x nv, for mj_fullM
            std::vector<double> jacp_, jacr_;               // 3 x nv, for mj_jac
            std::vector<double> vecIn_, vecOut_;            // nv, for mj_mulM
            int link7_ = -1;                                // body id of link7
    };

    // The Menagerie Panda directory this build was configured with (third_party/mujoco_menagerie/franka_emika_panda).
    std::string defaultPandaModelDir();

    // ---------------------------------------------------------------------------------------------------------
    // The hardware margins

    // The margins of the hardware section (design doc, Hardware; their bounds are HardwareBounds in se3pf/panda.hpp) at
    // one control instant. Every joint is measured against its own limit (PandaArm::dqMax, ddqMax, tauMax: the wrist
    // joints 5-7 have 12 N m against 87 N m, and their own velocity and acceleration limits), so the per-joint
    // fractions are kept: the largest is the margin, and the joint that sets it can change along a path.
    struct HardwareMargins {
        double jointMargin = 0.0;            // min_i distance of q_i to its limits, rad (PandaArm::jointMargin)
        int marginJoint = -1;                // the joint that sets it
        Vec7 dqFrac = Vec7::Zero();          // |qdot_i| / qdot_max,i                       (bound 50 %)
        Vec7 ddqFrac = Vec7::Zero();         // |qdd*_i| / qdd_max,i, commanded             (bound 30 %)
        Vec7 tauFrac = Vec7::Zero();         // |tau_i| / tau_max,i, delivered (G included) (bound 60 %)
        double tauRateCmd = 0.0;             // max_i |d(tau_i - G_i)/dt|, N m/s: the torque sent to libfranka
                                             // (bound 60 % of PandaArm::tauRateMax)
        double tauRate = 0.0;                // max_i |d tau_i/dt|, N m/s: the full torque (for comparison)
        double flangeSpeed = 0.0;            // |v| of the flange (J_b qdot, translational part), m/s (bound 0.5)
    };

    // The margins at (q, qdot) with the commanded acceleration qdd*, the delivered torque tau and its gravity part G,
    // and the torque rates as backward differences over the control period dt from the previous instant's tau and G.
    // (Before the first instant the arm rests under gravity compensation: tauPrev = GPrev = G(q0), so the first
    // commanded rate is |tau - G| / dt.) No allocation.
    HardwareMargins hardwareMargins(const Vec7& q, const Vec7& qd, const Vec7& qdd, const Vec7& tau, const Vec7& G,
                                    const Vec7& tauPrev, const Vec7& GPrev, double dt);

    // ---------------------------------------------------------------------------------------------------------
    // The closed loop

    // One run of the Panda example with the plant in MuJoCo. Every control period dt: Controller<PandaArm>::step at
    // the plant's (q, qdot) gives qdd*, a torque goes to the plant, and the plant moves on.
    struct MujocoRunOptions {
        enum class Start {
            Paper,       // at rest at Psi(0.1, xi0) with beta = 0 (paperStart), for `duration`
            SelfMotion   // at gamma(0.1) with beta = 0.2 and A(q0) qdot0 = (nu_d, 0, ..., 0), for one loop L / nu_d
        };
        enum class Torque {
            Exact,  // qdd* held for the period, and MujocoPanda::exactStep at every plant step: tau = M qdd* +
                    // C qdot + G - passive at the plant's state, so qacc = qdd* exactly (Euler integrator;
                    // panda_nohand.xml only, see exactStep)
            Hold    // tau = M qdd* + C qdot + G from the state at the start of the period, held for dt (what the
                    // robot does). The passive forces (damping) are not compensated: they stand for the robot's
                    // friction, which the controller does not model either
        };

        std::string plant;              // model file of the plant
        std::string model;              // model file for M and C qdot + G of the torque; empty: the plant itself
        bool plantGravity = false;      // with a model file: G from the plant instead (the robot adds G of its own
                                        // model, which includes the end effector configured in Desk)
        Start start = Start::Paper;
        Torque torque = Torque::Exact;
        double duration = PandaExample::T;  // s (Start::Paper only; SelfMotion runs one loop)
        double dt = PandaExample::dt;   // control period, s
        int substeps = 1;               // plant steps per period; the plant's timestep is h = dt / substeps
        MujocoPanda::Integrator integrator = MujocoPanda::Integrator::RK4;  // plant integrator for Hold
        bool damping = false;           // joint damping 1 N m s/rad in the plant (and in a separate model)
        double rampTime = 0.0;          // mu ramp T (Controller::setRampTime); 0 = off
        bool baseline = false;          // also run the same controller at the joint-acceleration level
        bool crossCheck = false;        // compare PandaArm::output with MuJoCo's kinematics at every step
        bool keepLog = true;            // keep one MujocoLogRow per period
    };

    // One row of the log, at a control instant t_k = k dt, before the plant moves on.
    struct MujocoLogRow {
        double t = 0.0;
        double eta = 0.0, etaDot = 0.0;
        Vec5 xi = Vec5::Zero();
        double xiNorm = 0.0, xiDesign = 0.0;     // |xi| and |xi0| (1 + 4t) e^{-4t}
        double beta = 0.0, betaDesign = 0.0;     // beta and beta0 (1 + 4t) e^{-4t}
        Vec7 q = Vec7::Zero(), qd = Vec7::Zero();
        Vec7 tau = Vec7::Zero();                 // ctrl over [t, t + dt): the full joint torque, gravity included
        Vec7 G = Vec7::Zero();                   // the gravity term inside tau (tau - G is what libfranka is sent)
        Vec7 qdd = Vec7::Zero();                 // qdd* of the controller
        double sigmaMin = 0.0;                   // sigma_min(A)
        double jointMargin = 0.0;                // PandaArm::jointMargin(q)
        double flangeSpeed = 0.0;                // |v| of the flange, m/s
        double tauRateCmd = 0.0;                 // max_i |d(tau_i - G_i)/dt|: the commanded torque on the robot
        double tauRate = 0.0;                    // max_i |d tau_i / dt|: the full torque
        int contacts = 0, constraintRows = 0;
        double ramp = 1.0;                       // r(t / T)
        double xiNormKin = 0.0, betaKin = 0.0;   // the joint-acceleration-level baseline (0 without it)
        double qErrKin = 0.0;                    // max_i |q_i - q_kin,i|
    };

    // The log and the summary numbers of one run. Peaks are over every control instant of the run.
    struct MujocoRun {
        MujocoRunOptions options;
        std::vector<MujocoLogRow> log;
        long steps = 0;                          // control periods taken
        bool ok = false;                         // false if a step failed (left the tube, singular A, NaN) or
                                                 // MuJoCo warned
        double failureTime = -1.0;
        double wallTime = 0.0;                   // s

        double xi0 = 0.0, beta0 = 0.0;           // start values of |xi| and beta
        double maxDeviation = 0.0, maxDeviationTime = 0.0;   // max ||xi| - design|
        double maxXi = 0.0;                      // max |xi|
        double maxBeta = 0.0;                    // max |beta|
        double maxBetaDeviation = 0.0;           // max |beta - beta0 (1 + 4t) e^{-4t}|
        double finalXi = 0.0, finalBeta = 0.0, finalEtaDot = 0.0;
        double steadyXi = 0.0, steadyBeta = 0.0; // max |xi| and |beta| over the last 5 s (all of a shorter run)
        double finalLiftError = 0.0;             // max_i |q_i - q*_i(eta)| after the last period
        double finalLiftError2 = 0.0;            // its 2-norm
        int maxIterations = 0;                   // Newton updates per projection

        // Margins of the hardware section (design doc), as peaks of HardwareMargins over the run
        double minJointMargin = 1e9;
        int minMarginJoint = -1;
        double maxDqFrac = 0.0;                  // max_i |qdot_i| / qdot_max,i
        double maxDdqFrac = 0.0;                 // max_i |qdd*_i| / qdd_max,i (commanded)
        double maxDdqFracTime = 0.0;
        double maxTauFrac = 0.0;                 // max_i |tau_i| / tau_max,i (delivered, gravity included)
        Vec7 peakDqFrac = Vec7::Zero();          // the same per joint: entry i is joint i's peak over the run
        Vec7 peakDdqFrac = Vec7::Zero();         // (so maxDqFrac = peakDqFrac.maxCoeff(), and so on)
        Vec7 peakTauFrac = Vec7::Zero();
        double maxTauRateCmd = 0.0;              // max |d(tau - G)/dt|, N m/s
        double maxTauRate = 0.0;                 // max |d tau/dt|, N m/s
        double maxFlangeSpeed = 0.0;             // m/s
        double minSigma = 1e9;                   // min sigma_min(A)
        int maxContacts = 0, maxConstraintRows = 0;
        long clampedSteps = 0;                   // periods with a torque outside the actuator range
        int warnings = 0;                        // MuJoCo warnings

        double maxQaccError = 0.0;               // Torque::Exact: max |qacc - qdd*| (round-off)

        // The baseline: the same controller with q += dt qdot + dt^2/2 qdd*, qdot += dt qdd* (design doc, Panda
        // simulation), from the same start.
        double maxBaselineQError = 0.0;          // max_i |q_i - q_kin,i| over the run
        double baselineMaxDeviation = 0.0, baselineMaxBeta = 0.0, baselineMaxBetaDeviation = 0.0;
        double baselineMaxXi = 0.0, baselineFinalXi = 0.0, baselineFinalEtaDot = 0.0;
        double baselineFinalLiftError = 0.0, baselineFinalLiftError2 = 0.0;
        double baselineMaxDqFrac = 0.0, baselineMaxDdqFrac = 0.0;

        // The kinematics cross-check (crossCheck): max differences of PandaArm::output and MuJoCo
        double maxPoseError = 0.0;               // entries of the 4 x 4 pose
        double maxJacobianError = 0.0;           // entries of J_b
    };

    // The controller of the paper's Panda example: PandaExample::gains() (poles at -4, k_t = 5, K_p,beta = 16,
    // K_d,beta = 8), PandaExample::reference() (nu_d = 0.08) and the lift's beta (Lemma 3); reset to eta-hat = s0 =
    // 0.1. The example must outlive it.
    Controller<PandaArm> exampleController(const PandaExample& example);

    // Runs the closed loop (offline: the log is reserved before the loop, and nothing in the loop allocates or
    // prints). The run has stepCount(T, dt) + 1 control periods, T = duration (Paper) or one loop L / nu_d
    // (SelfMotion). Throws std::invalid_argument for bad options (also a T or dt that is not finite, or T / dt above
    // PandaExample::maxSteps) and std::runtime_error if a model cannot be loaded.
    MujocoRun simulateMujoco(const PandaExample& example, const MujocoRunOptions& options);

}  // namespace se3pf
