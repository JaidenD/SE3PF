// Tests of the MuJoCo layer (sim/mujoco_panda.hpp, docs/MUJOCO.md): the Menagerie Panda in torque mode, and short
// closed loops of the Panda example at the torque level. A separate program, tests_mujoco, built only with
// -DSE3PF_WITH_MUJOCO=ON, so the default tests never need MuJoCo. About 2 s in Release. The closed loops run the
// controller (1-3 ms per step in Debug), so the default Debug run (about 23 s) hides the two slowest ones ([.slow],
// see SLOW_IN_DEBUG; about 45 s more in Debug).
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mujoco_panda.hpp"
#include "se3pf/panda.hpp"
#include "panda_fixture.hpp"
#include "test_helpers.hpp"

#include <mujoco/mujoco.h>
#include <unistd.h>  // getpid

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

using namespace se3pf;
using Catch::Matchers::WithinAbs;
using testing::maxDiff;
using testing::randn;
using testing::sci;
using Start = MujocoRunOptions::Start;
using Torque = MujocoRunOptions::Torque;
using B = HardwareBounds;

namespace {
    std::string modelFile(const char* name) { return defaultPandaModelDir() + "/" + name; }
    std::string noHand() { return modelFile("panda_nohand.xml"); }
    std::string withHand() { return modelFile("panda.xml"); }

    // Built once (panda_fixture.hpp; about 4 s in Debug, 0.03 s in Release)
    const PandaExample& example() { return testing::pandaExample(); }

    // Random configuration within the joint limits
    Vec7 randomConfiguration() { return testing::randomConfig(); }

    // Options of the two closed loops of the app: (A) exact model and (C) hardware-like
    MujocoRunOptions exactRun(double T)
    {
        MujocoRunOptions o;
        o.plant = noHand();
        o.torque = Torque::Exact;
        o.duration = T;
        o.baseline = true;
        o.crossCheck = true;
        return o;
    }

    MujocoRunOptions hardwareRun(double T)
    {
        MujocoRunOptions o;
        o.plant = noHand();
        o.torque = Torque::Hold;
        o.integrator = MujocoPanda::Integrator::RK4;
        o.damping = true;
        o.rampTime = 1.0;
        o.duration = T;
        return o;
    }
}  // namespace

// ---------------------------------------------------------------------------------------------------------------
// The plant

TEST_CASE("MujocoPanda: loads the arm with torque motors; options; the hand model; refuses what is not a Panda")
{
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<MujocoPanda>);
    STATIC_REQUIRE_FALSE(std::is_copy_assignable_v<MujocoPanda>);

    MujocoPanda arm(noHand());
    const mjModel* m = arm.model();
    CHECK_FALSE(arm.hasHand());
    CHECK(m->nq == 7);
    CHECK(m->nv == 7);
    CHECK(m->nu == 7);
    // ctrl = tau: unit gain, no bias, clamped to the Franka torque limits
    CHECK(arm.torqueLimit() == PandaArm::tauMax());
    for (int i = 0; i < 7; ++i) {
        CHECK(m->actuator_gaintype[i] == mjGAIN_FIXED);
        CHECK(m->actuator_gainprm[i * mjNGAIN] == 1.0);
        CHECK(m->actuator_biastype[i] == mjBIAS_NONE);
        CHECK(m->actuator_ctrlrange[2 * i] == -PandaArm::tauMax()[i]);
    }
    // Default options: 1 ms, Euler, damping off; Euler's damping explicit (if on); MuJoCo's automatic reset off
    CHECK(m->opt.timestep == 1e-3);
    CHECK(m->opt.integrator == mjINT_EULER);
    CHECK((m->opt.disableflags & mjDSBL_DAMPER) != 0);
    CHECK((m->opt.disableflags & mjDSBL_EULERDAMP) != 0);
    CHECK((m->opt.disableflags & mjDSBL_AUTORESET) != 0);
    CHECK(m->dof_damping[0] == 1.0);     // the model's damping, switched off by the flag
    CHECK(m->dof_armature[0] == 0.1);

    MujocoPanda::Options o;
    o.timestep = 0.5e-3;
    o.integrator = MujocoPanda::Integrator::RK4;
    o.damping = true;
    MujocoPanda rk(noHand(), o);
    CHECK(rk.model()->opt.timestep == 0.5e-3);
    CHECK(rk.model()->opt.integrator == mjINT_RK4);
    CHECK((rk.model()->opt.disableflags & mjDSBL_DAMPER) == 0);

    // panda.xml: two finger joints and the finger actuator, kept as a position servo; the hand starts open
    MujocoPanda hand(withHand());
    CHECK(hand.hasHand());
    CHECK(hand.model()->nq == 9);
    CHECK(hand.model()->nu == 8);
    CHECK(hand.model()->actuator_biastype[7] == mjBIAS_AFFINE);
    CHECK(hand.data()->qpos[7] == 0.04);
    CHECK(hand.data()->qpos[8] == 0.04);
    CHECK(hand.data()->ctrl[7] == 255.0);
    for (int i = 0; i < 7; ++i)
        CHECK(hand.data()->ctrl[i] == 0.0);  // the keyframe's servo targets would be torques now

    CHECK_THROWS_AS(MujocoPanda(noHand(), MujocoPanda::Options{0.0, MujocoPanda::Integrator::Euler, false}),
                    std::invalid_argument);
    CHECK_THROWS_AS(MujocoPanda(modelFile("no_such_file.xml")), std::runtime_error);

    // A valid MJCF model that is not a Panda: one pendulum with a motor (a file name of this process, so that test
    // programs run in parallel do not share it)
    const std::filesystem::path file =
        std::filesystem::temp_directory_path() / ("se3pf_not_a_panda_" + std::to_string(::getpid()) + ".xml");
    {
        std::ofstream f(file);
        f << "<mujoco><worldbody><body name=\"link1\"><joint name=\"j\"/><geom size=\"0.1\" mass=\"1\"/></body>"
             "</worldbody><actuator><motor joint=\"j\"/></actuator></mujoco>";
    }
    CHECK_THROWS_AS(MujocoPanda(file.string()), std::runtime_error);
    std::filesystem::remove(file);
}

TEST_CASE("MujocoPanda: flange pose and body Jacobian agree with PandaArm::output (FK 1e-12, J_b 1e-12)")
{
    MujocoPanda arm(noHand());
    MujocoPanda hand(withHand());  // same link frames; the hand body hangs below the flange

    // Anchors of the design doc: FK(0) = (0.088, 0, 0.926) with R = diag(1, -1, -1)
    Pose g0 = Pose::Identity();
    g0.block<3,3>(0,0) = Vector3d(1.0, -1.0, -1.0).asDiagonal();
    g0.block<3,1>(0,3) = Vector3d(0.088, 0.0, 0.926);
    CHECK(maxDiff(arm.flange(Vec7::Zero()), g0) < 1e-12);

    double poseErr = 0.0, jacErr = 0.0, handErr = 0.0;
    for (int k = 0; k < 200; ++k) {
        const Vec7 q = randomConfiguration();
        Mat67 Jb = Mat67::Zero(), JbMj = Mat67::Zero(), JbHand = Mat67::Zero();
        const Pose g = PandaArm::output(q, &Jb);
        // flange first, then compare its Jacobian: the order in which function arguments are evaluated is
        // unspecified, so std::max(maxDiff(flange(q, &J), g), maxDiff(J, Jb)) may read J before flange fills it.
        const Pose gMj = arm.flange(q, &JbMj);
        const Pose gHand = hand.flange(q, &JbHand);
        poseErr = std::max(poseErr, maxDiff(gMj, g));
        jacErr = std::max(jacErr, maxDiff(JbMj, Jb));
        handErr = std::max({handErr, maxDiff(gHand, g), maxDiff(JbHand, Jb)});
    }
    CAPTURE(sci(poseErr), sci(jacErr), sci(handErr));
    CHECK(poseErr < 1e-12);
    CHECK(jacErr < 1e-12);
    CHECK(handErr < 1e-12);

    // flange works in scratch data: the plant's state is untouched
    const Vec7 q = PandaArm::restPosture();
    arm.setState(q, Vec7::Zero());
    arm.flange(randomConfiguration());
    CHECK(arm.q() == q);
}

TEST_CASE("MujocoPanda: M is symmetric positive definite with the armature; bias at rest is G; passive = -D qdot")
{
    MujocoPanda::Options o;
    o.damping = true;
    MujocoPanda arm(noHand(), o);
    const Vec7 q = PandaArm::restPosture() + randn<Vec7>(0.2);
    const Vec7 qd = randn<Vec7>(0.5);

    arm.setState(q, Vec7::Zero());
    Mat<7, 7> M = Mat<7, 7>::Zero();
    Vec7 bias = Vec7::Zero(), passive = Vec7::Zero();
    arm.dynamics(M, bias, passive);
    CHECK(maxDiff(M, M.transpose()) < 1e-15);
    CHECK(M.diagonal().minCoeff() > 0.1);  // armature 0.1 kg m^2 on every joint
    CHECK(Eigen::SelfAdjointEigenSolver<Mat<7, 7>>(M).eigenvalues().minCoeff() > 0.0);
    CHECK(maxDiff(bias, arm.gravity(q)) < 1e-12);  // C qdot = 0 at rest
    CHECK(passive.isZero(0.0));

    // Moving: passive = -D qdot with D = 1 N m s/rad; the bias now has C qdot, and gravity(q) still ignores qdot
    arm.setState(q, qd);
    arm.dynamics(M, bias, passive);
    CHECK(maxDiff(passive, -qd) < 1e-12);
    CHECK(maxDiff(bias, arm.gravity(q)) > 1e-3);

    // With the damper off, no passive force
    MujocoPanda off(noHand());
    off.setState(q, qd);
    off.dynamics(M, bias, passive);
    CHECK(passive.isZero(0.0));

    // The dense M of dynamics (mj_fullM) against MuJoCo's own product M qdd (mj_mulM, used by exactStep): the torque
    // M qdd + bias - passive must be the one exactStep applies, off-diagonal coupling included.
    arm.setState(q, qd);
    arm.dynamics(M, bias, passive);
    const Vec7 qdd = randn<Vec7>(1.0);
    const Vec7 tau = M * qdd + bias - passive;
    CHECK(maxDiff(arm.exactStep(qdd), tau) < 1e-12);
}

TEST_CASE("MujocoPanda: exactStep gives qacc = qdd* and qdot += h qdd* to 1e-12, damping on or off (eq:fb)")
{
    for (const bool damping : {false, true}) {
        CAPTURE(damping);
        MujocoPanda::Options o;
        o.damping = damping;
        MujocoPanda arm(noHand(), o);
        arm.setState(PandaArm::restPosture(), randn<Vec7>(0.3));
        double accErr = 0.0, velErr = 0.0;
        for (int k = 0; k < 200; ++k) {
            const Vec7 qdd = randn<Vec7>(1.0);
            const Vec7 qd = arm.qd();
            const Vec7 tau = arm.exactStep(qdd);
            CHECK(arm.withinTorqueLimits(tau));
            accErr = std::max(accErr, maxDiff(arm.qacc(), qdd));
            // Euler: qdot_{k+1} = qdot_k + h qdd. The reported qacc alone would not show an implicit damping step
            // (MuJoCo's default for Euler, switched off with mjDSBL_EULERDAMP): it integrates (M + h D)^-1 (...),
            // which misses this update by about h * 0.03 rad/s^2 at these speeds.
            velErr = std::max(velErr, maxDiff(arm.qd(), qd + 1e-3 * qdd));
        }
        CAPTURE(sci(accErr), sci(velErr));
        CHECK(accErr < 1e-12);
        CHECK(velErr < 1e-12);
        CHECK(arm.contacts() == 0);
        CHECK(arm.constraintRows() == 0);
        CHECK(arm.warnings() == 0);
        CHECK_THAT(arm.time(), WithinAbs(0.2, 1e-12));
    }
}

TEST_CASE("MujocoPanda: tau = qfrc_bias holds the rest posture for 2 s; tau = 0 lets it fall")
{
    for (const auto integrator : {MujocoPanda::Integrator::Euler, MujocoPanda::Integrator::RK4}) {
        MujocoPanda::Options o;
        o.integrator = integrator;
        MujocoPanda arm(noHand(), o);
        const Vec7 qr = PandaArm::restPosture();
        arm.setState(qr, Vec7::Zero());
        Mat<7, 7> M = Mat<7, 7>::Zero();
        Vec7 bias = Vec7::Zero(), passive = Vec7::Zero();
        double drift = 0.0;
        for (int k = 0; k < 2000; ++k) {
            arm.dynamics(M, bias, passive);
            arm.holdStep(bias);
            drift = std::max(drift, (arm.q() - qr).cwiseAbs().maxCoeff());
        }
        CAPTURE(sci(drift));
        CHECK(drift < 1e-10);

        arm.setState(qr, Vec7::Zero());
        arm.holdStep(Vec7::Zero(), 100);  // 0.1 s without torque
        CHECK((arm.q() - qr).cwiseAbs().maxCoeff() > 1e-2);
    }

    // The torque check: MuJoCo would clamp a torque outside the range
    MujocoPanda arm(noHand());
    Vec7 tau = 0.99 * PandaArm::tauMax();
    CHECK(arm.withinTorqueLimits(tau));
    CHECK(arm.withinTorqueLimits(-tau));
    tau[4] = 12.5;
    CHECK_FALSE(arm.withinTorqueLimits(tau));
    tau[4] = std::nan("");
    CHECK_FALSE(arm.withinTorqueLimits(tau));
}

TEST_CASE("MujocoPanda: counts contacts (a folded arm) and constraint rows (a joint past its limit)")
{
    MujocoPanda arm(noHand());
    arm.setState(PandaArm::restPosture(), Vec7::Zero());
    CHECK(arm.contacts() == 0);
    CHECK(arm.constraintRows() == 0);

    // Joint 4 folded to -3.0 rad with joint 2 at 0.5 rad: link 7 hits the lower links
    Vec7 folded = Vec7::Zero();
    folded[1] = 0.5;
    folded[3] = -3.0;
    arm.setState(folded, Vec7::Zero());
    CHECK(arm.contacts() > 0);
    CHECK(arm.constraintRows() > 0);

    // Joint 4 past its upper limit -0.0698 rad: one limit row, no contact
    Vec7 q = PandaArm::restPosture();
    q[3] = -0.05;
    arm.setState(q, Vec7::Zero());
    CHECK(arm.contacts() == 0);
    CHECK(arm.constraintRows() == 1);
}

// ---------------------------------------------------------------------------------------------------------------
// The hardware margins

TEST_CASE("hardwareMargins: each joint against its own limit (each joint sets each margin once)")
{
    // Along the paper's path joint 2 sets the torque peak, so a wrong wrist limit (12 N m) would never show in a
    // closed loop. Here every joint in turn sits at 50 % of its own limit and the others at 20 % of theirs: one
    // limit for all joints (87 N m, say) would put a wrist joint at 7 % and pick the wrong joint.
    const Vec7 lo = PandaArm::qMin(), hi = PandaArm::qMax();
    const Vec7 dqMax = PandaArm::dqMax(), ddqMax = PandaArm::ddqMax(), tauMax = PandaArm::tauMax();
    const Vec7 G = randn<Vec7>(10.0);
    for (int j = 0; j < 7; ++j) {
        CAPTURE(j);
        Vec7 q = 0.5 * (lo + hi);   // every joint at least 1.5 rad from its limits ...
        q[j] = hi[j] - 0.05;        // ... except joint j
        Vec7 qd = 0.2 * dqMax, qdd = -0.2 * ddqMax, tau = 0.2 * tauMax;
        qd[j] = -0.5 * dqMax[j];
        qdd[j] = 0.5 * ddqMax[j];
        tau[j] = -0.5 * tauMax[j];
        Vec7 expected = Vec7::Constant(0.2);
        expected[j] = 0.5;

        const HardwareMargins h = hardwareMargins(q, qd, qdd, tau, G, tau, G, 1e-3);
        CHECK_THAT(h.jointMargin, WithinAbs(0.05, 1e-12));
        CHECK(h.marginJoint == j);
        CHECK(maxDiff(h.dqFrac, expected) < 1e-15);
        CHECK(maxDiff(h.ddqFrac, expected) < 1e-15);
        CHECK(maxDiff(h.tauFrac, expected) < 1e-15);
        int i = -1;
        h.tauFrac.maxCoeff(&i);
        CHECK(i == j);
        CHECK(h.tauRateCmd == 0.0);  // tau and G as at the previous instant
        CHECK(h.tauRate == 0.0);
    }
}

TEST_CASE("hardwareMargins: torque rates of tau - G (what libfranka is sent) and of tau; flange speed")
{
    const double dt = 1e-3;
    const Vec7 q = PandaArm::restPosture();
    const Vec7 GPrev = randn<Vec7>(10.0), tauPrev = GPrev + randn<Vec7>(1.0);

    // In one period the gravity term of joint 2 grows by 10 N m/s while the commanded part of joint 6 grows by
    // 3 N m/s: the full torque changes at 10 N m/s, the commanded torque at 3.
    Vec7 dG = Vec7::Zero(), dCmd = Vec7::Zero();
    dG[1] = 10.0 * dt;
    dCmd[5] = -3.0 * dt;
    const Vec7 G = GPrev + dG, tau = tauPrev + dG + dCmd;
    HardwareMargins h = hardwareMargins(q, Vec7::Zero(), Vec7::Zero(), tau, G, tauPrev, GPrev, dt);
    CAPTURE(sci(h.tauRateCmd), sci(h.tauRate));
    CHECK_THAT(h.tauRateCmd, WithinAbs(3.0, 1e-9));
    CHECK_THAT(h.tauRate, WithinAbs(10.0, 1e-9));

    // The first instant, after rest under gravity compensation (tauPrev = GPrev = G(q0)): the commanded torque
    // jumps from 0 to tau - G
    h = hardwareMargins(q, Vec7::Zero(), Vec7::Zero(), tau, G, G, G, dt);
    CHECK_THAT(h.tauRateCmd, WithinAbs((tau - G).cwiseAbs().maxCoeff() / dt, 1e-9));

    // Flange speed: |v| of the body twist J_b qdot. qdot = J_b^+ (omega, v) with |v| = 0.5 m/s and a rotation
    // that must not count.
    Mat67 Jb = Mat67::Zero();
    PandaArm::output(q, &Jb);
    Vec6 zeta = Vec6::Zero();
    zeta << 0.4, -0.2, 0.3, 0.3, 0.0, -0.4;
    const Vec7 qd = Jb.transpose() * (Jb * Jb.transpose()).ldlt().solve(zeta);
    h = hardwareMargins(q, qd, Vec7::Zero(), tau, G, tau, G, dt);
    CHECK_THAT(h.flangeSpeed, WithinAbs(0.5, 1e-12));
}

// ---------------------------------------------------------------------------------------------------------------
// The closed loop

TEST_CASE("simulateMujoco: argument checks")
{
    const PandaExample& ex = example();
    MujocoRunOptions o = exactRun(0.01);
    o.dt = 0.0;
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::invalid_argument);
    o = exactRun(0.01);
    o.substeps = 0;
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::invalid_argument);
    o = exactRun(0.01);
    o.model = noHand();  // the exact model is the plant's own
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::invalid_argument);
    o = exactRun(0.01);
    o.plant = withHand();  // the exact step is exact for the arm alone
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::invalid_argument);
    o = hardwareRun(0.01);
    o.plantGravity = true;  // needs a separate model
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::invalid_argument);
    o = hardwareRun(0.01);
    o.plant.clear();
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::invalid_argument);
    o = hardwareRun(0.01);
    o.plant = modelFile("no_such_file.xml");
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::runtime_error);
    // A run length that is not finite or has more than PandaExample::maxSteps periods (stepCount)
    o = hardwareRun(std::numeric_limits<double>::infinity());
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::invalid_argument);
    o = hardwareRun(1e15);
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::invalid_argument);
    o = exactRun(0.0);
    o.start = Start::SelfMotion;
    o.dt = 1e-320;  // one loop L / nu_d would be 1e321 periods
    CHECK_THROWS_AS(simulateMujoco(ex, o), std::invalid_argument);
}

// The design doc's Panda run (Closed loops: max deviation 2.43e-5, max |beta| 1.5e-5) at the torque level, first
// 2 s, which hold both peaks (at t = 0.29 s and 0.95 s in MuJoCo, 0.41 s at the joint-acceleration level).
// Measured: max ||xi| - design| 3.89e-5, max |beta| 3.58e-5. They exceed the joint-acceleration-level values by
// MuJoCo's semi-implicit Euler, which moves q by h^2 qdd* more per step than the exact update (an O(h) difference:
// halving h halves |q_mj - q_kin|, checked below). The bounds 5e-5 leave 25 % margin over the measured values; a
// torque error, such as a missing C qdot or a missing armature, gives xi errors far above them.
TEST_CASE("simulateMujoco: 2 s, exact model from the paper start: |xi| within 5e-5 of the design, |beta| < 5e-5")
{
    const PandaExample& ex = example();
    const MujocoRun run = simulateMujoco(ex, exactRun(2.0));
    REQUIRE(run.ok);
    REQUIRE(run.steps == 2001);
    REQUIRE(run.log.size() == 2001);

    CAPTURE(sci(run.maxDeviation), sci(run.maxBeta), sci(run.maxQaccError));
    CHECK(run.maxQaccError < 1e-12);
    CHECK(run.maxDeviation < 5e-5);
    CHECK(run.maxBeta < 5e-5);
    CHECK(run.maxPoseError < 1e-12);
    CHECK(run.maxJacobianError < 1e-12);
    CHECK(run.maxContacts == 0);
    CHECK(run.maxConstraintRows == 0);
    CHECK(run.clampedSteps == 0);
    CHECK(run.warnings == 0);
    CHECK(run.maxIterations == 1);

    // The baseline is the design doc's run: its peaks fall in the first 2 s
    CAPTURE(sci(run.baselineMaxDeviation), sci(run.baselineMaxBeta));
    CHECK_THAT(run.baselineMaxDeviation, WithinAbs(2.43e-5, 0.01e-5));
    CHECK_THAT(run.baselineMaxBeta, WithinAbs(1.52e-5, 0.01e-5));

    // The first command, at rest (no ramp): 49 % of a joint's acceleration limit (design doc, Hardware)
    CHECK_THAT(run.maxDdqFrac, WithinAbs(0.494, 0.002));
    CHECK(run.maxDdqFracTime == 0.0);

    // The log against the summary
    double maxDev = 0.0, maxBeta = 0.0;
    for (std::size_t k = 0; k < run.log.size(); ++k) {
        const MujocoLogRow& r = run.log[k];
        REQUIRE(r.t == k * 1e-3);
        REQUIRE(r.xiNorm == r.xi.norm());
        maxDev = std::max(maxDev, std::abs(r.xiNorm - r.xiDesign));
        maxBeta = std::max(maxBeta, std::abs(r.beta));
    }
    CHECK(run.maxDeviation == maxDev);
    CHECK(run.maxBeta == maxBeta);
    CHECK_THAT(run.log.front().xiNorm, WithinAbs(PandaExample::xi0().norm(), 1e-9));
    CHECK(run.log.front().q == ex.paperStart());
}

TEST_CASE("simulateMujoco: halving MuJoCo's step halves |q_mj - q_kin| (an O(h) integrator effect)", SLOW_IN_DEBUG)
{
    const PandaExample& ex = example();
    MujocoRunOptions o = exactRun(1.0);
    o.crossCheck = false;
    o.keepLog = false;
    const MujocoRun h1 = simulateMujoco(ex, o);
    o.substeps = 2;
    const MujocoRun h2 = simulateMujoco(ex, o);
    REQUIRE(h1.ok);
    REQUIRE(h2.ok);
    const double ratio = h1.maxBaselineQError / h2.maxBaselineQError;
    CAPTURE(sci(h1.maxBaselineQError), sci(h2.maxBaselineQError), ratio);
    CHECK(h1.maxBaselineQError > 1e-5);  // Euler's h^2 qdd* per step is visible ...
    CHECK(ratio > 1.8);                  // ... and first order in h
    CHECK(ratio < 2.2);
    CHECK(h2.maxQaccError < 1e-12);      // the exact torque at every substep
}

// Zero-order hold (run B of the app): tau = M qdd* + C qdot + G from the state at the start of each 1 ms period, held
// while RK4 integrates; damping off. The held torque is exact only at the start of the period, so the outputs
// leave their designed responses by an O(dt) amount: over 2 s, max ||xi| - design| 3.2e-5 and max |beta| 2.6e-4
// (the app shows both halve with dt). The bounds leave about 50 % margin; a torque without M, or with M's
// diagonal only, misses them by far.
TEST_CASE("simulateMujoco: 2 s, zero-order hold: |xi| within 5e-5 of the design, |beta| < 4e-4")
{
    const PandaExample& ex = example();
    MujocoRunOptions o = hardwareRun(2.0);
    o.damping = false;
    o.rampTime = 0.0;
    const MujocoRun run = simulateMujoco(ex, o);
    REQUIRE(run.ok);
    CAPTURE(sci(run.maxDeviation), sci(run.maxBeta));
    CHECK(run.maxDeviation < 5e-5);
    CHECK(run.maxBeta < 4e-4);
    CHECK(run.maxBeta > 1e-4);  // the hold is visible: beta is 15 times its value with the exact model
    CHECK(run.maxContacts == 0);
    CHECK(run.maxConstraintRows == 0);
}

namespace {
    // The peaks of the summary, recomputed from the log: every margin of the hardware section from q, qdot, qdd*,
    // tau and G alone, per joint and overall, without hardwareMargins. The torque rates are backward differences;
    // before the first period the arm rests under gravity compensation, tau = G = G(q0) (the first row's G).
    void checkPeaks(const MujocoRun& run, double dt)
    {
        REQUIRE(run.log.size() == static_cast<std::size_t>(run.steps));
        REQUIRE(!run.log.empty());
        double margin = 1e9, speed = 0.0, xi = 0.0, beta = 0.0, sigma = 1e9;
        double rateCmd = 0.0, rate = 0.0, rateCmdRowErr = 0.0, rateRowErr = 0.0;
        Vec7 dq = Vec7::Zero(), ddq = Vec7::Zero(), tau = Vec7::Zero();
        for (std::size_t k = 0; k < run.log.size(); ++k) {
            const MujocoLogRow& r = run.log[k];
            margin = std::min(margin, PandaArm::jointMargin(r.q));
            for (int i = 0; i < 7; ++i) {
                dq[i] = std::max(dq[i], std::abs(r.qd[i]) / PandaArm::dqMax()[i]);
                ddq[i] = std::max(ddq[i], std::abs(r.qdd[i]) / PandaArm::ddqMax()[i]);
                tau[i] = std::max(tau[i], std::abs(r.tau[i]) / PandaArm::tauMax()[i]);
            }
            Mat67 Jb = Mat67::Zero();
            PandaArm::output(r.q, &Jb);
            const Vec6 zeta = Jb * r.qd;
            speed = std::max(speed, zeta.tail<3>().norm());
            xi = std::max(xi, r.xiNorm);
            beta = std::max(beta, std::abs(r.beta));
            sigma = std::min(sigma, r.sigmaMin);

            const Vec7 tauPrev = k > 0 ? run.log[k - 1].tau : run.log[0].G;
            const Vec7 GPrev = k > 0 ? run.log[k - 1].G : run.log[0].G;
            const double rowRateCmd = ((r.tau - r.G) - (tauPrev - GPrev)).cwiseAbs().maxCoeff() / dt;
            const double rowRate = (r.tau - tauPrev).cwiseAbs().maxCoeff() / dt;
            rateCmd = std::max(rateCmd, rowRateCmd);
            rate = std::max(rate, rowRate);
            rateCmdRowErr = std::max(rateCmdRowErr, std::abs(r.tauRateCmd - rowRateCmd));
            rateRowErr = std::max(rateRowErr, std::abs(r.tauRate - rowRate));
        }
        CHECK(run.minJointMargin == margin);
        CHECK(maxDiff(run.peakDqFrac, dq) < 1e-15);
        CHECK(maxDiff(run.peakDdqFrac, ddq) < 1e-15);
        CHECK(maxDiff(run.peakTauFrac, tau) < 1e-15);
        CHECK_THAT(run.maxDqFrac, WithinAbs(dq.maxCoeff(), 1e-15));
        CHECK_THAT(run.maxDdqFrac, WithinAbs(ddq.maxCoeff(), 1e-15));
        CHECK_THAT(run.maxTauFrac, WithinAbs(tau.maxCoeff(), 1e-15));
        CHECK_THAT(run.maxFlangeSpeed, WithinAbs(speed, 1e-15));
        CHECK(run.maxXi == xi);
        CHECK(run.maxBeta == beta);
        CHECK(run.minSigma == sigma);
        CAPTURE(sci(rateCmdRowErr), sci(rateRowErr), rateCmd, rate);
        CHECK_THAT(run.maxTauRateCmd, WithinAbs(rateCmd, 1e-9));
        CHECK_THAT(run.maxTauRate, WithinAbs(rate, 1e-9));
        CHECK(rateCmdRowErr < 1e-9);
        CHECK(rateRowErr < 1e-9);
    }

    // The held torque of every logged period, recomputed with models of its own: tau = M qdd* + C qdot + G of the
    // model file at the logged (q, qdot), with G from the plant file instead if plantGravity, and no passive force
    // (the damping is not compensated). Also the logged G. Largest differences in N m.
    void recomputeHoldTorque(const MujocoRun& run, double& tauErr, double& gErr)
    {
        const MujocoRunOptions& o = run.options;
        MujocoPanda model(o.model.empty() ? o.plant : o.model);
        MujocoPanda plant(o.plant);
        Mat<7, 7> M = Mat<7, 7>::Zero();
        Vec7 bias = Vec7::Zero(), passive = Vec7::Zero();
        tauErr = gErr = 0.0;
        for (const MujocoLogRow& r : run.log) {
            model.setState(r.q, r.qd);
            model.dynamics(M, bias, passive);
            const Vec7 Gmodel = model.gravity(r.q);
            const Vec7 G = o.plantGravity ? plant.gravity(r.q) : Gmodel;
            const Vec7 tau = M * r.qdd + bias - Gmodel + G;  // eq:fb, with the gravity term exchanged if asked
            tauErr = std::max(tauErr, maxDiff(r.tau, tau));
            gErr = std::max(gErr, maxDiff(r.G, G));
        }
    }

    // How much of the joint damping acts on the arm, from the log of a hold run with the plant panda_nohand.xml:
    // over each period the plant's mean acceleration (qdot_{k+1} - qdot_k) / dt exceeds the commanded qdd*_k by
    // the damping's own acceleration M^-1 (-D qdot_k) (tau holds M qdd* + C qdot + G and nothing else), up to the
    // O(dt) change of M, C qdot + G and qdot within the period. Returns the least-squares factor c of that excess
    // on M^-1 (-D qdot): 1 if the damping acts and is not compensated, 0 if the plant has no damping or the torque
    // cancels it. Also the largest residual and the largest damping acceleration, rad/s^2.
    void dampingFelt(const MujocoRun& run, double& factor, double& residual, double& largest)
    {
        MujocoPanda::Options withDamping;
        withDamping.damping = true;
        MujocoPanda arm(run.options.plant, withDamping);
        Mat<7, 7> M = Mat<7, 7>::Zero();
        Vec7 bias = Vec7::Zero(), passive = Vec7::Zero();
        double num = 0.0, den = 0.0;
        residual = largest = 0.0;
        const double dt = run.options.dt;
        for (std::size_t k = 0; k + 1 < run.log.size(); ++k) {
            const MujocoLogRow& r = run.log[k];
            arm.setState(r.q, r.qd);
            arm.dynamics(M, bias, passive);  // passive = -D qdot
            const Vec7 damping = M.ldlt().solve(passive);
            const Vec7 excess = (run.log[k + 1].qd - r.qd) / dt - r.qdd;
            num += excess.dot(damping);
            den += damping.dot(damping);
            residual = std::max(residual, maxDiff(excess, damping));
            largest = std::max(largest, damping.cwiseAbs().maxCoeff());
        }
        factor = den > 0.0 ? num / den : 0.0;
    }
}  // namespace

// The hardware configuration (design doc, Hardware; HardwareBounds): tau from the start of each 1 ms period,
// held; mu ramp over 1 s; the model's joint damping on and not compensated. Measured over 2 s: joint margin 0.56 rad,
// |qd| 13 %, |qdd*| 21 % (at t = 0.67 s), |tau| 22 %, torque rate 1.3 N m/s, flange speed 0.065 m/s, |xi| 0.0312
// (the start), |beta| 0.014 rad, sigma_min(A) 0.047. Without the ramp the first command violates two margins.
TEST_CASE("simulateMujoco: 2 s hardware-like loop (hold, ramp, damping): every margin holds; without the ramp not")
{
    const PandaExample& ex = example();
    const MujocoRun run = simulateMujoco(ex, hardwareRun(2.0));
    REQUIRE(run.ok);
    CAPTURE(run.minJointMargin, run.maxDqFrac, run.maxDdqFrac, run.maxTauFrac, run.maxTauRateCmd,
            run.maxFlangeSpeed, run.maxXi, run.maxBeta, run.minSigma, run.finalEtaDot);
    CHECK(run.minJointMargin >= B::jointMargin);
    CHECK(run.maxDqFrac <= B::dqFrac);
    CHECK(run.maxDdqFrac <= B::ddqFrac);
    CHECK(run.maxTauFrac <= B::tauFrac);
    CHECK(run.maxTauRateCmd <= B::tauRateFrac * PandaArm::tauRateMax);
    CHECK(run.maxFlangeSpeed <= B::flangeSpeed);
    CHECK(run.maxXi <= B::xi);
    CHECK(run.maxBeta <= B::beta);
    CHECK(run.minSigma >= B::sigmaMin);
    CHECK(run.maxContacts == 0);
    CHECK(run.clampedSteps == 0);
    CHECK(run.warnings == 0);
    checkPeaks(run, 1e-3);

    // The commanded torque changes slowly with the ramp: 1.27 N m/s measured, far inside 600 N m/s. A rate taken
    // against a stale previous torque (or against tau - G of the first period) is hundreds of N m/s.
    CHECK(run.maxTauRateCmd < 2.0);

    // The held torque is eq:fb with the arm's own model and nothing for the damping
    double tauErr = 0.0, gErr = 0.0;
    recomputeHoldTorque(run, tauErr, gErr);
    CAPTURE(sci(tauErr), sci(gErr));
    CHECK(tauErr < 1e-12);
    CHECK(gErr < 1e-12);

    // The damping acts on the plant and is not compensated, which is what makes this run hardware-like. Directly:
    // the plant accelerates by qdd* + M^-1 (-D qdot) (factor 1; 0 without damping or with it compensated) ...
    double factor = 0.0, residual = 0.0, largest = 0.0;
    dampingFelt(run, factor, residual, largest);
    CAPTURE(factor, sci(residual), sci(largest));
    CHECK_THAT(factor, WithinAbs(1.0, 0.01));
    CHECK(largest > 0.5);
    CHECK(residual < 0.05 * largest);
    // ... and in the outputs: beta and eta_dot leave their designs (without the damping: |beta| 2.6e-4 and
    // eta_dot 0.080 at 2 s; with it 1.39e-2 and 0.053)
    CHECK(run.maxBeta > 5e-3);
    CHECK(run.finalEtaDot < 0.07);

    // The ramp: the first command at rest is exactly zero, so the torque starts at G and the arm does not move
    const MujocoLogRow& first = run.log.front();
    CHECK(first.ramp == 0.0);
    CHECK(first.qdd.isZero(0.0));
    CHECK(first.tauRateCmd < 1e-9);

    // Without the ramp: 49 % of the acceleration limit at t = 0, and a torque-rate jump above 600 N m/s
    MujocoRunOptions o = hardwareRun(0.005);
    o.rampTime = 0.0;
    const MujocoRun noRamp = simulateMujoco(ex, o);
    REQUIRE(noRamp.ok);
    CAPTURE(noRamp.maxDdqFrac, noRamp.maxTauRateCmd);
    CHECK(noRamp.maxDdqFrac > B::ddqFrac);
    CHECK(noRamp.maxTauRateCmd > B::tauRateFrac * PandaArm::tauRateMax);
    checkPeaks(noRamp, 1e-3);
}

// Runs D and D-G of the app, first second: the plant carries the Franka Hand (panda.xml, fingers held open), and M
// and C qdot + G come from a second MujocoPanda on panda_nohand.xml, set to the plant's arm state every period.
// Without the hand in G (D) the arm sags while the ramp is still small and |xi| leaves its 0.04 margin (0.21 at
// 0.72 s); with the hand's weight in G (D-G, as libfranka adds it when the hand is configured) every margin holds.
TEST_CASE("simulateMujoco: 1 s with the hand as payload: D sags out of the margins, D-G (hand in G) holds them")
{
    const PandaExample& ex = example();
    MujocoRunOptions o = hardwareRun(1.0);
    o.plant = withHand();
    o.model = noHand();
    const MujocoRun D = simulateMujoco(ex, o);
    o.plantGravity = true;
    const MujocoRun DG = simulateMujoco(ex, o);
    REQUIRE(D.ok);
    REQUIRE(DG.ok);
    CAPTURE(D.maxXi, DG.maxXi, DG.maxBeta, DG.maxDqFrac, DG.maxDdqFrac, DG.maxTauFrac, DG.maxTauRateCmd);

    CHECK(D.maxXi > 0.1);
    CHECK(DG.minJointMargin >= B::jointMargin);
    CHECK(DG.maxDqFrac <= B::dqFrac);
    CHECK(DG.maxDdqFrac <= B::ddqFrac);
    CHECK(DG.maxTauFrac <= B::tauFrac);
    CHECK(DG.maxTauRateCmd < 2.0);
    CHECK(DG.maxFlangeSpeed <= B::flangeSpeed);
    CHECK(DG.maxXi <= B::xi);
    CHECK(DG.maxBeta <= B::beta);
    CHECK(DG.maxContacts == 0);
    CHECK(DG.maxConstraintRows >= 1);  // the finger coupling of panda.xml
    CHECK(DG.clampedSteps == 0);
    checkPeaks(DG, 1e-3);

    // The torque of both runs, recomputed: arm-only M and C qdot + G at the plant's state, and in D-G the gravity
    // term of the plant (with the hand) in place of the arm's
    for (const MujocoRun* run : {&D, &DG}) {
        CAPTURE(run->options.plantGravity);
        double tauErr = 0.0, gErr = 0.0;
        recomputeHoldTorque(*run, tauErr, gErr);
        CAPTURE(sci(tauErr), sci(gErr));
        CHECK(tauErr < 1e-12);
        CHECK(gErr < 1e-12);
    }
    // The hand's weight in G is felt at the joints: up to several N m (the logged G of the two runs at the start)
    CHECK(maxDiff(DG.log.front().G, D.log.front().G) > 1.0);
}

TEST_CASE("simulateMujoco: the self-motion run: beta follows its design, xi stays near 0", SLOW_IN_DEBUG)
{
    // Run A2 of the app (one loop): beta from 0.2 along beta0 (1 + 4t) e^{-4t}, xi near 0, eta_dot = nu_d.
    const PandaExample& ex = example();
    MujocoRunOptions o = exactRun(0.0);
    o.start = Start::SelfMotion;
    o.crossCheck = false;
    const MujocoRun run = simulateMujoco(ex, o);
    REQUIRE(run.ok);
    REQUIRE(run.steps == 9627);  // one loop: k dt < L / nu_d
    CAPTURE(sci(run.maxXi), sci(run.maxBetaDeviation), sci(run.finalLiftError));
    CHECK_THAT(run.log.front().beta, WithinAbs(0.2, 1e-10));
    CHECK_THAT(run.log.front().etaDot, WithinAbs(PandaExample::nuD, 1e-12));
    CHECK(run.maxXi < 5e-5);              // measured 1.9e-5 (joint-acceleration level 1.2e-6, + O(h))
    CHECK(run.maxBetaDeviation < 3e-4);   // measured 2.4e-4 (joint-acceleration level 1.51e-4, + O(h))
    CHECK(run.finalLiftError < 1e-4);     // measured 5.6e-5 (joint-acceleration level 9.9e-6, + O(h))
    CHECK(run.baselineMaxBetaDeviation < 2e-4);
    CHECK(run.baselineFinalLiftError2 < 2e-5);
}
