// Tests of the hardware experiments on the MuJoCo Panda (sim/mujoco_experiments.hpp, docs/HARDWARE.md): MuJoCo's
// dynamics as a DynamicsModel, the preflight of every experiment with realistic torques, and every experiment's
// closed loop at the torque level as the robot driver runs it (held torque, libfranka's filter, uncompensated joint
// damping), with a push and with an unmodeled hand. Part of tests_mujoco.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mujoco_experiments.hpp"
#include "mujoco_panda.hpp"
#include "se3pf/experiments.hpp"
#include "se3pf/preflight.hpp"
#include "panda_fixture.hpp"
#include "test_helpers.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace se3pf;
using testing::maxDiff;
using testing::sci;

namespace {
    std::string modelFile(const char* name) { return defaultPandaModelDir() + "/" + name; }
    std::string noHand() { return modelFile("panda_nohand.xml"); }
    std::string withHand() { return modelFile("panda.xml"); }
    const PandaExample& example() { return testing::pandaExample(); }

    const ExperimentSpec& spec(const char* name)
    {
        const ExperimentSpec* e = findExperiment(name);
        REQUIRE(e != nullptr);
        return *e;
    }

    void checkPeaksWithin(const RunPeaks& p, const SafetyLimits& l)
    {
        CHECK(p.minJointMargin >= l.jointMargin);
        CHECK(p.maxDqFrac <= l.dqFrac);
        CHECK(p.maxDdqFrac <= l.ddqFrac);
        CHECK(p.maxTauFrac <= l.tauFrac);
        CHECK(p.maxTauRate <= l.tauRateFrac * PandaArm::tauRateMax);
        CHECK(p.maxFlangeSpeed <= l.flangeSpeed);
        CHECK(p.maxXi <= l.xi);
        CHECK(p.maxBeta <= l.beta);
        CHECK(p.minSigma >= l.sigmaMin);
    }

    MujocoExperimentOptions armPlant()
    {
        MujocoExperimentOptions o;
        o.plant = noHand();
        return o;
    }
}  // namespace

TEST_CASE("MujocoDynamics: MuJoCo's M, C qd and G, for the arm and for the arm with the hand")
{
    for (const std::string& file : {noHand(), withHand()}) {
        CAPTURE(file);
        const MujocoDynamics dyn(file);
        MujocoPanda arm(file);
        for (int k = 0; k < 20; ++k) {
            const Vec7 q = testing::randomConfig(0.2);
            const Vec7 dq = testing::randn<Vec7>(0.5);
            arm.setState(q, dq);
            Mat<7, 7> M = Mat<7, 7>::Zero();
            Vec7 bias = Vec7::Zero(), passive = Vec7::Zero();
            arm.dynamics(M, bias, passive);
            const Vec7 G = arm.gravity(q);
            CHECK(maxDiff(dyn.mass(q), M) < 1e-12);
            CHECK(maxDiff(dyn.gravity(q), G) < 1e-12);
            CHECK(maxDiff(Vec7(dyn.coriolis(q, dq) + dyn.gravity(q)), bias) < 1e-12);
            CHECK(dyn.coriolis(q, Vec7::Zero()).cwiseAbs().maxCoeff() < 1e-12);
            // C qd is quadratic in qd
            CHECK(maxDiff(dyn.coriolis(q, Vec7(2.0 * dq)), Vec7(4.0 * dyn.coriolis(q, dq))) < 1e-10);
            const Mat<7, 7> Md = dyn.mass(q);
            CHECK(maxDiff(Md, Mat<7, 7>(Md.transpose())) < 1e-12);
            CHECK(Eigen::SelfAdjointEigenSolver<Mat<7, 7>>(Md).eigenvalues().minCoeff() >= 0.1 - 1e-9);  // armature
        }
    }
    // The hand adds weight: G differs from the arm's at the rest posture (joint 2 and 4 carry it)
    const Vec7 qr = PandaArm::restPosture();
    CHECK(std::abs(MujocoDynamics(withHand()).gravity(qr)[1] - MujocoDynamics(noHand()).gravity(qr)[1]) > 1.0);
}

TEST_CASE("preflight with MuJoCo's model: every experiment passes, for the arm and for the arm with the hand")
{
    const MujocoDynamics arm(noHand()), hand(withHand());
    for (const ExperimentSpec& e : experiments()) {
        CAPTURE(e.name);
        const ExperimentSetup s(example(), e);
        for (const DynamicsModel* dyn : {static_cast<const DynamicsModel*>(&arm), static_cast<const DynamicsModel*>(&hand)}) {
            const PreflightReport r = preflight(s, dyn, s.monitor());
            INFO(r.summary());
            CHECK(r.ok);
            CHECK(r.torqueChecked);
            checkPeaksWithin(r.peaks, r.limits);
            // Realistic torques: gravity alone takes about a fifth of joint 2's or joint 4's limit
            CHECK(r.peaks.maxTauFrac > 0.2);
            CHECK(r.peaks.maxTauFrac < 0.35);
            // The mu ramp keeps the commanded torque rate tiny (MuJoCo run C: 1.3 N m/s)
            CHECK(r.peaks.maxTauRate < 2.0);
        }
    }
}

TEST_CASE("The experiments on the MuJoCo Panda (torque level, hold, libfranka's filter, damping): all end at rest "
          "within every margin",
          SLOW_IN_DEBUG)
{
    for (const ExperimentSpec& e : experiments()) {
        CAPTURE(e.name);
        const ExperimentSetup s(example(), e);
        const MujocoExperimentRun run = simulateExperimentMujoco(s, s.monitor(), armPlant());
        CAPTURE(run.tEnd, run.steadyXi, run.steadyBeta, run.cruiseEtaDot, run.peaks.maxBeta, run.peaks.maxDdqFrac);
        CHECK(run.finished);
        CHECK(run.violation == Violation::None);
        CHECK(run.warnings == 0);
        CHECK(run.maxContacts == 0);
        CHECK(run.clampedPeriods == 0);
        checkPeaksWithin(run.peaks, s.monitor().limits());
        CHECK(run.endSpeed < RunConstants::restSpeed);
        CHECK(run.tEnd >= e.duration);
        CHECK(run.tEnd < e.duration + 1.0);
        // Every row is a following row; the log matches the peaks
        CHECK(static_cast<long>(run.log.size()) == run.peaks.rows);
        // The damping is felt while moving (docs/MUJOCO.md, run C): the path speed lags, beta and xi leave the lift
        if (e.mode == TangentialReference::Mode::Velocity) {
            CHECK(run.cruiseEtaDot < 0.9 * e.speed.nuMax);
            CHECK(run.cruiseEtaDot > 0.6 * e.speed.nuMax);
            CHECK(run.peaks.maxBeta > 5e-3);
            CHECK(run.peaks.maxBeta < 0.6 * HardwareBounds::beta);
        }
    }
}

TEST_CASE("The paper experiment on the MuJoCo Panda with the Franka Hand configured (run D-G): every margin",
          SLOW_IN_DEBUG)
{
    const ExperimentSetup s(example(), spec("paper"));
    MujocoExperimentOptions o;
    o.plant = withHand();
    o.model = noHand();     // M and C qd without the hand ...
    o.plantGravity = true;  // ... G with it, as libfranka adds it with the load configured in Desk
    const MujocoExperimentRun run = simulateExperimentMujoco(s, s.monitor(), o);
    CAPTURE(run.peaks.maxTauFrac, run.steadyXi);
    CHECK(run.finished);
    CHECK(run.violation == Violation::None);
    checkPeaksWithin(run.peaks, s.monitor().limits());
    CHECK(run.peaks.maxTauFrac > 0.27);  // the hand's weight: 29.5 % on joint 2 (25.1 % without it)
}

TEST_CASE("A push on joint 4 during hold on the MuJoCo Panda: the safety stop engages and holds the arm at rest")
{
    const ExperimentSetup s(example(), spec("hold"));
    MujocoExperimentOptions o = armPlant();
    o.pushJoint = 3;
    o.pushTorque = 12.0;
    o.pushStart = 3.0;
    o.pushDuration = 0.3;
    const MujocoExperimentRun run = simulateExperimentMujoco(s, s.monitor(), o);
    CAPTURE(run.tStop, run.tEnd, run.endSpeed, run.maxStopDrift);
    CHECK_FALSE(run.finished);
    CHECK(run.violation != Violation::None);
    CHECK(run.tStop > 3.0);
    CHECK(run.tStop < 3.1);
    CHECK(run.stopSettled);
    CHECK(run.endSpeed < RunConstants::stopSpeed);
    CHECK(run.tEnd - run.tStop < RunConstants::stopTime);
    // The push (12 N m against K_s = 300 N m/rad) deflects the hold while it lasts; the hold brings the arm back
    CHECK(run.maxStopDrift < 0.1);
    // The stop's law after the violation: the damped hold at the configuration of the violation
    const ExperimentLogRow* first = nullptr;
    long stopRows = 0;
    for (const ExperimentLogRow& row : run.log)
        if (row.phase == 1) {
            if (first == nullptr)
                first = &row;
            ++stopRows;
            CHECK(maxDiff(row.tau, SafetyStop::torque(row.q, row.dq, first->q)) < 1e-12);
        }
    CHECK(stopRows > 10);
    REQUIRE(first != nullptr);
    CHECK(maxDiff(run.log.back().q, first->q) < 5e-3);  // back where it was stopped
    CHECK(run.warnings == 0);
}

TEST_CASE("An unmodeled hand on the MuJoCo Panda (run D): the monitor stops the sagging arm and the stop holds it")
{
    // The hand (0.76 kg in the Menagerie model) in the plant, not in the robot's G: with joint damping alone the
    // flange fell 0.3 m in the 2 s of the stop; the damped hold keeps every joint within 0.03 rad of where it stopped.
    for (const char* name : {"hold", "paper", "self-motion"}) {
        CAPTURE(name);
        const ExperimentSetup s(example(), spec(name));
        MujocoExperimentOptions o;
        o.plant = withHand();
        o.model = noHand();
        o.plantGravity = false;  // the hand missing in Desk: G without it
        const MujocoExperimentRun run = simulateExperimentMujoco(s, s.monitor(), o);
        CAPTURE(run.tStop, toString(run.violation), run.maxStopDrift, run.endSpeed, run.tEnd);
        CHECK_FALSE(run.finished);
        CHECK(run.violation != Violation::None);
        CHECK(run.tStop < 0.3);
        double maxXi = 0.0;
        for (const ExperimentLogRow& row : run.log)
            if (row.phase == 0)
                maxXi = std::max(maxXi, row.xiNorm);
        CHECK(maxXi <= HardwareBounds::xi);  // stopped before |xi| left its bound (run D without a monitor: 0.215)
        CHECK(run.maxStopDrift < 0.03);      // MuJoCo: 0.021 (hold), 0.022 (paper)
        CHECK(run.stopSettled);
        CHECK(run.endSpeed < RunConstants::stopSpeed);
        CHECK(run.tEnd - run.tStop < 1.5);
        CHECK(run.warnings == 0);
    }
}

TEST_CASE("The safety stop from 50 % of every joint's velocity limit on the MuJoCo Panda stops within the joint "
          "margin and settles")
{
    // The worst case of a JointLimit stop: the monitor allows |qd_i| up to 50 % of the limit, and stops 0.10 rad
    // from a limit. Damping alone (D = 10) travelled 0.17-0.23 rad here; the damped hold must stay within the margin.
    for (const ExperimentSpec& e : experiments()) {
        const ExperimentSetup s(example(), e);
        for (int sign : {1, -1}) {
            Vec7 qd0 = 0.5 * PandaArm::dqMax();
            for (int i = 0; i < 7; ++i)
                if ((i % 2 == 0) != (sign > 0))
                    qd0[i] = -qd0[i];
            for (const std::string& plant : {noHand(), withHand()}) {
                CAPTURE(e.name, sign, plant);
                const MujocoStopRun r = simulateStopMujoco(plant, "", s.start(), qd0);
                CAPTURE(r.duration, r.endSpeed, r.drift.transpose(), r.maxTorque);
                CHECK(r.drift.maxCoeff() < HardwareBounds::jointMargin);  // MuJoCo: at most 0.07 rad
                CHECK(r.settled);
                CHECK(r.duration < 1.5);
                CHECK(r.maxTorque < 45.0);  // MuJoCo: 39 N m, below the clamp (52 N m on joints 1-4)
                CHECK(r.minJointMargin > 0.0);
                CHECK(r.warnings == 0);
            }
        }
    }
    // The unmodeled hand, from rest: the hold deflects by less than 0.02 rad
    const MujocoStopRun hand = simulateStopMujoco(withHand(), noHand(), example().start(example().redundancy(),
                                                                                        Vec5::Zero(), 0.0),
                                                  Vec7::Zero());
    CAPTURE(hand.drift.transpose());
    CHECK(hand.drift.maxCoeff() < 0.02);
    CHECK(hand.settled);
}
