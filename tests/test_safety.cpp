// Tests of the online safety monitor and the approach motion (se3pf/safety.hpp; design doc, Hardware;
// docs/HARDWARE.md): every check, every joint against its own limit, the boundaries (a value at its bound passes, one
// beyond it fails), NaN and infinity, and that no experiment can loosen a margin of HardwareBounds.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "se3pf/panda.hpp"
#include "se3pf/safety.hpp"
#include "test_helpers.hpp"

#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>

using namespace se3pf;
using Catch::Matchers::WithinAbs;
using testing::maxDiff;
using testing::sci;
using B = HardwareBounds;

namespace {
    const double kNaN = std::numeric_limits<double>::quiet_NaN();
    const double kInf = std::numeric_limits<double>::infinity();

    // A configuration well inside every joint range (smallest margin 0.79 rad, joint 6)
    Vec7 restQ() { return PandaArm::restPosture(); }

    Vec7 unit(int i, double value)
    {
        Vec7 v = Vec7::Zero();
        v[i] = value;
        return v;
    }

    // Every field of a SafetyLimits, to loop over them
    double* field(SafetyLimits& l, int k)
    {
        double* f[] = {&l.jointMargin, &l.dqFrac, &l.ddqFrac, &l.tauFrac, &l.tauRateFrac,
                       &l.flangeSpeed, &l.xi, &l.beta, &l.sigmaMin};
        return f[k];
    }
    double value(const SafetyLimits& l, int k)
    {
        SafetyLimits copy = l;
        return *field(copy, k);
    }
    constexpr int kFields = 9;
    // jointMargin and sigmaMin are lower bounds ("at least"); the others are upper bounds
    bool isLowerBound(int k) { return k == 0 || k == 8; }

    // Flange speed |v| of the twist J_b(q) qd
    double flangeSpeed(const Vec7& q, const Vec7& dq)
    {
        Mat67 Jb = Mat67::Zero();
        PandaArm::output(q, &Jb);
        return Vec6(Jb * dq).tail<3>().norm();
    }
}  // namespace

TEST_CASE("Violation: every value has its own description")
{
    std::set<std::string> seen;
    const Violation all[] = {Violation::None,         Violation::Projection,  Violation::JointLimit,
                             Violation::JointVelocity, Violation::Acceleration, Violation::Torque,
                             Violation::TorqueRate,   Violation::FlangeSpeed, Violation::Tube,
                             Violation::Redundancy,   Violation::Conditioning, Violation::NonFinite,
                             Violation::Timeout};
    for (Violation v : all) {
        const std::string s = toString(v);
        CAPTURE(s);
        CHECK(!s.empty());
        CHECK(s != "unknown");
        seen.insert(s);
    }
    CHECK(seen.size() == sizeof(all) / sizeof(all[0]));
}

TEST_CASE("SafetyMonitor: the default limits are HardwareBounds")
{
    const SafetyMonitor m;
    const SafetyLimits& l = m.limits();
    CHECK(l.jointMargin == B::jointMargin);
    CHECK(l.dqFrac == B::dqFrac);
    CHECK(l.ddqFrac == B::ddqFrac);
    CHECK(l.tauFrac == B::tauFrac);
    CHECK(l.tauRateFrac == B::tauRateFrac);
    CHECK(l.flangeSpeed == B::flangeSpeed);
    CHECK(l.xi == B::xi);
    CHECK(l.beta == B::beta);
    CHECK(l.sigmaMin == B::sigmaMin);
    CHECK_FALSE(m.clamped());
    // The design doc's values, so that a change of HardwareBounds itself shows up here
    CHECK(B::jointMargin == 0.10);
    CHECK(B::dqFrac == 0.5);
    CHECK(B::ddqFrac == 0.3);
    CHECK(B::tauFrac == 0.6);
    CHECK(B::tauRateFrac == 0.6);
    CHECK(B::flangeSpeed == 0.5);
    CHECK(B::xi == 0.04);
    CHECK(B::beta == 0.1);
    CHECK(B::sigmaMin == 1e-3);
}

TEST_CASE("SafetyMonitor: a margin can not be loosened; a tighter value is kept; invalid values throw")
{
    const SafetyLimits hw;
    for (int k = 0; k < kFields; ++k) {
        CAPTURE(k);
        const double bound = value(hw, k);
        // Looser: a smaller lower bound or a larger upper bound is clamped to the hardware bound
        SafetyLimits looser;
        *field(looser, k) = isLowerBound(k) ? 0.0 : 10.0 * bound;
        const SafetyMonitor ml(looser);
        CHECK(ml.clamped());
        CHECK(value(ml.limits(), k) == bound);
        // Tighter: kept as asked
        SafetyLimits tighter;
        const double t = isLowerBound(k) ? 2.0 * bound : 0.5 * bound;
        *field(tighter, k) = t;
        const SafetyMonitor mt(tighter);
        CHECK_FALSE(mt.clamped());
        CHECK(value(mt.limits(), k) == t);
        // NaN, infinity and negative values are refused
        for (double bad : {kNaN, kInf, -kInf, -1e-9}) {
            SafetyLimits b;
            *field(b, k) = bad;
            CHECK_THROWS_AS(SafetyMonitor(b), std::invalid_argument);
        }
    }

    // A loosened monitor really checks the hardware bound: |qd_1| = 0.6 dqMax fails with dqFrac = 0.9 requested.
    SafetyLimits loose;
    loose.dqFrac = 0.9;
    loose.ddqFrac = 0.9;
    loose.xi = 1.0;
    const SafetyMonitor m(loose);
    CHECK(m.checkState(restQ(), unit(6, 0.6 * PandaArm::dqMax()[6])) == Violation::JointVelocity);
    CHECK(m.checkAcceleration(unit(0, 0.31 * PandaArm::ddqMax()[0])) == Violation::Acceleration);
    CHECK(m.checkOutputs((Vec5() << 0.05, 0, 0, 0, 0).finished(), 0.0, 1.0) == Violation::Tube);

    // A tighter monitor checks the tighter value
    SafetyLimits tight;
    tight.dqFrac = 0.2;
    tight.xi = 0.01;
    const SafetyMonitor mt(tight);
    CHECK(mt.checkState(restQ(), unit(0, 0.25 * PandaArm::dqMax()[0])) == Violation::JointVelocity);
    CHECK(mt.checkOutputs((Vec5() << 0.011, 0, 0, 0, 0).finished(), 0.0, 1.0) == Violation::Tube);
}

TEST_CASE("SafetyMonitor::checkOutputs: |xi|, |beta| and sigma_min(A) at and beyond their bounds; NaN")
{
    const SafetyMonitor m;
    const Vec5 zero = Vec5::Zero();
    CHECK(m.checkOutputs(zero, 0.0, 1.0) == Violation::None);

    // |xi| = 0.04 exactly passes (along each axis), just beyond fails
    for (int i = 0; i < 5; ++i) {
        CAPTURE(i);
        Vec5 xi = Vec5::Zero();
        xi[i] = -B::xi;
        CHECK(m.checkOutputs(xi, 0.0, 1.0) == Violation::None);
        xi[i] = B::xi * (1.0 + 1e-12);
        CHECK(m.checkOutputs(xi, 0.0, 1.0) == Violation::Tube);
    }
    // The norm is checked, not the entries: (0.03, 0.03) has |xi| = 0.0424
    CHECK(m.checkOutputs((Vec5() << 0.03, 0.03, 0, 0, 0).finished(), 0.0, 1.0) == Violation::Tube);

    CHECK(m.checkOutputs(zero, B::beta, 1.0) == Violation::None);
    CHECK(m.checkOutputs(zero, -B::beta, 1.0) == Violation::None);
    CHECK(m.checkOutputs(zero, B::beta * (1.0 + 1e-12), 1.0) == Violation::Redundancy);
    CHECK(m.checkOutputs(zero, -B::beta * (1.0 + 1e-12), 1.0) == Violation::Redundancy);

    CHECK(m.checkOutputs(zero, 0.0, B::sigmaMin) == Violation::None);
    CHECK(m.checkOutputs(zero, 0.0, B::sigmaMin * (1.0 - 1e-12)) == Violation::Conditioning);
    CHECK(m.checkOutputs(zero, 0.0, 0.0) == Violation::Conditioning);

    // Non-finite values in any argument
    for (double bad : {kNaN, kInf}) {
        Vec5 xi = Vec5::Zero();
        xi[3] = bad;
        CHECK(m.checkOutputs(xi, 0.0, 1.0) == Violation::NonFinite);
        CHECK(m.checkOutputs(zero, bad, 1.0) == Violation::NonFinite);
        CHECK(m.checkOutputs(zero, 0.0, bad) == Violation::NonFinite);
    }
}

TEST_CASE("SafetyMonitor::checkState: every joint at its position margin on both sides; NaN")
{
    const SafetyMonitor m;
    const Vec7 lo = PandaArm::qMin(), hi = PandaArm::qMax();
    REQUIRE(m.checkState(restQ(), Vec7::Zero()) == Violation::None);
    for (int i = 0; i < 7; ++i) {
        CAPTURE(i);
        Vec7 q = restQ();
        q[i] = lo[i] + B::jointMargin + 1e-9;
        CHECK(m.checkState(q, Vec7::Zero()) == Violation::None);
        q[i] = lo[i] + B::jointMargin - 1e-9;
        CHECK(m.checkState(q, Vec7::Zero()) == Violation::JointLimit);
        q[i] = hi[i] - B::jointMargin - 1e-9;
        CHECK(m.checkState(q, Vec7::Zero()) == Violation::None);
        q[i] = hi[i] - B::jointMargin + 1e-9;
        CHECK(m.checkState(q, Vec7::Zero()) == Violation::JointLimit);
        q[i] = hi[i] + 1.0;  // outside the range
        CHECK(m.checkState(q, Vec7::Zero()) == Violation::JointLimit);
    }
    for (double bad : {kNaN, kInf, -kInf}) {
        for (int i = 0; i < 7; ++i) {
            Vec7 q = restQ();
            q[i] = bad;
            CHECK(m.checkState(q, Vec7::Zero()) == Violation::NonFinite);
            CHECK(m.checkState(restQ(), unit(i, bad)) == Violation::NonFinite);
        }
    }
}

TEST_CASE("SafetyMonitor::checkState: every joint against its own velocity limit; the flange speed")
{
    const SafetyMonitor m;
    const Vec7 dqMax = PandaArm::dqMax();
    for (int i = 0; i < 7; ++i) {
        CAPTURE(i);
        // At 50 % of its own limit the velocity check passes (the flange speed may still stop it); beyond, it fails
        for (double sign : {1.0, -1.0}) {
            CHECK(m.checkState(restQ(), unit(i, sign * B::dqFrac * dqMax[i])) != Violation::JointVelocity);
            CHECK(m.checkState(restQ(), unit(i, sign * B::dqFrac * dqMax[i] * (1.0 + 1e-12)))
                  == Violation::JointVelocity);
        }
    }
    // Joints 5-7 have 2.61 rad/s, joints 1-4 2.175: 1.2 rad/s is above 50 % for joint 4, not for joint 5
    CHECK(m.checkState(restQ(), unit(3, 1.2)) == Violation::JointVelocity);
    CHECK(m.checkState(restQ(), unit(4, 1.2)) != Violation::JointVelocity);

    // Flange speed: scale a joint velocity direction so that |v| is 0.5 m/s, with every joint below 50 %
    const Vec7 u = (Vec7() << 0.5, 1.0, 0.5, -1.0, 0.3, 0.5, 0.2).finished();
    const double perUnit = flangeSpeed(restQ(), u);
    const double s = B::flangeSpeed / perUnit;
    REQUIRE((s * u).cwiseAbs().cwiseQuotient(dqMax).maxCoeff() < B::dqFrac);  // the velocity check passes
    CAPTURE(perUnit, s);
    CHECK(m.checkState(restQ(), (s * (1.0 - 1e-9)) * u) == Violation::None);
    CHECK(m.checkState(restQ(), (s * (1.0 + 1e-9)) * u) == Violation::FlangeSpeed);
    // The angular part does not count: joint 7 turns the flange about its own axis without moving its origin
    CHECK(flangeSpeed(restQ(), unit(6, 1.0)) < 1e-12);
    CHECK(m.checkState(restQ(), unit(6, 1.3)) == Violation::None);
}

TEST_CASE("SafetyMonitor::checkAcceleration: every joint at 30 % of its own limit; NaN")
{
    const SafetyMonitor m;
    const Vec7 ddqMax = PandaArm::ddqMax();
    CHECK(m.checkAcceleration(Vec7::Zero()) == Violation::None);
    for (int i = 0; i < 7; ++i) {
        CAPTURE(i);
        for (double sign : {1.0, -1.0}) {
            CHECK(m.checkAcceleration(unit(i, sign * B::ddqFrac * ddqMax[i])) == Violation::None);
            CHECK(m.checkAcceleration(unit(i, sign * B::ddqFrac * ddqMax[i] * (1.0 + 1e-12)))
                  == Violation::Acceleration);
        }
        CHECK(m.checkAcceleration(unit(i, kNaN)) == Violation::NonFinite);
        CHECK(m.checkAcceleration(unit(i, -kInf)) == Violation::NonFinite);
    }
    // Joint 2 has 7.5 rad/s^2, joint 6 has 20: 3 rad/s^2 is 40 % for joint 2 and 15 % for joint 6
    CHECK(m.checkAcceleration(unit(1, 3.0)) == Violation::Acceleration);
    CHECK(m.checkAcceleration(unit(5, 3.0)) == Violation::None);
}

TEST_CASE("SafetyMonitor::checkTorque: delivered torque per joint, commanded torque rate, dt; NaN")
{
    const SafetyMonitor m;
    const Vec7 tauMax = PandaArm::tauMax();
    const double dt = 1e-3;
    const Vec7 zero = Vec7::Zero();
    CHECK(m.checkTorque(zero, zero, zero, dt) == Violation::None);

    for (int i = 0; i < 7; ++i) {
        CAPTURE(i);
        for (double sign : {1.0, -1.0}) {
            // Delivered torque at 60 % of its own limit (52.2 N m for joints 1-4, 7.2 N m for 5-7)
            const Vec7 at = unit(i, sign * B::tauFrac * tauMax[i]);
            CHECK(m.checkTorque(at, zero, zero, dt) == Violation::None);
            CHECK(m.checkTorque(Vec7(at * (1.0 + 1e-12)), zero, zero, dt) == Violation::Torque);
            // Commanded rate: 600 N m/s is 0.6 N m in one period, measured against the previous command
            const Vec7 prev = unit(i, 3.0);
            const Vec7 rateAt = prev + unit(i, sign * B::tauRateFrac * PandaArm::tauRateMax * dt * 0.999);
            CHECK(m.checkTorque(zero, rateAt, prev, dt) == Violation::None);
            CHECK(m.checkTorque(zero, unit(i, sign * B::tauRateFrac * PandaArm::tauRateMax * dt), zero, dt)
                  == Violation::None);
            const Vec7 rateBeyond = prev + unit(i, sign * B::tauRateFrac * PandaArm::tauRateMax * dt * 1.001);
            CHECK(m.checkTorque(zero, rateBeyond, prev, dt) == Violation::TorqueRate);
            // The same command after a period of 2 ms is half the rate
            CHECK(m.checkTorque(zero, rateBeyond, prev, 2.0 * dt) == Violation::None);
        }
        CHECK(m.checkTorque(unit(i, kNaN), zero, zero, dt) == Violation::NonFinite);
        CHECK(m.checkTorque(zero, unit(i, kInf), zero, dt) == Violation::NonFinite);
        CHECK(m.checkTorque(zero, zero, unit(i, kNaN), dt) == Violation::NonFinite);
    }
    // 10 N m is fine for joint 1 (87 N m) and too much for joint 5 (12 N m)
    CHECK(m.checkTorque(unit(0, 10.0), zero, zero, dt) == Violation::None);
    CHECK(m.checkTorque(unit(4, 10.0), zero, zero, dt) == Violation::Torque);
    // The torque margin is on the delivered torque (tau + G), the rate on the commanded one (tau): a large G with a
    // steady command passes the rate, a jump of the command fails it
    CHECK(m.checkTorque(unit(1, 40.0), unit(1, 0.1), unit(1, 0.1), dt) == Violation::None);
    CHECK(m.checkTorque(unit(1, 40.0), unit(1, 1.0), zero, dt) == Violation::TorqueRate);
    // A period that is not positive or not finite fails closed
    for (double badDt : {0.0, -1e-3, kNaN, kInf})
        CHECK(m.checkTorque(zero, zero, zero, badDt) == Violation::NonFinite);
}

TEST_CASE("ApproachMotion: within 20 % / 10 % / 5 % of the velocity, acceleration and jerk limits, at least 3 s")
{
    const Vec7 q0 = PandaArm::restPosture();
    const Vec7 jerkMax = PandaArm::jerkMax(), dqMax = PandaArm::dqMax(), ddqMax = PandaArm::ddqMax();
    const Vec7 moves[] = {(Vec7() << 0.6, -0.3, 0.8, 0.4, -0.5, 0.3, 1.0).finished(),
                          (Vec7() << 0.01, 0, 0, 0, 0, 0, 0).finished(), Vec7::Zero(),
                          (Vec7() << 0, 1.2, 0, 0, 0, 0, 0).finished()};
    for (const Vec7& d : moves) {
        const ApproachMotion m(q0, q0 + d);
        const double T = m.duration();
        CAPTURE(T);
        CHECK(T >= ApproachMotion::kMinTime);
        double v = 0.0, a = 0.0, j = 0.0;
        for (double t = 0.0; t <= T; t += 1e-3) {
            v = std::max(v, m.velocity(t).cwiseAbs().cwiseQuotient(dqMax).maxCoeff());
            a = std::max(a, m.acceleration(t).cwiseAbs().cwiseQuotient(ddqMax).maxCoeff());
            j = std::max(j, m.jerk(t).cwiseAbs().cwiseQuotient(jerkMax).maxCoeff());
        }
        CHECK(v <= ApproachMotion::kVelocityFraction + 1e-12);
        CHECK(a <= ApproachMotion::kAccelerationFraction + 1e-12);
        CHECK(j <= ApproachMotion::kJerkFraction + 1e-12);
        // Starts and ends at rest, ends exactly at the goal (and stays there)
        CHECK(maxDiff(m.position(0.0), q0) == 0.0);
        CHECK(maxDiff(m.position(T), Vec7(q0 + d)) == 0.0);
        CHECK(maxDiff(m.position(T + 5.0), Vec7(q0 + d)) == 0.0);
        CHECK(m.velocity(0.0).cwiseAbs().maxCoeff() == 0.0);
        CHECK(m.velocity(T).cwiseAbs().maxCoeff() < 1e-15);
        CHECK(m.acceleration(0.0).cwiseAbs().maxCoeff() == 0.0);
        CHECK(m.acceleration(T).cwiseAbs().maxCoeff() < 1e-15);
        CHECK(m.jerk(T + 1e-3).cwiseAbs().maxCoeff() == 0.0);
        CHECK(m.maxJointDistance() == d.cwiseAbs().maxCoeff());
        // Derivatives against central differences in the middle of the motion
        const double t = 0.37 * T, h = 1e-5;
        CHECK(maxDiff(Vec7((m.position(t + h) - m.position(t - h)) / (2 * h)), m.velocity(t)) < 1e-8);
        CHECK(maxDiff(Vec7((m.velocity(t + h) - m.velocity(t - h)) / (2 * h)), m.acceleration(t)) < 1e-7);
        CHECK(maxDiff(Vec7((m.acceleration(t + h) - m.acceleration(t - h)) / (2 * h)), m.jerk(t)) < 1e-6);
    }
    // The velocity bound sets the time of a long move of joint 2: 15/8 * 1.2 / (0.2 * 2.175) = 5.17 s
    CHECK_THAT(ApproachMotion(q0, q0 + moves[3]).duration(), WithinAbs(1.875 * 1.2 / (0.2 * 2.175), 1e-12));
    CHECK_THROWS_AS(ApproachMotion(q0, unit(2, kNaN)), std::invalid_argument);
}
