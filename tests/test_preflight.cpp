// Tests of the hardware experiments and their preflight (se3pf/experiments.hpp, se3pf/preflight.hpp;
// docs/HARDWARE.md): the catalog, the speed profile, the setup (start configuration, controller), the shared control
// cycle ExperimentLoop, the friction feedforward (CTRL-1), and the preflight: every experiment passes, kinematic
// only and with a dynamics model, and a bad experiment fails with the right reason.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "se3pf/experiments.hpp"
#include "se3pf/preflight.hpp"
#include "se3pf/safety.hpp"
#include "panda_fixture.hpp"
#include "test_helpers.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace se3pf;
using Catch::Matchers::WithinAbs;
using testing::maxDiff;
using testing::pandaExample;
using testing::sci;
using B = HardwareBounds;

namespace {
    // A simple stand-in for franka::Model: a constant diagonal inertia plus a payload of 1 kg at the flange
    // (M = diag(m) + J_v^T J_v), no Coriolis terms, and the payload's weight as G. Smooth, positive definite, and with
    // torques of the Panda's order (joint 2 carries about 6 N m of G at the rest posture).
    class FakeDynamics final : public DynamicsModel {
        public:
            Mat<7, 7> mass(const Vec7& q) const override
            {
                Mat67 Jb = Mat67::Zero();
                PandaArm::output(q, &Jb);
                const Mat<3, 7> Jv = Jb.bottomRows<3>();
                const Vec7 rotor = (Vec7() << 0.6, 1.2, 0.4, 0.6, 0.05, 0.05, 0.02).finished();
                return Mat<7, 7>(rotor.asDiagonal()) + Jv.transpose() * Jv;
            }
            Vec7 coriolis(const Vec7&, const Vec7&) const override { return Vec7::Zero(); }
            Vec7 gravity(const Vec7& q) const override
            {
                // G = J_v^T (R^T m g e_z): the weight of 1 kg at the flange origin, in the flange frame
                Mat67 Jb = Mat67::Zero();
                const Pose g = PandaArm::output(q, &Jb);
                const Vector3d fBody = g.block<3, 3>(0, 0).transpose() * Vector3d(0.0, 0.0, 9.81);
                return Jb.bottomRows<3>().transpose() * fBody;
            }
    };

    // An inertia-only model: M = m I, no Coriolis, no gravity
    class ScaledInertia final : public DynamicsModel {
        public:
            explicit ScaledInertia(double m) : m_(m) {}
            Mat<7, 7> mass(const Vec7&) const override { return m_ * Mat<7, 7>::Identity(); }
            Vec7 coriolis(const Vec7&, const Vec7&) const override { return Vec7::Zero(); }
            Vec7 gravity(const Vec7&) const override { return Vec7::Zero(); }

        private:
            double m_;
    };

    const ExperimentSpec& spec(const char* name)
    {
        const ExperimentSpec* e = findExperiment(name);
        REQUIRE(e != nullptr);
        return *e;
    }

    // The peaks of a report must respect the limits it was checked against
    void checkWithinLimits(const PreflightReport& r)
    {
        const RunPeaks& p = r.peaks;
        const SafetyLimits& l = r.limits;
        CHECK(p.minJointMargin >= l.jointMargin);
        CHECK(p.maxDqFrac <= l.dqFrac);
        CHECK(p.maxDdqFrac <= l.ddqFrac);
        CHECK(p.maxFlangeSpeed <= l.flangeSpeed);
        CHECK(p.maxXi <= l.xi);
        CHECK(p.maxBeta <= l.beta);
        CHECK(p.minSigma >= l.sigmaMin);
        if (r.torqueChecked) {
            CHECK(p.maxTauFrac <= l.tauFrac);
            CHECK(p.maxTauRate <= l.tauRateFrac * PandaArm::tauRateMax);
        }
    }
}  // namespace

// ---------------------------------------------------------------------------------------------------------------
// The catalog and the speed profile

TEST_CASE("experiments: the catalog, in the order the lab runs it; each spec is valid and never loosens a margin")
{
    const std::vector<ExperimentSpec>& all = experiments();
    REQUIRE(all.size() == 5);
    const char* order[] = {"hold", "loops", "recover", "paper", "self-motion"};
    for (int i = 0; i < 5; ++i) {
        CAPTURE(i);
        CHECK(all[i].name == order[i]);
        CHECK(findExperiment(order[i]) == &all[i]);
    }
    CHECK(findExperiment("no-such-experiment") == nullptr);
    CHECK(findExperiment("") == nullptr);

    for (const ExperimentSpec& e : all) {
        CAPTURE(e.name);
        CHECK_FALSE(e.title.empty());
        CHECK_FALSE(e.purpose.empty());
        CHECK_FALSE(e.watch.empty());
        // Never looser than the hardware bounds: the monitor need not clamp anything
        CHECK_FALSE(SafetyMonitor(e.limits).clamped());
        // Starts inside the tube and redundancy margins, with room to spare
        CHECK(e.xi0.norm() <= 0.8 * B::xi);
        CHECK(std::abs(e.beta0) <= 0.04);
        // The mu ramp (Rules for the 1 kHz loop: hardware only, T = 1-2 s)
        CHECK(e.rampTime >= 1.0);
        CHECK(e.rampTime <= 2.0);
        // The paper's poles at -4 and k_t = 5, except loops: poles at -6 and k_t = 10 (CTRL-1, 2026-10-19)
        const double lambda = e.name == "loops" ? 6.0 : 4.0;
        CHECK(e.gains.Kp == lambda * lambda);
        CHECK(e.gains.Kd == 2.0 * lambda);
        CHECK(e.gains.KpBeta == lambda * lambda);
        CHECK(e.gains.KdBeta == 2.0 * lambda);
        // Integral action only in loops (CTRL-1)
        if (e.name == "loops") {
            CHECK(e.gains.kiEta == 25.0);
            CHECK(e.gains.Ki == 36.0);
            CHECK(e.gains.KiBeta == 36.0);
            CHECK(e.gains.integralMaxEta == 1.0);
            CHECK(e.gains.integralMaxXi == 1.0);
            CHECK(e.gains.integralMaxBeta == 3.0);
        } else {
            CHECK(e.gains.kiEta == 0.0);
            CHECK(e.gains.Ki == 0.0);
            CHECK(e.gains.KiBeta == 0.0);
        }
        if (e.mode == TangentialReference::Mode::Velocity) {
            CHECK(e.gains.kt == (e.name == "loops" ? 15.0 : 5.0));
            CHECK(e.speed.end() <= e.duration);
            CHECK(e.speed.rampUp >= 2.0);  // smooth start
            CHECK(e.speed.rampDown >= 4.0);  // ends at rest
            CHECK(e.speed.nuMax <= PandaExample::nuD);
        }
    }

    // The specifics of each experiment (docs/HARDWARE.md)
    const PandaExample& ex = pandaExample();
    const double L = ex.tube().length();
    CHECK(spec("hold").mode == TangentialReference::Mode::Position);
    CHECK(spec("hold").duration == 10.0);
    CHECK(spec("hold").xi0.norm() == 0.0);
    CHECK(spec("loops").speed.nuMax == 0.05);
    CHECK(std::abs(spec("loops").speed.distance() - 2.0 * L) < 1e-3);  // two loops: ends where it starts
    CHECK(spec("loops").duration == 38.0);
    CHECK(maxDiff(spec("recover").xi0, Vec5(0.5 * PandaExample::xi0())) == 0.0);
    CHECK(maxDiff(spec("paper").xi0, PandaExample::xi0()) == 0.0);
    for (const char* name : {"recover", "paper"}) {
        CAPTURE(name);
        CHECK(spec(name).speed.nuMax == PandaExample::nuD);
        CHECK(spec(name).speed.rampUp == 2.0);
        CHECK(spec(name).speed.rampDown == 4.0);
        CHECK(spec(name).speed.end() == 20.0);  // ramps down over the last 4 s of 20 s
        CHECK(spec(name).duration == 20.0);
        CHECK(spec(name).beta0 == 0.0);
    }
    CHECK(spec("self-motion").mode == TangentialReference::Mode::Position);
    CHECK(spec("self-motion").beta0 == 0.04);
    CHECK(spec("self-motion").duration == 6.0);

    // The friction feedforward (CTRL-1): on in loops at the Gaz et al. 2019 levels (joint 6 raised), off everywhere else (the runs of
    // 2026-10-19 before it)
    for (const ExperimentSpec& e : all) {
        CAPTURE(e.name);
        CHECK(e.friction.active() == (e.name == "loops"));
        CHECK(e.friction.velocityScale == (e.name == "loops" ? 0.005 : 0.01));
    }
    Vec7 expected = FrictionFeedforward::gaz2019().coulomb;
    expected[0] = 0.70;  // joint 1: stuck under 0.58 N m (loops_7) and 0.80 N m (loops_6)
    expected[4] = 0.60;  // the wrist: below each joint's breakaway measured on 2026-10-19 (loops_2 to loops_5)
    expected[5] = 0.35;
    expected[6] = 0.50;
    CHECK(maxDiff(spec("loops").friction.coulomb, expected) == 0.0);
    CHECK(maxDiff(spec("loops").friction.viscous, (Vec7() << 0, 0, 0, 0, 0.5, 0.5, 0.5).finished()) == 0.0);
    CHECK(spec("loops").friction.coulomb.maxCoeff() < 0.75);  // joint 1: 0.70 N m, the largest
    CHECK(spec("loops").friction.velocityScale == 0.005);
    // The rotor inertia: 0.3 kg m^2 on the wrist of loops, none elsewhere
    for (const ExperimentSpec& e : all) {
        CAPTURE(e.name);
        CHECK(e.rotorInertia.isZero() == (e.name != "loops"));
    }
    CHECK(maxDiff(spec("loops").rotorInertia, (Vec7() << 0, 0, 0, 0, 0.5, 0.5, 0.5).finished()) == 0.0);
}

TEST_CASE("SpeedProfile: values, an exact time derivative, continuity at the joints, zero outside")
{
    const SpeedProfile p{0.08, 2.0, 14.0, 4.0};
    double nu = 0.0, dnu = 0.0;
    p.eval(-1.0, nu, dnu);
    CHECK(nu == 0.0);
    CHECK(dnu == 0.0);
    p.eval(0.0, nu, dnu);
    CHECK(nu == 0.0);
    CHECK(dnu == 0.0);
    p.eval(1.0, nu, dnu);
    CHECK_THAT(nu, WithinAbs(0.04, 1e-15));  // half way up: r(1/2) = 1/2
    CHECK_THAT(dnu, WithinAbs(0.08 * 1.875 / 2.0, 1e-15));
    p.eval(10.0, nu, dnu);
    CHECK(nu == 0.08);
    CHECK(dnu == 0.0);
    p.eval(18.0, nu, dnu);
    CHECK_THAT(nu, WithinAbs(0.04, 1e-15));
    CHECK_THAT(dnu, WithinAbs(-0.08 * 1.875 / 4.0, 1e-15));
    p.eval(20.0, nu, dnu);
    CHECK(nu == 0.0);
    CHECK(dnu == 0.0);
    p.eval(std::nan(""), nu, dnu);
    CHECK(nu == 0.0);
    CHECK(dnu == 0.0);
    // distance: 0.08 (1 + 14 + 2) = 1.36
    CHECK_THAT(p.distance(), WithinAbs(1.36, 1e-14));

    // d nu / dt against a central difference, everywhere including near the joints of the pieces
    const double h = 1e-6;
    for (double t = 0.05; t < 20.0; t += 0.1) {
        double np = 0.0, nm = 0.0, d0 = 0.0, unused = 0.0;
        p.eval(t + h, np, unused);
        p.eval(t - h, nm, unused);
        p.eval(t, nu, d0);
        CAPTURE(t);
        CHECK(std::abs((np - nm) / (2.0 * h) - d0) < 1e-8);
    }
    // nu and d nu / dt are continuous at the joints 2, 16 and 20 s
    for (double tj : {2.0, 16.0, 20.0}) {
        double a = 0.0, da = 0.0, b = 0.0, db = 0.0;
        p.eval(tj - 1e-9, a, da);
        p.eval(tj + 1e-9, b, db);
        CAPTURE(tj);
        CHECK(std::abs(a - b) < 1e-12);
        CHECK(std::abs(da - db) < 1e-9);
    }
    // The integral of nu is the distance (trapezoidal rule, 1 ms)
    double integral = 0.0, prev = 0.0;
    for (int k = 1; k <= 21000; ++k) {
        p.eval(1e-3 * k, nu, dnu);
        integral += 0.5e-3 * (prev + nu);
        prev = nu;
    }
    CHECK_THAT(integral, WithinAbs(p.distance(), 1e-8));
}

// ---------------------------------------------------------------------------------------------------------------
// The friction feedforward (CTRL-1)

TEST_CASE("FrictionFeedforward: the law, its limits, Gaz et al. 2019 and scaling")
{
    FrictionFeedforward f;
    CHECK_FALSE(f.active());
    CHECK(f.torque(Vec7::Constant(1.0)).isZero());  // off: nothing, whatever v
    f.coulomb = (Vec7() << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7).finished();
    f.viscous = (Vec7() << 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.5).finished();
    f.velocityScale = 0.01;
    CHECK(f.active());
    // At rest nothing; odd in v; the sign of the motion
    CHECK(f.torque(Vec7::Zero()).isZero());
    const Vec7 v = (Vec7() << 0.05, -0.05, 0.1, -0.2, 1.0, -1.0, 0.03).finished();
    CHECK(maxDiff(f.torque(v), Vec7(-f.torque(Vec7(-v)))) < 1e-15);
    for (int i = 0; i < 7; ++i) {
        CAPTURE(i);
        CHECK_THAT(f.torque(v)[i], WithinAbs(f.coulomb[i] * std::tanh(v[i] / 0.01) + f.viscous[i] * v[i], 1e-15));
        CHECK(f.torque(v)[i] * v[i] > 0.0);
    }
    // Saturated far above v0: Fc sign(v) + Fv v (tanh(50) = 1 to double precision)
    CHECK(maxDiff(f.torque(Vec7::Constant(0.5)), Vec7(f.coulomb + 0.5 * f.viscous)) < 1e-15);
    CHECK(maxDiff(f.torque(Vec7::Constant(-0.5)), Vec7(-f.coulomb - 0.5 * f.viscous)) < 1e-15);
    // At v = v0 76 % of Fc, at 3 v0 99.5 %: the transition is +-3 v0 wide
    CHECK_THAT(f.torque(Vec7::Constant(0.01))[0], WithinAbs(0.1 * std::tanh(1.0), 1e-15));
    CHECK(std::tanh(1.0) > 0.76);
    CHECK(std::tanh(3.0) > 0.995);
    // Bounded by Fc + Fv |v| for every finite v; a non-finite v gives a non-finite torque on that joint only
    // (ExperimentLoop turns it into NonFinite; 0 * inf is NaN, so an infinite v never saturates quietly)
    for (double x : {1e-6, 1e-3, 0.02, 0.3, 7.0, 1e6}) {
        const Vec7 vx = Vec7::Constant(x);
        CHECK((f.torque(vx).array() <= (f.coulomb + x * f.viscous).array() + 1e-15).all());
        CHECK(f.torque(vx).allFinite());
    }
    Vec7 bad = v;
    bad[2] = std::nan("");
    bad[0] = std::numeric_limits<double>::infinity();
    CHECK(std::isnan(f.torque(bad)[2]));
    CHECK_FALSE(std::isfinite(f.torque(bad)[0]));
    CHECK(std::isfinite(f.torque(bad)[1]));
    CHECK(std::isfinite(f.torque(bad)[6]));
    // Scaling: Fc and Fv, not v0; 0 turns it off
    const FrictionFeedforward h = f.scaled(0.5);
    CHECK(maxDiff(h.coulomb, Vec7(0.5 * f.coulomb)) == 0.0);
    CHECK(maxDiff(h.viscous, Vec7(0.5 * f.viscous)) == 0.0);
    CHECK(h.velocityScale == 0.01);
    CHECK_FALSE(f.scaled(0.0).active());
    CHECK(maxDiff(h.torque(v), Vec7(0.5 * f.torque(v))) < 1e-15);

    // Gaz et al. 2019: Fc = phi1 / 2, the mean of the two asymptotes of their sigmoid (recomputed here from the
    // published phi1, phi2, phi3), no viscous term, v0 = 0.01 rad/s
    const FrictionFeedforward g = FrictionFeedforward::gaz2019();
    CHECK(g.active());
    CHECK(g.viscous.isZero());
    CHECK(g.velocityScale == 0.01);
    const double phi1[7] = {0.54615, 0.87224, 0.64068, 1.2794, 0.83904, 0.30301, 0.56489};
    const double phi2[7] = {5.1181, 9.0657, 10.136, 5.5903, 8.3469, 17.133, 10.336};
    const double phi3[7] = {0.039533, 0.025882, -0.04607, 0.036194, 0.026226, -0.021047, 0.0035526};
    const double expected[7] = {0.2731, 0.4361, 0.3203, 0.6397, 0.4195, 0.1515, 0.2824};
    for (int i = 0; i < 7; ++i) {
        CAPTURE(i);
        const double sigma = 1.0 / (1.0 + std::exp(-phi2[i] * phi3[i]));
        const double offset = phi1[i] * sigma;  // the sigmoid's value subtracted at qd = 0
        const double plus = phi1[i] - offset;    // qd -> +inf
        const double minus = offset;             // |tau_f| as qd -> -inf
        CHECK_THAT(g.coulomb[i], WithinAbs(0.5 * (plus + minus), 1e-12));
        CHECK_THAT(g.coulomb[i], WithinAbs(expected[i], 1e-4));
        CHECK(std::abs(plus - minus) <= 0.5 * g.coulomb[i]);  // the asymmetry dropped: at most +-25 %
        // The sigmoid's slope at rest is its transition width, not a viscous coefficient: at 0.05 rad/s it gives
        // only 12 to 42 % of Fc (the stalled joints need the full level from the start)
        const double at005 = phi1[i] / (1.0 + std::exp(-phi2[i] * (0.05 + phi3[i]))) - offset;
        CHECK(at005 < 0.45 * g.coulomb[i]);
    }
}

TEST_CASE("ExperimentSetup refuses an invalid friction feedforward")
{
    const PandaExample& ex = pandaExample();
    ExperimentSpec bad = spec("loops");
    bad.friction.coulomb[4] = -0.1;
    CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
    bad = spec("loops");
    bad.friction.coulomb[0] = std::nan("");
    CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
    bad = spec("loops");
    bad.friction.viscous[6] = -1.0;
    CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
    bad = spec("loops");
    bad.friction.viscous[1] = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
    for (double v0 : {0.0, -0.01, std::nan(""), std::numeric_limits<double>::infinity()}) {
        CAPTURE(v0);
        bad = spec("loops");
        bad.friction.velocityScale = v0;
        CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
    }
    // Off (every other experiment) and a wide transition are valid
    ExperimentSpec ok = spec("paper");
    ok.friction.velocityScale = 10.0;
    CHECK_NOTHROW(ExperimentSetup(ex, ok));
}

TEST_CASE("ExperimentLoop: tau_f follows the reference velocity along the lift and counts in the torque checks")
{
    const PandaExample& ex = pandaExample();
    const ExperimentSetup s(ex, spec("loops"));
    const FrictionFeedforward& f = s.spec().friction;
    REQUIRE(f.active());
    ModelTerms m;
    m.M = 2.0 * Mat<7, 7>::Identity();
    m.c = Vec7::Constant(0.1);
    m.G = Vec7::Constant(1.0);  // small: the artificial state below (at rest, 0.05 behind nu_d) commands a few N m

    // referenceSpeed: the profile's nu_d for a velocity reference, 0 for a position reference
    double nu = 0.0, dnu = 0.0;
    spec("loops").speed.eval(10.0, nu, dnu);
    CHECK(nu == 0.05);
    CHECK(s.referenceSpeed(10.0) == nu);
    CHECK(s.referenceSpeed(0.0) == 0.0);
    CHECK(s.referenceSpeed(-1.0) == 0.0);
    CHECK(s.referenceSpeed(100.0) == 0.0);
    spec("loops").speed.eval(2.0, nu, dnu);
    CHECK(s.referenceSpeed(2.0) == nu);
    CHECK(ExperimentSetup(ex, spec("hold")).referenceSpeed(3.0) == 0.0);
    CHECK(ExperimentSetup(ex, spec("self-motion")).referenceSpeed(3.0) == 0.0);

    // At t = 0 the profile is at rest: no feedforward, and the first torque is c alone (the mu ramp zeroes qdd*)
    {
        ExperimentLoop loop(s, s.monitor());
        ExperimentLoop::Cycle c;
        loop.follow(0.0, s.start(), Vec7::Zero(), &m, c);
        CHECK(c.violation == Violation::None);
        CHECK(c.tauF.isZero());
        CHECK(maxDiff(c.tau, m.c) < 1e-12);
    }
    // Cruising (nu_d = 0.05 at t = 10 s) from rest on the lift at s = 0.12: v = q*'(0.12) nu_d, tau_f = Fc tanh(v
    // / v0) with the sign of q*', and tau = M qdd* + c + tau_f
    Vec7 qs = Vec7::Zero(), dqs = Vec7::Zero();
    ex.lift().eval(0.12, qs, &dqs);
    const Vec7 v = dqs * 0.05;
    {
        ExperimentLoop loop(s, s.monitor());
        ExperimentLoop::Cycle c;
        // A fresh loop compares its first torque with zero (the arm rests under gravity compensation before the
        // loop); the torque of this artificial state (at rest, 0.05 behind nu_d) exceeds 0.6 N m, so warm up with
        // one cycle and check the second
        loop.follow(10.0, qs, Vec7::Zero(), &m, c);
        CHECK((c.violation == Violation::None || c.violation == Violation::TorqueRate));
        loop.follow(10.0, qs, Vec7::Zero(), &m, c);
        REQUIRE(c.violation == Violation::None);
        CHECK_THAT(c.status.eta, WithinAbs(0.12, 1e-9));
        CHECK(maxDiff(c.tauF, f.torque(v)) < 1e-12);
        CHECK(c.tauF.cwiseAbs().maxCoeff() > 0.05);  // felt
        for (int i = 0; i < 7; ++i) {
            CAPTURE(i);
            CHECK(c.tauF[i] * dqs[i] >= 0.0);
            CHECK(std::abs(c.tauF[i]) <= f.coulomb[i] + f.viscous[i] * std::abs(v[i]) + 1e-15);
        }
        CHECK(maxDiff(c.tau, Vec7(2.0 * c.qdd + s.spec().rotorInertia.cwiseProduct(c.qdd) + m.c + c.tauF)) < 1e-12);
        CHECK(maxDiff(c.tauG, Vec7(c.tau + m.G)) < 1e-12);
        // Kinematic only (no model): no torque, but tau_f is still reported (the preflight's log)
        ExperimentLoop kinematic(s, s.monitor());
        kinematic.follow(10.0, qs, Vec7::Zero(), nullptr, c);
        CHECK(c.violation == Violation::None);
        CHECK(c.tau.isZero());
        CHECK(c.tauG.isZero());
        CHECK(maxDiff(c.tauF, f.torque(v)) < 1e-12);
    }
    // The rotor inertia: tau = (M + diag(Jr)) qdd* + c + tau_f; the preflight's motion is unchanged
    {
        ExperimentSpec e = spec("loops");
        e.rotorInertia = (Vec7() << 0.05, 0.0, 0.0, 0.0, 0.05, 0.1, 0.1).finished();
        const ExperimentSetup es(ex, e);
        ExperimentLoop loop(es, es.monitor());
        ExperimentLoop::Cycle c;
        loop.follow(10.0, qs, Vec7::Zero(), &m, c);
        loop.follow(10.0, qs, Vec7::Zero(), &m, c);
        REQUIRE(c.violation == Violation::None);
        CHECK(c.qdd.cwiseAbs().maxCoeff() > 1e-3);
        CHECK(maxDiff(c.tau, Vec7(2.0 * c.qdd + e.rotorInertia.cwiseProduct(c.qdd) + m.c + c.tauF)) < 1e-12);
        ExperimentSpec plain = spec("loops");
        plain.rotorInertia.setZero();
        const ExperimentSetup ps(ex, plain);
        ExperimentLoop loop2(ps, ps.monitor());
        ExperimentLoop::Cycle c2;
        loop2.follow(10.0, qs, Vec7::Zero(), &m, c2);
        loop2.follow(10.0, qs, Vec7::Zero(), &m, c2);
        CHECK(maxDiff(c2.qdd, c.qdd) == 0.0);  // the same qdd*: only the torque differs
        CHECK(maxDiff(c2.tau, Vec7(c.tau - e.rotorInertia.cwiseProduct(c.qdd))) < 1e-12);
        // Invalid values are refused
        e.rotorInertia[4] = -0.1;
        CHECK_THROWS_AS(ExperimentSetup(ex, e), std::invalid_argument);
        e.rotorInertia[4] = std::nan("");
        CHECK_THROWS_AS(ExperimentSetup(ex, e), std::invalid_argument);
    }
    // Off (the other experiments): tau_f = 0 and tau = M qdd* + c, as before CTRL-1
    {
        const ExperimentSetup p(ex, spec("paper"));
        ExperimentLoop loop(p, p.monitor());
        ExperimentLoop::Cycle c;
        loop.follow(10.0, qs, Vec7::Zero(), &m, c);  // warm-up (the rate check against zero), as above
        loop.follow(10.0, qs, Vec7::Zero(), &m, c);
        REQUIRE(c.violation == Violation::None);
        CHECK(c.tauF.isZero());
        CHECK(maxDiff(c.tau, Vec7(2.0 * c.qdd + m.c)) < 1e-12);
    }
    // A position reference holds a point: v = 0 and tau_f = 0 however large the parameters
    {
        ExperimentSpec h = spec("hold");
        h.friction = FrictionFeedforward::gaz2019();
        const ExperimentSetup hs(ex, h);
        ExperimentLoop loop(hs, hs.monitor());
        ExperimentLoop::Cycle c;
        loop.follow(3.0, hs.start(), Vec7::Constant(0.01), &m, c);
        CHECK(c.violation == Violation::None);
        CHECK(c.tauF.isZero());
        CHECK(maxDiff(c.tau, Vec7(2.0 * c.qdd + m.c)) < 1e-12);
    }
    // tau_f counts in the torque check: 8 N m of Fc on joint 6 exceeds 60 % of its 12 N m once the lift moves it
    ModelTerms light;
    light.M = 2.0 * Mat<7, 7>::Identity();
    {
        ExperimentSpec e = spec("loops");
        e.friction.coulomb[5] = 8.0;
        e.friction.velocityScale = 1e-6;  // saturated at any speed
        REQUIRE(dqs[5] != 0.0);
        const ExperimentSetup es(ex, e);
        ExperimentLoop loop(es, es.monitor());
        ExperimentLoop::Cycle c;
        loop.follow(10.0, qs, Vec7::Zero(), &light, c);
        CHECK(c.violation == Violation::Torque);
        CHECK_THAT(std::abs(c.tauF[5]), WithinAbs(8.0 + e.friction.viscous[5] * std::abs(v[5]), 1e-9));
        CHECK(maxDiff(c.tau, Vec7(2.0 * c.qdd + e.rotorInertia.cwiseProduct(c.qdd) + c.tauF)) < 1e-12);
        CHECK(std::abs(c.tau[5]) > 0.6 * 12.0);  // beyond the 60 % margin of joint 6
    }
    // ... and in the torque-rate check: a profile that starts at nu_d at t = 0 with Fc = 1 N m makes the first
    // command jump from 0 to 1 N m: 1000 N m/s. The catalog's ramped profile starts at nu_d = 0, so its tau_f does.
    {
        ExperimentSpec e = spec("loops");
        e.speed = {0.05, 0.0, 34.8, 0.0};
        e.friction.coulomb = Vec7::Constant(1.0);
        e.friction.velocityScale = 1e-6;
        const ExperimentSetup es(ex, e);
        ExperimentLoop loop(es, es.monitor());
        ExperimentLoop::Cycle c;
        loop.follow(0.0, es.start(), Vec7::Zero(), &light, c);
        CHECK(c.violation == Violation::TorqueRate);
        CHECK(c.tauF.cwiseAbs().maxCoeff() > 0.999);
        ExperimentSpec ramped = e;
        ramped.speed = spec("loops").speed;
        const ExperimentSetup rs(ex, ramped);
        ExperimentLoop loop2(rs, rs.monitor());
        loop2.follow(0.0, rs.start(), Vec7::Zero(), &light, c);
        CHECK(c.violation == Violation::None);
        CHECK(c.tauF.isZero());
    }
}

TEST_CASE("preflight with the friction feedforward: loops passes with a model, the log and the peaks report tau_f, "
          "the motion is unchanged",
          SLOW_IN_DEBUG)
{
    const PandaExample& ex = pandaExample();
    const FakeDynamics dyn;
    const ExperimentSetup s(ex, spec("loops"));
    const FrictionFeedforward& f = s.spec().friction;
    std::vector<ExperimentLogRow> log;
    const PreflightReport r = preflight(s, &dyn, s.monitor(), &log);
    INFO(r.summary());
    CHECK(r.ok);
    checkWithinLimits(r);
    // The peak: at least the largest Fc (joint 4, 0.64 N m), at most Fc + Fv |v| with the lift's speeds below
    // 0.25 rad/s at nu_d (joint 5: 0.60 + 0.5 |v|); the summary names it
    CHECK(r.peaks.maxTauF >= f.coulomb.maxCoeff() - 1e-6);  // tanh saturates
    CHECK(r.peaks.maxTauF <= (f.coulomb + 0.25 * f.viscous).maxCoeff());
    char peak[64];
    std::snprintf(peak, sizeof(peak), "tau_f <= %.2f N m (joint %d)", r.peaks.maxTauF, r.peaks.maxTauFJoint + 1);
    CHECK(r.summary().find(peak) != std::string::npos);
    // Smooth: the torque rate stays far below the 600 N m/s margin (the sign changes go through tanh)
    CAPTURE(r.peaks.maxTauRate);
    CHECK(r.peaks.maxTauRate < 20.0);
    // Every row: tau_f is the law at q*'(eta) nu_d(t), zero while the profile is at rest, and tau includes it
    REQUIRE(log.size() > 35000);
    CHECK(log.front().tauF.isZero());
    int felt = 0;
    for (std::size_t k = 0; k < log.size(); k += 97) {
        const ExperimentLogRow& row = log[k];
        Vec7 qs = Vec7::Zero(), dqs = Vec7::Zero();
        ex.lift().eval(row.eta, qs, &dqs);
        const Vec7 expected = f.torque(Vec7(dqs * s.referenceSpeed(row.t)));
        CHECK(maxDiff(row.tauF, expected) < 1e-12);
        CHECK(maxDiff(row.tau, Vec7(dyn.mass(row.q) * row.qdd + s.spec().rotorInertia.cwiseProduct(row.qdd)
                                    + dyn.coriolis(row.q, row.dq) + row.tauF))
              < 1e-12);
        if (row.t > 4.0 && row.t < 30.8 && row.tauF.cwiseAbs().maxCoeff() > 0.1)
            ++felt;
        if (row.t >= 34.8)
            CHECK(row.tauF.isZero());
    }
    CHECK(felt > 200);
    // The feedforward does not change the predicted motion (the preflight integrates qdd*): without it the same
    // q, qdd and the same torque less tau_f
    ExperimentSpec offSpec = spec("loops");
    offSpec.friction = FrictionFeedforward();
    const ExperimentSetup off(ex, offSpec);
    std::vector<ExperimentLogRow> logOff;
    const PreflightReport ro = preflight(off, &dyn, off.monitor(), &logOff);
    REQUIRE(ro.ok);
    REQUIRE(logOff.size() == log.size());
    CHECK(ro.peaks.maxTauF == 0.0);
    CHECK(ro.summary().find("tau_f <= 0.00 N m") != std::string::npos);
    for (std::size_t k = 0; k < log.size(); k += 211) {
        CHECK(maxDiff(logOff[k].q, log[k].q) == 0.0);
        CHECK(maxDiff(logOff[k].qdd, log[k].qdd) == 0.0);
        CHECK(maxDiff(logOff[k].tau, Vec7(log[k].tau - log[k].tauF)) < 1e-12);
        CHECK(logOff[k].tauF.isZero());
    }
    CHECK(ro.peaks.maxDdqFrac == r.peaks.maxDdqFrac);
    CHECK(ro.peaks.maxXi == r.peaks.maxXi);
}

// ---------------------------------------------------------------------------------------------------------------
// The setup and the control cycle

TEST_CASE("ExperimentSetup: the start has the experiment's outputs; the controller is configured as specified")
{
    const PandaExample& ex = pandaExample();
    for (const ExperimentSpec& e : experiments()) {
        CAPTURE(e.name);
        const ExperimentSetup s(ex, e);
        CHECK(s.s0() == 0.1);
        CHECK(s.duration() == e.duration);
        CHECK(s.timeout() == e.duration + 5.0);
        const Controller<PandaArm> c = s.controller();
        CHECK(c.rampTime() == e.rampTime);
        CHECK_THAT(c.prediction(), WithinAbs(0.1, 1e-15));
        CHECK(c.gains().Kp == e.gains.Kp);
        CHECK(c.gains().kdEta == e.gains.kdEta);
        CHECK(c.reference().mode == e.mode);
        CHECK(c.redundancy() == &ex.redundancy());
        Controller<PandaArm>::OutputState out;
        REQUIRE(c.outputs(s.start(), 0.1, out));
        CHECK_THAT(out.eta, WithinAbs(0.1, 1e-10));
        CHECK(maxDiff(out.xi, e.xi0) < 1e-10);
        CHECK_THAT(out.beta, WithinAbs(e.beta0, 1e-10));
        // The reference: eta_d = s0 (position), or the speed profile
        double ref = 0.0, d1 = 0.0, d2 = 0.0;
        const TangentialReference r = s.reference();
        if (e.mode == TangentialReference::Mode::Position) {
            r.position(3.0, ref, d1, d2);
            CHECK(ref == 0.1);
            CHECK(d1 == 0.0);
            CHECK(d2 == 0.0);
        } else {
            double nu = 0.0, dnu = 0.0;
            e.speed.eval(5.0, nu, dnu);
            r.velocity(5.0, 0.3, ref, d1, d2);
            CHECK(ref == nu);
            CHECK(d1 == dnu);
            CHECK(d2 == 0.0);
        }
        CHECK_FALSE(s.monitor().clamped());
    }
    // The paper start is PandaExample's
    CHECK(maxDiff(ExperimentSetup(ex, spec("paper")).start(), ex.paperStart()) < 1e-12);

    // Invalid specs are refused
    ExperimentSpec bad = spec("loops");
    bad.duration = 30.0;  // the speed profile lasts 34.8 s
    CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
    bad = spec("hold");
    bad.duration = 0.0;
    CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
    bad.duration = std::nan("");
    CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
    bad = spec("hold");
    bad.rampTime = -1.0;
    CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
    bad = spec("loops");
    bad.speed.rampUp = -1.0;
    CHECK_THROWS_AS(ExperimentSetup(ex, bad), std::invalid_argument);
}

TEST_CASE("ExperimentLoop: the cycle's order of checks, the torque, the end condition, the timeout, the damping stop")
{
    const PandaExample& ex = pandaExample();
    const ExperimentSetup s(ex, spec("hold"));
    const Vec7 q0 = s.start();

    // At rest at the start: no command, no violation; not finished before T, finished after it
    {
        ExperimentLoop loop(s, s.monitor());
        ExperimentLoop::Cycle c;
        loop.follow(0.0, q0, Vec7::Zero(), nullptr, c);
        CHECK(c.violation == Violation::None);
        CHECK_FALSE(c.finished);
        CHECK(c.qdd.cwiseAbs().maxCoeff() < 1e-12);  // mu ramp: the first command is zero
        CHECK(c.tau.cwiseAbs().maxCoeff() == 0.0);  // no model: no torque
        loop.follow(s.duration() - 1e-3, q0, Vec7::Zero(), nullptr, c);
        CHECK_FALSE(c.finished);
        loop.follow(s.duration(), q0, Vec7::Zero(), nullptr, c);
        CHECK(c.finished);
        // Moving: not finished
        loop.follow(s.duration() + 0.1, q0, Vec7::Constant(0.01), nullptr, c);
        CHECK_FALSE(c.finished);
        CHECK(c.violation == Violation::None);
        // Timeout after T + 5 s
        loop.follow(s.timeout() + 1e-3, q0, Vec7::Constant(0.01), nullptr, c);
        CHECK(c.violation == Violation::Timeout);
        CHECK_FALSE(c.finished);
    }
    // With a model: tau = M qdd* + c, tauG = tau + G; the rate is checked against the previous cycle's tau
    {
        ExperimentLoop loop(s, s.monitor());
        ModelTerms m;
        m.M = 2.0 * Mat<7, 7>::Identity();
        m.c = Vec7::Constant(0.1);
        m.G = Vec7::Constant(3.0);
        ExperimentLoop::Cycle c;
        loop.follow(0.0, q0, Vec7::Zero(), &m, c);
        CHECK(c.violation == Violation::None);
        CHECK(maxDiff(c.tau, Vec7(2.0 * c.qdd + m.c)) < 1e-15);
        CHECK(maxDiff(c.tauG, Vec7(c.tau + m.G)) < 1e-15);
        // A jump of c by 1 N m in one period is 1000 N m/s
        m.c = Vec7::Constant(1.1);
        loop.follow(1e-3, q0, Vec7::Zero(), &m, c);
        CHECK(c.violation == Violation::TorqueRate);
        // A G that the joint cannot deliver within 60 %
        ExperimentLoop loop2(s, s.monitor());
        m.c.setZero();
        m.G = Vec7::Zero();
        m.G[5] = 8.0;  // joint 6: 8 > 7.2 N m
        loop2.follow(0.0, q0, Vec7::Zero(), &m, c);
        CHECK(c.violation == Violation::Torque);
    }
    // Order: the state checks come after the outputs, the acceleration after the state
    {
        ExperimentLoop loop(s, s.monitor());
        ExperimentLoop::Cycle c;
        Vec7 dq = Vec7::Zero();
        dq[3] = 0.5 * 2.175 * 1.01;  // above 50 % of joint 4's limit, and above the hold's tighter 20 %
        loop.follow(0.5, q0, dq, nullptr, c);
        CHECK(c.violation == Violation::JointVelocity);
        Vec7 qn = q0;
        qn[2] = std::nan("");
        ExperimentLoop loop2(s, s.monitor());
        loop2.follow(0.5, qn, Vec7::Zero(), nullptr, c);
        CHECK(c.violation == Violation::NonFinite);
        // Far from the path the projection fails, or |xi| is too large
        ExperimentLoop loop3(s, s.monitor());
        Vec7 qf = q0;
        qf[0] += 0.5;
        loop3.follow(0.5, qf, Vec7::Zero(), nullptr, c);
        CHECK((c.violation == Violation::Projection || c.violation == Violation::Tube));
    }
    // observe: the outputs without a command, for the log of the stop
    {
        ExperimentLoop loop(s, s.monitor());
        ExperimentLoop::Status st;
        REQUIRE(loop.observe(q0, Vec7::Zero(), st));
        CHECK_THAT(st.eta, WithinAbs(0.1, 1e-10));
        CHECK(st.xi.norm() < 1e-10);
        CHECK(st.etaDot == 0.0);
        Vec7 qn = q0;
        qn[1] = std::nan("");
        CHECK_FALSE(loop.observe(qn, Vec7::Zero(), st));
    }
}

TEST_CASE("SafetyStop: the damped hold at q_stop, its clamp, non-finite values and its end")
{
    const Vec7 K = RunConstants::stopStiffness();
    const Vec7 D = RunConstants::stopDamping();
    CHECK(maxDiff(D, (Vec7() << 30, 30, 30, 30, 3, 3, 1.5).finished()) == 0.0);
    CHECK(maxDiff(K, (Vec7() << 300, 300, 300, 300, 50, 50, 20).finished()) == 0.0);
    // At most the gains of libfranka's joint impedance example (examples/joint_impedance_control.cpp), on the Panda
    CHECK((K.array() <= (Vec7() << 600, 600, 600, 600, 250, 150, 50).finished().array()).all());
    CHECK((D.array() <= (Vec7() << 50, 50, 50, 50, 30, 25, 15).finished().array()).all());
    const Vec7 limit = 0.6 * PandaArm::tauMax();
    CHECK(RunConstants::stopTorqueFrac == HardwareBounds::tauFrac);
    const Vec7 qs = PandaArm::restPosture();
    const Vec7 dq = (Vec7() << 0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7).finished();
    const Vec7 dqs = (Vec7() << 0.01, -0.02, 0.03, -0.04, 0.05, -0.06, 0.07).finished();

    SECTION("the law: tau = -K_s (q - q_stop) - D qd")
    {
        CHECK(maxDiff(SafetyStop::torque(qs, dq, qs), Vec7(-D.cwiseProduct(dq))) < 1e-15);
        const Vec7 q = qs + dqs;
        CHECK(maxDiff(SafetyStop::torque(q, Vec7::Zero(), qs), Vec7(-K.cwiseProduct(dqs))) < 1e-12);
        CHECK(maxDiff(SafetyStop::torque(q, dq, qs), Vec7(-K.cwiseProduct(dqs) - D.cwiseProduct(dq))) < 1e-12);
        // The hold: a 0.73 kg weight at 0.5 m (3.6 N m on joint 2) deflects it by less than 0.02 rad
        CHECK(0.73 * 9.81 * 0.5 / K[1] < 0.02);
        // At 50 % of every velocity limit and at q_stop the damping alone stays inside the clamp
        const Vec7 dqHalf = 0.5 * PandaArm::dqMax();
        CHECK((SafetyStop::torque(qs, dqHalf, qs).cwiseAbs().array() < limit.array()).all());
    }
    SECTION("each joint clamped to 60 % of its torque limit")
    {
        for (int i = 0; i < 7; ++i) {
            CAPTURE(i);
            Vec7 q = qs;
            q[i] += 1.0;  // far beyond any deflection: K_s alone asks for 300 N m on joints 1-4
            const Vec7 tau = SafetyStop::torque(q, Vec7::Zero(), qs);
            CHECK(tau[i] == -limit[i]);
            Vec7 v = Vec7::Zero();
            v[i] = -1e3;
            CHECK(SafetyStop::torque(qs, v, qs)[i] == limit[i]);
        }
    }
    SECTION("non-finite values: damping only without q, stiffness only without qd, never a non-finite torque")
    {
        Vec7 q = qs + dqs;
        q[2] = std::nan("");
        Vec7 v = dq;
        v[4] = std::numeric_limits<double>::infinity();
        v[2] = std::nan("");
        Vec7 qStop = qs;
        qStop[0] = std::nan("");
        const Vec7 tau = SafetyStop::torque(q, v, qStop);
        CHECK(tau.allFinite());
        CHECK(tau[0] == -D[0] * v[0]);                              // q_stop not finite: damping only
        CHECK(tau[2] == 0.0);                                       // neither q nor qd finite
        CHECK_THAT(tau[4], WithinAbs(-K[4] * dqs[4], 1e-12));       // qd not finite: stiffness only
        CHECK_THAT(tau[1], WithinAbs(-K[1] * dqs[1] - D[1] * v[1], 1e-12));
    }
    SECTION("the stop latches q_stop and ends after stopSettle at rest, or after stopTime")
    {
        SafetyStop stop;
        bool done = false;
        CHECK_FALSE(stop.active());
        CHECK(stop.step(0.0, qs, dq, done).isZero());  // not active: nothing, done
        CHECK(done);

        stop.begin(1.0, qs);
        CHECK(stop.active());
        CHECK(stop.startTime() == 1.0);
        CHECK(maxDiff(stop.holdConfiguration(), qs) == 0.0);
        // Moving: not done
        CHECK(maxDiff(stop.step(1.0, qs, dq, done), Vec7(-D.cwiseProduct(dq))) < 1e-15);
        CHECK_FALSE(done);
        // At rest from t = 1.2 s: done only after 100 ms of rest
        int k = 200;
        for (; k < 400; ++k) {
            stop.step(1.0 + 1e-3 * k, qs, Vec7::Constant(1e-3), done);
            if (done)
                break;
        }
        CHECK((k == 300 || k == 301));  // 100 ms after the rest began (floating-point time)
        // A motion in between restarts the rest timer
        SafetyStop s2;
        s2.begin(0.0, qs);
        s2.step(0.00, qs, Vec7::Zero(), done);
        s2.step(0.05, qs, Vec7::Zero(), done);
        CHECK_FALSE(done);
        s2.step(0.06, qs, Vec7::Constant(0.01), done);  // |qd| = 0.026 rad/s: moving again
        s2.step(0.12, qs, Vec7::Zero(), done);
        CHECK_FALSE(done);  // 0.06 s of rest after the motion, not 0.12 s
        s2.step(0.23, qs, Vec7::Zero(), done);
        CHECK(done);
        // Never at rest (or a NaN velocity): done at stopTime
        SafetyStop s3;
        s3.begin(2.0, qs);
        s3.step(2.0 + RunConstants::stopTime - 1e-3, qs, dq, done);
        CHECK_FALSE(done);
        s3.step(2.0 + RunConstants::stopTime, qs, Vec7::Constant(std::nan("")), done);
        CHECK(done);
        // A joint not finite when the stop began is held where it is measured again
        Vec7 qn = qs;
        qn[3] = std::nan("");
        SafetyStop s5;
        s5.begin(0.0, qn);
        CHECK(s5.step(0.0, qn, Vec7::Zero(), done).allFinite());
        CHECK(std::isnan(s5.holdConfiguration()[3]));
        Vec7 q1 = qs;
        q1[3] += 0.01;
        s5.step(1e-3, q1, Vec7::Zero(), done);
        CHECK(s5.holdConfiguration()[3] == q1[3]);
        Vec7 q2 = q1;
        q2[3] += 0.01;
        CHECK_THAT(s5.step(2e-3, q2, Vec7::Zero(), done)[3], WithinAbs(-K[3] * 0.01, 1e-9));
        // A NaN time ends it
        SafetyStop s4;
        s4.begin(0.0, qs);
        s4.step(std::nan(""), qs, dq, done);
        CHECK(done);
    }
}

// ---------------------------------------------------------------------------------------------------------------
// The preflight

TEST_CASE("preflight: every experiment passes, kinematic only; the peaks of the design doc", SLOW_IN_DEBUG)
{
    const PandaExample& ex = pandaExample();
    for (const ExperimentSpec& e : experiments()) {
        CAPTURE(e.name);
        const ExperimentSetup s(ex, e);
        std::vector<ExperimentLogRow> log;
        const PreflightReport r = preflight(s, nullptr, s.monitor(), &log);
        INFO(r.summary());
        CHECK(r.ok);
        CHECK(r.reason.empty());
        CHECK(r.violation == Violation::None);
        CHECK_FALSE(r.torqueChecked);
        CHECK(r.peaks.maxTauFrac == 0.0);
        CHECK(r.converged);
        CHECK(r.summary().find("NOT checked") != std::string::npos);
        checkWithinLimits(r);
        // Ends at rest just after T (every reference ends smoothly)
        CHECK(r.simulated >= e.duration);
        CHECK(r.simulated <= e.duration + 0.5);
        CHECK(log.back().dq.norm() < RunConstants::restSpeed);
        // The start outputs
        CHECK_THAT(r.startXi, WithinAbs(e.xi0.norm(), 1e-10));
        CHECK_THAT(r.startBeta, WithinAbs(e.beta0, 1e-10));
        // The log has one row per simulated period, and its peaks are the report's
        CHECK(static_cast<long>(log.size()) == r.peaks.rows);
        CHECK_THAT(log.back().t, WithinAbs(r.simulated, 1e-12));
        RunPeaks again;
        for (const ExperimentLogRow& row : log)
            again.add(row);
        CHECK(again.maxXi == r.peaks.maxXi);
        CHECK(again.maxDdqFrac == r.peaks.maxDdqFrac);
        CHECK(again.minJointMargin == r.peaks.minJointMargin);
        // The arm follows the lift: beta stays below 1e-3 unless it starts off it
        if (e.beta0 == 0.0)
            CHECK(r.peaks.maxBeta < 1e-3);
    }
    // The paper's start with the ramp (design doc, Hardware: 11.5 % with nu_d = 0.08 from t = 0; the experiment's
    // nu_d ramp over 2 s lowers it), joint margin 0.359 rad at joint 4 as in the paper's run
    const PreflightReport p = preflight(ExperimentSetup(ex, spec("paper")), nullptr, SafetyMonitor());
    CAPTURE(p.peaks.maxDdqFrac, p.peaks.maxDqFrac, p.peaks.minJointMargin);
    CHECK(p.peaks.maxDdqFrac < 0.115);
    CHECK(p.peaks.maxDdqFrac > 0.05);
    CHECK_THAT(p.peaks.maxXi, WithinAbs(0.0312, 1e-4));
    CHECK_THAT(p.peaks.minJointMargin, WithinAbs(0.359, 1e-3));
    CHECK(p.peaks.minMarginJoint == 3);
    // The loops end where they start (two loops), on the lift
    std::vector<ExperimentLogRow> loopsLog;
    const ExperimentSetup loops(ex, spec("loops"));
    REQUIRE(preflight(loops, nullptr, loops.monitor(), &loopsLog).ok);
    CHECK(maxDiff(loopsLog.back().q, loops.start()) < 2e-3);
}

TEST_CASE("preflight: every experiment passes with a dynamics model; the torques are M qdd* + C qd and G",
          SLOW_IN_DEBUG)
{
    const PandaExample& ex = pandaExample();
    const FakeDynamics dyn;
    for (const ExperimentSpec& e : experiments()) {
        CAPTURE(e.name);
        const ExperimentSetup s(ex, e);
        std::vector<ExperimentLogRow> log;
        const PreflightReport r = preflight(s, &dyn, s.monitor(), &log);
        INFO(r.summary());
        CHECK(r.ok);
        CHECK(r.torqueChecked);
        CHECK(r.summary().find("NOT checked") == std::string::npos);
        checkWithinLimits(r);
        CHECK(r.peaks.maxTauFrac > 0.0);
        // The torques of the log, recomputed from the model at a few rows (plus the friction feedforward and the
        // rotor inertia: loops)
        for (std::size_t k = 0; k < log.size(); k += log.size() / 7) {
            const ExperimentLogRow& row = log[k];
            const Vec7 tau = dyn.mass(row.q) * row.qdd + e.rotorInertia.cwiseProduct(row.qdd)
                             + dyn.coriolis(row.q, row.dq) + row.tauF;
            CHECK(maxDiff(row.tau, tau) < 1e-12);
            CHECK(maxDiff(row.tauG, Vec7(tau + dyn.gravity(row.q))) < 1e-12);
        }
        // The mu ramp: the first command is zero, so the first torque rate is zero
        CHECK(log.front().tau.cwiseAbs().maxCoeff() < 1e-12);
    }
}

TEST_CASE("preflight: a bad experiment fails with the right reason, at the right time")
{
    const PandaExample& ex = pandaExample();

    SECTION("a start offset outside the tube margin: Tube at t = 0")
    {
        ExperimentSpec e = spec("paper");
        e.xi0 = 1.5 * PandaExample::xi0();  // |xi0| = 0.047 > 0.04
        const ExperimentSetup s(ex, e);
        const PreflightReport r = preflight(s, nullptr, s.monitor());
        CHECK_FALSE(r.ok);
        CHECK(r.violation == Violation::Tube);
        CHECK(r.tFail == 0.0);
        CHECK(r.reason == toString(Violation::Tube));
        CHECK(r.summary().find("FAILED") == 0);
    }
    SECTION("a start beta outside the redundancy margin: Redundancy at t = 0")
    {
        ExperimentSpec e = spec("self-motion");
        e.beta0 = 0.06;
        const ExperimentSetup s(ex, e);
        const PreflightReport r = preflight(s, nullptr, s.monitor());
        CHECK_FALSE(r.ok);
        CHECK(r.violation == Violation::Redundancy);
        CHECK(r.tFail == 0.0);
    }
    SECTION("huge gains: the commanded acceleration")
    {
        ExperimentSpec e = spec("paper");
        e.gains = Gains::fromLambda(30.0);
        const ExperimentSetup s(ex, e);
        const PreflightReport r = preflight(s, nullptr, s.monitor());
        CHECK_FALSE(r.ok);
        CHECK(r.violation == Violation::Acceleration);
        CHECK(r.tFail > 0.0);
        CHECK(r.tFail < 1.0);  // during the ramp
    }
    SECTION("the paper's start without the mu ramp: 40 % (49 % at constant nu_d) of an acceleration limit at t = 0")
    {
        ExperimentSpec e = spec("paper");
        e.rampTime = 0.0;
        const ExperimentSetup s(ex, e);
        std::vector<ExperimentLogRow> log;
        const PreflightReport r = preflight(s, nullptr, s.monitor(), &log);
        CHECK_FALSE(r.ok);
        CHECK(r.violation == Violation::Acceleration);
        CHECK(r.tFail == 0.0);
        REQUIRE(log.size() == 1);
        // 40.4 %: the experiment's nu_d starts at 0; with the paper's constant nu_d = 0.08 from t = 0 the
        // tangential law adds k_t nu_d = 0.4 and the first command is the design doc's 49.4 %
        CHECK_THAT(log[0].qdd.cwiseAbs().cwiseQuotient(PandaArm::ddqMax()).maxCoeff(), WithinAbs(0.404, 1e-3));
        e.speed = {PandaExample::nuD, 0.0, 20.0, 0.0};
        const PreflightReport c = preflight(ExperimentSetup(ex, e), nullptr, SafetyMonitor(), &log);
        CHECK(c.violation == Violation::Acceleration);
        REQUIRE(log.size() == 1);
        CHECK_THAT(log[0].qdd.cwiseAbs().cwiseQuotient(PandaArm::ddqMax()).maxCoeff(), WithinAbs(0.494, 1e-3));
    }
    SECTION("no mu ramp at half the offset: the torque jump fails the torque-rate margin")
    {
        ExperimentSpec e = spec("recover");
        e.rampTime = 0.0;  // the first command, 25 % of an acceleration limit, passes; its torque jumps from 0
        const ExperimentSetup s(ex, e);
        const ScaledInertia dyn(0.5);
        const PreflightReport r = preflight(s, &dyn, s.monitor());
        CHECK_FALSE(r.ok);
        CHECK(r.violation == Violation::TorqueRate);
        CHECK(r.tFail == 0.0);
        // Kinematic only it passes the first period: the rate is a torque check
        const PreflightReport k = preflight(s, nullptr, s.monitor());
        CHECK((k.ok || k.tFail > 0.0));
    }
    SECTION("a heavy model: the delivered torque")
    {
        const ExperimentSetup s(ex, spec("paper"));
        const ScaledInertia dyn(60.0);  // 60 kg m^2 on every joint: the wrist torques exceed 7.2 N m
        const PreflightReport r = preflight(s, &dyn, s.monitor());
        CHECK_FALSE(r.ok);
        CHECK((r.violation == Violation::Torque || r.violation == Violation::TorqueRate));
    }
    SECTION("a reference that does not end at rest: Timeout at T + 5 s")
    {
        ExperimentSpec e = spec("loops");
        e.speed = {0.05, 2.0, 1.0, 0.0};  // stops abruptly at 3 s ...
        e.gains.kt = 0.01;                // ... and the tangential law takes 100 s to brake
        e.gains.kiEta = 0.0;              // (without the integral action, which would brake it)
        e.duration = 3.0;
        const ExperimentSetup s(ex, e);
        const PreflightReport r = preflight(s, nullptr, s.monitor());
        CHECK_FALSE(r.ok);
        CHECK(r.violation == Violation::Timeout);
        CHECK(r.tFail > s.timeout());
        CHECK(r.tFail < s.timeout() + 2e-3);
    }
    SECTION("a tighter monitor fails where the hardware bounds pass")
    {
        const ExperimentSetup s(ex, spec("paper"));
        SafetyLimits tight;
        tight.ddqFrac = 0.05;  // the paper's run peaks at about 10 %
        const PreflightReport r = preflight(s, nullptr, SafetyMonitor(tight));
        CHECK_FALSE(r.ok);
        CHECK(r.violation == Violation::Acceleration);
        CHECK(r.limits.ddqFrac == 0.05);
    }
}

TEST_CASE("preflight from a given start: the measured start after the approach")
{
    // The robot driver runs the preflight again from the configuration the arm reached (apps/robot/robot_driver.cpp)
    const PandaExample& ex = pandaExample();
    const ExperimentSetup s(ex, spec("self-motion"));
    const PreflightReport base = preflight(s, nullptr, s.monitor());
    REQUIRE(base.ok);
    // The setup's own start: the same report
    const Vec7 q0 = s.start();
    const PreflightReport same = preflight(s, nullptr, s.monitor(), nullptr, &q0);
    CHECK(same.ok);
    CHECK(same.peaks.maxBeta == base.peaks.maxBeta);
    CHECK(same.simulated == base.simulated);
    // A start 3e-4 rad off on every joint (a realistic approach): passes, from that start
    const Vec7 near = q0 + Vec7::Constant(3e-4);
    const PreflightReport n = preflight(s, nullptr, s.monitor(), nullptr, &near);
    CHECK(n.ok);
    CHECK(n.startBeta != base.startBeta);
    // A start beyond the margin: the same pose with beta = 0.055 instead of 0.04 (Redundancy at t = 0)
    const Vec7 far = ex.start(ex.redundancy(), Vec5::Zero(), 0.055);
    const PreflightReport f = preflight(s, nullptr, s.monitor(), nullptr, &far);
    CHECK_FALSE(f.ok);
    CHECK(f.violation == Violation::Redundancy);
    CHECK(f.tFail == 0.0);
    // A non-finite start fails closed
    Vec7 bad = q0;
    bad[3] = std::nan("");
    const PreflightReport b = preflight(s, nullptr, s.monitor(), nullptr, &bad);
    CHECK_FALSE(b.ok);
    CHECK(b.violation == Violation::NonFinite);
}
