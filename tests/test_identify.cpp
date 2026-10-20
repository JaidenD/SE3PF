// Tests of the joint identification (se3pf/identify.hpp; CTRL-1): the pulse sequence on a simulated joint with
// Coulomb friction and a known inertia, the breakaway and inertia recovered by the fit, the limits of every pulse,
// the blend that keeps the torque rate inside its margin, and the refusals.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "se3pf/experiments.hpp"
#include "se3pf/identify.hpp"
#include "se3pf/panda.hpp"
#include "test_helpers.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace se3pf;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
    // A 7-joint plant with diagonal inertia M, Coulomb friction Fc (regularized below vEps) and no gravity, held by
    // the identification's own hold torque: semi-implicit Euler at 1 kHz, as the mock robot. Returns the log.
    struct SimResult {
        std::vector<IdentifyLogRow> rows;
        double maxRate = 0.0;   // max |d tau / dt| of the command, N m/s
        double maxTravel = 0.0; // max |q_j - q_j at the pulse start| over every pulse, rad
        double maxSpeed = 0.0;
        bool done = false;
        double T = 0.0;
    };
    SimResult simulate(JointIdentification& ident, const Vec7& q0, const Vec7& M, const Vec7& Fc, double maxT)
    {
        SimResult r;
        Vec7 q = q0, dq = Vec7::Zero(), tauPrev = Vec7::Zero();
        const double dt = 1e-3, vEps = 1e-4;
        double qPulse0 = 0.0;
        int lastJoint = -1;
        for (long k = 0; k * dt < maxT; ++k) {
            const double t = k * dt;
            bool done = false;
            const Vec7 tau = ident.step(t, q, dq, done);
            if (k > 0)
                r.maxRate = std::max(r.maxRate, (tau - tauPrev).cwiseAbs().maxCoeff() / dt);
            tauPrev = tau;
            IdentifyLogRow row;
            row.t = t;
            row.phase = static_cast<int>(ident.phase());
            row.joint = ident.activeJoint() + 1;
            row.direction = ident.direction();
            row.pulse = ident.pulse();
            row.pulseTorque = ident.pulseTorque();
            row.tau = tau;
            row.q = q;
            row.dq = dq;
            row.tauJ = tau;  // a sensor that sees the command (no gravity here)
            r.rows.push_back(row);
            if (ident.activeJoint() != lastJoint || ident.phase() == IdentifyPhase::Hold) {
                lastJoint = ident.activeJoint();
                if (lastJoint >= 0)
                    qPulse0 = q[lastJoint];
            }
            if (lastJoint >= 0 && ident.phase() != IdentifyPhase::Hold) {
                r.maxTravel = std::max(r.maxTravel, std::abs(q[lastJoint] - qPulse0));
                r.maxSpeed = std::max(r.maxSpeed, std::abs(dq[lastJoint]));
            }
            r.T = t;
            if (done) {
                r.done = true;
                break;
            }
            // The plant: M qdd = tau - Fc sign(dq), implicit in the friction (as mock/franka/robot.h)
            for (int i = 0; i < 7; ++i) {
                const double D = Fc[i] / std::max(std::abs(dq[i]), vEps);
                const double ddq = (tau[i] - D * dq[i]) / (M[i] + dt * D);
                dq[i] += dt * ddq;
                q[i] += dt * dq[i];
            }
        }
        return r;
    }
}  // namespace

TEST_CASE("JointIdentification: breakaway and inertia of a simulated joint, within the pulse limits")
{
    IdentifySettings s;
    const Vec7 q0 = PandaArm::restPosture();
    const Vec7 M = (Vec7() << 1.0, 1.7, 1.1, 1.0, 0.25, 0.2, 0.15).finished();
    const Vec7 Fc = (Vec7() << 0.8, 1.2, 0.7, 0.6, 0.7, 0.4, 0.6).finished();
    JointIdentification ident(s, {4, 6, 1}, q0);  // joints 5, 7 and 2
    CHECK(ident.duration() > 10.0);
    const SimResult r = simulate(ident, q0, M, Fc, ident.duration() + 5.0);
    REQUIRE(r.done);
    CHECK(r.T < ident.duration());
    CAPTURE(r.maxRate, r.maxTravel, r.maxSpeed);
    CHECK(r.maxRate < 0.6 * PandaArm::tauRateMax);  // the blend keeps the command inside the torque-rate margin
    CHECK(r.maxRate > 20.0);                         // ... but the steps are steps
    CHECK(r.maxTravel <= s.travelLimit + 0.02);      // a pulse ends at the travel limit (plus one period and the blend)
    CHECK(r.maxSpeed <= s.speedLimit + 0.1);
    // Every joint ends held near where it started, the hold re-latched after each pulse
    const IdentifyLogRow& last = r.rows.back();
    CHECK(last.phase == 0);
    CHECK((last.q - q0).cwiseAbs().maxCoeff() < 6.0 * s.travelLimit);
    CHECK(last.dq.norm() < 0.02);

    const std::vector<JointIdentificationResult> fit = fitIdentification(r.rows, s);
    REQUIRE(fit.size() == 7);
    for (int j : {4, 6, 1}) {
        CAPTURE(j);
        for (int d = 0; d < 2; ++d) {
            CAPTURE(d);
            REQUIRE(fit[j].found[d]);
            // The raw breakaway is the ramp's value when the joint reached motionSpeed: above Fc by the lag
            // sqrt(2 M v_m r) at most (plus a period); the corrected one is Fc
            const double lag = std::sqrt(2.0 * M[j] * s.motionSpeed * s.rampRate);
            CHECK(fit[j].breakawayRaw[d] >= Fc[j] - 0.01);
            CHECK(fit[j].breakawayRaw[d] <= Fc[j] + lag + 0.03);
            CHECK(std::abs(ident.breakaway(j, d) - fit[j].breakawayRaw[d]) < s.rampRate * (s.motionHold + 0.002));
            REQUIRE(fit[j].inertiaFound[d]);
            CHECK_THAT(fit[j].breakawayCmd[d], WithinAbs(Fc[j], 0.05));
            CHECK_THAT(fit[j].inertiaCmd[d], WithinRel(M[j], 0.15));
            CHECK_THAT(fit[j].inertiaSensor[d], WithinRel(M[j], 0.15));
            CHECK_THAT(fit[j].kineticFriction[d], WithinAbs(Fc[j], 0.15));
        }
    }
    // Joints not in the list: nothing
    for (int j : {0, 2, 3, 5})
        CHECK_FALSE((fit[j].found[0] || fit[j].found[1] || fit[j].inertiaFound[0]));
}

TEST_CASE("JointIdentification: a joint whose friction exceeds the cap is reported as not found, the steps skipped")
{
    IdentifySettings s;
    const Vec7 q0 = PandaArm::restPosture();
    Vec7 M = Vec7::Constant(0.5), Fc = Vec7::Constant(0.3);
    Fc[6] = 5.0;  // above the 1.5 N m cap of joint 7
    JointIdentification ident(s, {6}, q0);
    const SimResult r = simulate(ident, q0, M, Fc, ident.duration() + 5.0);
    REQUIRE(r.done);
    CHECK(std::isnan(ident.breakaway(6, 0)));
    CHECK(std::isnan(ident.breakaway(6, 1)));
    CHECK(r.maxTravel < 1e-3);  // it never moved
    double maxPulse = 0.0;
    int steps = 0;
    for (const IdentifyLogRow& row : r.rows) {
        maxPulse = std::max(maxPulse, std::abs(row.pulseTorque));
        steps += row.phase == 2 ? 1 : 0;
    }
    CHECK_THAT(maxPulse, WithinAbs(s.torqueCap[6], s.rampRate * 0.0015));  // the cap is reached within a period
    CHECK(steps == 0);
    const std::vector<JointIdentificationResult> fit = fitIdentification(r.rows, s);
    CHECK_FALSE(fit[6].found[0]);
    CHECK_FALSE(fit[6].inertiaFound[0]);
    CHECK(identificationSummary(fit).find("not found") != std::string::npos);
}

TEST_CASE("JointIdentification: refusals, non-finite measurements, the hold law")
{
    IdentifySettings s;
    const Vec7 q0 = PandaArm::restPosture();
    CHECK_THROWS_AS(JointIdentification(s, {}, q0), std::invalid_argument);
    CHECK_THROWS_AS(JointIdentification(s, {7}, q0), std::invalid_argument);
    CHECK_THROWS_AS(JointIdentification(s, {-1}, q0), std::invalid_argument);
    IdentifySettings bad = s;
    bad.rampRate = 0.0;
    CHECK_THROWS_AS(JointIdentification(bad, {0}, q0), std::invalid_argument);
    bad = s;
    bad.stepExtra[1] = bad.stepExtra[0];
    CHECK_THROWS_AS(JointIdentification(bad, {0}, q0), std::invalid_argument);
    bad = s;
    bad.torqueCap[3] = -1.0;
    CHECK_THROWS_AS(JointIdentification(bad, {0}, q0), std::invalid_argument);
    Vec7 qn = q0;
    qn[2] = std::nan("");
    CHECK_THROWS_AS(JointIdentification(s, {0}, qn), std::invalid_argument);

    // The first command is zero (the hold is blended in from zero over s.blend); during the settle the command is
    // the damped hold at q0; a NaN measurement gives a finite hold command
    JointIdentification ident(s, {0}, q0);
    bool done = false;
    const Vec7 dq = Vec7::Constant(0.01);
    CHECK(ident.step(0.0, q0, dq, done).isZero());
    CHECK(testing::maxDiff(ident.step(s.blend, q0, dq, done), SafetyStop::torque(q0, dq, q0)) < 1e-12);
    CHECK_FALSE(done);
    CHECK(ident.step(0.1, qn, dq, done).allFinite());
    CHECK(ident.step(0.2, q0, Vec7::Constant(std::nan("")), done).allFinite());
    // The ramp starts after the settle, from zero, at the rate, on the joint, in the + direction first
    Vec7 tau = ident.step(s.settle + 0.001, q0, Vec7::Zero(), done);
    CHECK(ident.phase() == IdentifyPhase::Ramp);
    CHECK(ident.activeJoint() == 0);
    CHECK(ident.direction() == 1);
    for (int k = 2; k < 200; ++k)
        tau = ident.step(s.settle + 0.001 * k, q0, Vec7::Zero(), done);
    CHECK(ident.pulseTorque() > 0.0);
    CHECK_THAT(ident.pulseTorque(), WithinAbs(s.rampRate * 0.198, 1e-9));  // the pulse began one period in
    CHECK_THAT(tau[0], WithinAbs(ident.pulseTorque(), 1e-9));  // after the blend, the command is the pulse
}
