// The paper's rigid-body example end to end: the setup, and the closed loops of the design doc's Closed loops
// table. The expected values come from the reference implementation (pf_sim rigid-body), which this
// implementation matches to its printed precision.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "se3pf/rigid_body_example.hpp"
#include "se3pf/se3.hpp"
#include "test_helpers.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

using namespace Eigen;
using namespace se3pf;
using Catch::Matchers::WithinAbs;
using testing::maxDiff;
using testing::sci;
using Reference = RigidBodyExample::Reference;

// A 20 s closed loop takes about 16 s in a Debug build and 0.1 s in Release. The tests tagged SLOW_IN_DEBUG
// (test_helpers.hpp) are hidden in Debug and run by default in Release. The default Debug run keeps the 20 s pose
// loop, and the first second of the velocity loop, which holds that run's largest deviation and largest input.

namespace {
    // Whether ex.tube() compiles for an expression ex of type T (not for a temporary: deleted &&-overload)
    template <class T, class = void>
    struct HasTube : std::false_type {};
    template <class T>
    struct HasTube<T, std::void_t<decltype(std::declval<T>().tube())>> : std::true_type {};

    // Built once and shared by the test cases below: the example (its tube) and the runs, each computed the first
    // time a test needs it.
    const RigidBodyExample& example()
    {
        static const RigidBodyExample ex;
        return ex;
    }

    const RigidBodyRun& poseRun()
    {
        static const RigidBodyRun run = simulate(example(), Reference::Position);
        return run;
    }

    const RigidBodyRun& velocityRun()
    {
        static const RigidBodyRun run = simulate(example(), Reference::Velocity);
        return run;
    }

    // The velocity reference for 1 s only
    const RigidBodyRun& velocityFirstSecond()
    {
        static const RigidBodyRun run = simulate(example(), Reference::Velocity, 1.0);
        return run;
    }

    // The summary numbers again, from the log itself
    void checkLog(const RigidBodyRun& run, double T)
    {
        const std::size_t rows = static_cast<std::size_t>(std::lround(T / run.dt)) + 1;
        REQUIRE(run.log.size() == rows);
        double maxDev = 0.0, maxU = 0.0;
        for (std::size_t k = 0; k < run.log.size(); ++k) {
            const RigidBodyLogRow& r = run.log[k];
            REQUIRE(r.t == k * run.dt);
            REQUIRE(r.xiNorm == r.xi.norm());
            REQUIRE(r.xiDesign == RigidBodyExample::designedXiNorm(r.t));
            maxDev = std::max(maxDev, std::abs(r.xiNorm - r.xiDesign));
            maxU = std::max(maxU, r.u.cwiseAbs().maxCoeff());
        }
        CHECK(run.maxDeviation == maxDev);
        CHECK(run.maxInput == maxU);
        CHECK(run.finalXiNorm == run.log.back().xiNorm);
        CHECK(run.finalTanError == run.log.back().eTan);
    }
}

TEST_CASE("RigidBodyExample: the paper's setup (Appendix: Rigid body)")
{
    const RigidBodyExample& ex = example();
    const double L = 2.0 * std::sqrt(3.0) * EIGEN_PI;
    CHECK_THAT(ex.path().length(), WithinAbs(L, 1e-8));
    CHECK_THAT(ex.tube().length(), WithinAbs(L, 1e-8));
    CHECK(&ex.tube().path() == &ex.path());  // the tube refers to the example's own path
    CHECK(ex.metric().matrix() == Mat6::Identity());
    CHECK_THAT(ex.tube().diagnostics().angles[0], WithinAbs(2.166040, 1e-5));  // the rigid-body holonomy

    const Mat6 M = ex.body().mass();
    CHECK(M.diagonal() == (Vec6() << 1.0, 1.2, 0.8, 2.0, 2.0, 2.0).finished());

    // Start: Psi(0.1, xi0) is 25.45 deg and 0.4389 m from gamma(0.1); its outputs are (0.1, xi0)
    const Pose g0 = ex.startPose();
    double deg = 0.0, metres = 0.0;
    testing::poseError(ex.path().pose(RigidBodyExample::kS0), g0, deg, metres);
    CHECK_THAT(deg, WithinAbs(25.45, 0.005));
    CHECK_THAT(metres, WithinAbs(0.4389, 5e-5));
    CHECK_THAT(RigidBodyExample::xi0().norm(), WithinAbs(0.6245, 5e-5));
    double eta = -1.0;
    Vec5 xi = Vec5::Zero();
    REQUIRE(ex.tube().outputs(g0, 0.1, eta, xi));
    CHECK_THAT(eta, WithinAbs(0.1, 1e-9));
    CHECK(maxDiff(xi, RigidBodyExample::xi0()) < 1e-9);

    // Gains: lambda = 4, k_t = 5, k_p = 9, k_d = 6
    const Gains k = RigidBodyExample::gains();
    CHECK(k.Kp == 16.0);
    CHECK(k.Kd == 8.0);
    CHECK(k.kt == 5.0);
    CHECK(k.kpEta == 9.0);
    CHECK(k.kdEta == 6.0);

    // Designed response from rest: |xi0| at t = 0, flat start, (1 + 4t) e^{-4t}
    CHECK(RigidBodyExample::designedXiNorm(0.0) == RigidBodyExample::xi0().norm());
    CHECK_THAT(RigidBodyExample::designedXiNorm(0.5), WithinAbs(0.6245 * 3.0 * std::exp(-2.0), 1e-4));
    CHECK(std::abs(RigidBodyExample::designedXiNorm(1e-4) - RigidBodyExample::designedXiNorm(0.0)) < 1e-7);

    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<RigidBodyExample>);
    STATIC_REQUIRE_FALSE(std::is_move_constructible_v<RigidBodyExample>);
    STATIC_REQUIRE_FALSE(std::is_copy_assignable_v<RigidBodyExample>);
    STATIC_REQUIRE(HasTube<const RigidBodyExample&>::value);
    STATIC_REQUIRE_FALSE(HasTube<RigidBodyExample>::value);  // no tube from a temporary example
}

TEST_CASE("RigidBodyExample: the two references and their derivatives")
{
    const double L = example().tube().length(), h = 1e-5;

    const TangentialReference pose = example().reference(Reference::Position);
    REQUIRE(pose.mode == TangentialReference::Mode::Position);
    REQUIRE(pose.valid());
    auto etaD = [&](double t, double& e, double& de, double& dde) { pose.position(t, e, de, dde); };
    double e = 0.0, de = 0.0, dde = 0.0;
    etaD(0.0, e, de, dde);
    CHECK_THAT(e, WithinAbs(L / 4.0, 1e-15));  // eta_d(0) = L/4, at rest
    CHECK(de == 0.0);
    etaD(10.0, e, de, dde);
    CHECK_THAT(e, WithinAbs(L / 4.0 + 2.0 * L / 5.0, 1e-14));  // the far end, at t = 10
    for (double t : {1.3, 7.7, 16.0}) {
        double ep, dep, ddep, em, dem, ddem;
        etaD(t + h, ep, dep, ddep);
        etaD(t - h, em, dem, ddem);
        etaD(t, e, de, dde);
        CAPTURE(t);
        CHECK_THAT(de, WithinAbs((ep - em) / (2.0 * h), 1e-8));
        CHECK_THAT(dde, WithinAbs((dep - dem) / (2.0 * h), 1e-8));
    }

    const TangentialReference vel = example().reference(Reference::Velocity);
    REQUIRE(vel.mode == TangentialReference::Mode::Velocity);
    REQUIRE(vel.valid());
    double nu = 0.0, dnuDt = 0.0, dnuDeta = 1.0;
    vel.velocity(0.0, 3.0, nu, dnuDt, dnuDeta);
    CHECK(nu == 0.5);
    CHECK(dnuDeta == 0.0);
    for (double t : {1.3, 7.7, 16.0}) {
        double nup, num, a, b;
        vel.velocity(t + h, 2.0, nup, a, b);
        vel.velocity(t - h, 2.0, num, a, b);
        vel.velocity(t, 2.0, nu, dnuDt, dnuDeta);
        CAPTURE(t);
        CHECK_THAT(nu, WithinAbs(0.5 + 0.2 * std::sin(0.5 * t), 1e-15));
        CHECK_THAT(dnuDt, WithinAbs((nup - num) / (2.0 * h), 1e-9));
    }
}

TEST_CASE("simulate: the pose reference, 20 s (Closed loops: 9.69e-4, 1.8e-6, -1.7e-6, 29.1)")
{
    const RigidBodyRun& run = poseRun();
    REQUIRE(run.ok);
    CHECK(run.failureTime == -1.0);
    CHECK(run.reference == Reference::Position);
    checkLog(run, 20.0);

    CHECK_THAT(run.startAttitudeDeg, WithinAbs(25.45, 0.005));
    CHECK_THAT(run.startPosition, WithinAbs(0.4389, 5e-5));

    // |xi| follows the designed response |xi0| (1 + 4t) e^{-4t} to within 1e-3 (the paper), and the deviation is
    // the reference's 9.69e-4: a wrong term in A, Adot qdot or u = M qdd* + C zeta changes it.
    CAPTURE(sci(run.maxDeviation), run.maxDeviationTime, sci(run.finalXiNorm), sci(run.finalTanError), run.maxInput);
    CHECK(run.maxDeviation < 1e-3);
    CHECK(std::abs(run.maxDeviation - 9.69e-4) < 0.05 * 9.69e-4);
    CHECK_THAT(run.maxDeviationTime, WithinAbs(0.43, 0.011));
    CHECK_THAT(run.finalXiNorm, WithinAbs(1.78e-6, 0.05e-6));
    CHECK_THAT(run.finalTanError, WithinAbs(-1.69e-6, 0.05e-6));
    CHECK_THAT(run.maxInput, WithinAbs(29.11, 0.01));
    CHECK(run.maxInputTime == 0.0);
    CHECK(run.maxIterations == 1);  // the prediction is good enough for one Newton update per cycle

    // The first row: at rest at the start, eta_d(0) = L/4
    const RigidBodyLogRow& first = run.log.front();
    const double L = example().tube().length();
    CHECK_THAT(first.eta, WithinAbs(0.1, 1e-9));
    CHECK(first.etaDot == 0.0);
    CHECK_THAT(first.ref, WithinAbs(L / 4.0, 1e-14));
    CHECK_THAT(first.eTan, WithinAbs(0.1 - L / 4.0, 1e-9));
    CHECK(maxDiff(first.xi, RigidBodyExample::xi0()) < 1e-9);

    // At the end eta_d(20) = L/4 again, and eta has followed it
    CHECK_THAT(run.log.back().eta, WithinAbs(L / 4.0, 1e-5));
}

TEST_CASE("simulate: the velocity reference, 20 s (Closed loops: 9.51e-4, 2.8e-6, -3.2e-6, 14.5)", SLOW_IN_DEBUG)
{
    const RigidBodyRun& run = velocityRun();
    REQUIRE(run.ok);
    CHECK(run.reference == Reference::Velocity);
    checkLog(run, 20.0);

    CAPTURE(sci(run.maxDeviation), run.maxDeviationTime, sci(run.finalXiNorm), sci(run.finalTanError), run.maxInput);
    CHECK(run.maxDeviation < 1e-3);
    CHECK(std::abs(run.maxDeviation - 9.51e-4) < 0.05 * 9.51e-4);
    CHECK_THAT(run.maxDeviationTime, WithinAbs(0.408, 0.011));
    CHECK_THAT(run.finalXiNorm, WithinAbs(2.81e-6, 0.05e-6));
    CHECK_THAT(run.finalTanError, WithinAbs(-3.25e-6, 0.05e-6));
    CHECK_THAT(run.maxInput, WithinAbs(14.48, 0.01));
    CHECK(run.maxInputTime == 0.0);
    CHECK(run.maxIterations == 1);

    // e_v = eta_dot - nu_d starts at -0.5 (at rest) and ends near zero; nu_d(20) = 0.5 + 0.2 sin(10)
    CHECK_THAT(run.log.front().eTan, WithinAbs(-0.5, 1e-15));
    CHECK_THAT(run.log.back().ref, WithinAbs(0.5 + 0.2 * std::sin(10.0), 1e-14));
    CHECK_THAT(run.log.back().etaDot, WithinAbs(run.log.back().ref, 1e-5));

    // Its first second is the 1 s run that stands in for it in Debug, row for row
    const RigidBodyRun& first = velocityFirstSecond();
    REQUIRE(run.log.size() > first.log.size());
    double diff = 0.0;
    for (std::size_t k = 0; k < first.log.size(); ++k) {
        const RigidBodyLogRow& a = run.log[k];
        const RigidBodyLogRow& b = first.log[k];
        diff = std::max({diff, std::abs(a.eta - b.eta), std::abs(a.etaDot - b.etaDot), maxDiff(a.xi, b.xi),
                         maxDiff(a.u, b.u)});
    }
    CHECK(diff == 0.0);
}

TEST_CASE("simulate: the velocity reference, first second (Closed loops: 9.51e-4 at t = 0.41 s, 14.5 at t = 0)")
{
    // The part of the 20 s run that holds its largest deviation and its largest input; in Debug it stands in for
    // the 20 s run, which is hidden there.
    const RigidBodyRun& run = velocityFirstSecond();
    REQUIRE(run.ok);
    CHECK(run.reference == Reference::Velocity);
    checkLog(run, 1.0);

    CAPTURE(sci(run.maxDeviation), run.maxDeviationTime, run.maxInput);
    CHECK(run.maxDeviation < 1e-3);
    CHECK(std::abs(run.maxDeviation - 9.51e-4) < 0.05 * 9.51e-4);
    CHECK_THAT(run.maxDeviationTime, WithinAbs(0.408, 0.011));
    CHECK_THAT(run.maxInput, WithinAbs(14.48, 0.01));
    CHECK(run.maxInputTime == 0.0);
    CHECK(run.maxIterations == 1);
    CHECK_THAT(run.log.front().eTan, WithinAbs(-0.5, 1e-15));
}

TEST_CASE("simulate: the substeps only refine the plant's integration: 1, 4 and 16 substeps agree")
{
    // u is held over each control period whatever the number of substeps, and RK4 with substeps of length
    // h = dt / substeps is accurate to O(h^4), so after 1 s the three runs agree to about 1e-12. A substep of the
    // wrong length would integrate the wrong time per period and change the run at O(1).
    const RigidBodyRun& four = velocityFirstSecond();
    REQUIRE(four.ok);
    CHECK(four.substeps == 4);
    const RigidBodyLogRow& b = four.log.back();
    CHECK(b.eta > 0.5);  // it has moved: eta went from 0.1 to about 0.55
    for (int substeps : {1, 16}) {
        const RigidBodyRun run = simulate(example(), Reference::Velocity, 1.0, 0.002, substeps);
        REQUIRE(run.ok);
        CHECK(run.substeps == substeps);
        REQUIRE(run.log.size() == four.log.size());
        const RigidBodyLogRow& a = run.log.back();
        const double errEta = std::abs(a.eta - b.eta);
        const double errEtaDot = std::abs(a.etaDot - b.etaDot);
        const double errXi = maxDiff(a.xi, b.xi);
        CAPTURE(substeps, a.eta, b.eta, sci(errEta), sci(errEtaDot), sci(errXi));
        CHECK(errEta < 1e-10);
        CHECK(errEtaDot < 1e-10);
        CHECK(errXi < 1e-10);
    }
}

TEST_CASE("simulate: the deviation is a sample-and-hold effect: halving dt halves it", SLOW_IN_DEBUG)
{
    // Design doc, Three kinds of test: 9.51e-4 at 2 ms, 4.75e-4 at 1 ms. If halving dt does not halve the
    // deviation, the controller has a bug that the sampling does not explain.
    const RigidBodyRun fine = simulate(example(), Reference::Velocity, 20.0, 0.001, 4);
    REQUIRE(fine.ok);
    checkLog(fine, 20.0);
    const double ratio = velocityRun().maxDeviation / fine.maxDeviation;
    CAPTURE(sci(fine.maxDeviation), ratio);
    CHECK(std::abs(fine.maxDeviation - 4.75e-4) < 0.05 * 4.75e-4);
    CHECK(ratio > 1.8);
    CHECK(ratio < 2.2);
    CHECK_THAT(fine.maxInput, WithinAbs(14.48, 0.01));  // the first command does not depend on dt
}

TEST_CASE("simulate: argument checks; T = 0 gives the first command only")
{
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(simulate(example(), Reference::Position, -1.0), std::invalid_argument);
    CHECK_THROWS_AS(simulate(example(), Reference::Position, 1.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(simulate(example(), Reference::Position, 1.0, 0.002, 0), std::invalid_argument);
    CHECK_THROWS_AS(simulate(example(), Reference::Position, nan), std::invalid_argument);
    CHECK_THROWS_AS(simulate(example(), Reference::Velocity, 1.0, nan), std::invalid_argument);
    // Not finite, or more than kMaxSteps periods: an invalid_argument before any allocation (not a length_error
    // from reserving 5e17 rows, or an undefined cast of inf to long)
    STATIC_REQUIRE(RigidBodyExample::kMaxSteps == 10'000'000);
    CHECK_THROWS_AS(simulate(example(), Reference::Position, inf), std::invalid_argument);
    CHECK_THROWS_AS(simulate(example(), Reference::Position, 1e15), std::invalid_argument);
    CHECK_THROWS_AS(simulate(example(), Reference::Velocity, 1.0, 1e-320), std::invalid_argument);

    const RigidBodyRun one = simulate(example(), Reference::Position, 0.0);
    REQUIRE(one.ok);
    REQUIRE(one.log.size() == 1);
    // The first command at rest: u = M qdd* (C zeta = 0), the same as in the 20 s run.
    CHECK(maxDiff(one.log[0].u, poseRun().log[0].u) == 0.0);
    CHECK_THAT(one.maxInput, WithinAbs(29.11, 0.01));
}
