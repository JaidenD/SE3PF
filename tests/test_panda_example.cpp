// The paper's Panda example: the setup and its start configurations, and the closed loops of the design doc's
// Closed loops table (Panda rows). The expected values of the closed loops come from the reference library (pf, in
// ~/Desktop/panda_path_following/include) with the paper's rest-posture lift condition added, run with this
// example's gains (16/8, k_t = 5), no ramp and the same held-acceleration loop: the oracle program
// tools/oracle/panda_oracle.cpp (its header says how to build and run it). This implementation matches it to 6-7
// significant digits. pf_sim's panda demos use other gains, a ramp on mu and other
// lifts, so they do not reproduce these numbers.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "se3pf/panda_example.hpp"
#include "se3pf/se3.hpp"
#include "panda_fixture.hpp"
#include "test_helpers.hpp"

#include <algorithm>
#include <chrono>
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
using testing::poseError;
using testing::sci;
using RunKind = PandaExample::RunKind;

// A Panda control step takes about 2 ms in a Debug build and 25 us in Release, so the 20 s paper run (20001 steps)
// and the one-loop self-motion run (9627 steps) would take about a minute in Debug. The tests tagged SLOW_IN_DEBUG
// (test_helpers.hpp) are hidden in Debug and run by default in Release. The default Debug run keeps the first second
// of the paper run and the first half second of the self-motion run, which hold their largest deviations from the
// designed responses.

namespace {
    // Built once for the whole test program (panda_fixture.hpp)
    const PandaExample& example() { return testing::pandaExample(); }

    // The closed loops, each computed the first time a test needs it
    const PandaRun& paperRun()
    {
        static const PandaRun run = simulate(example(), RunKind::Paper);
        return run;
    }

    const PandaRun& paperFirstSecond()
    {
        static const PandaRun run = simulate(example(), RunKind::Paper, 1.0);
        return run;
    }

    const PandaRun& selfMotionRun()
    {
        static const PandaRun run = simulate(example(), RunKind::SelfMotion);
        return run;
    }

    const PandaRun& selfMotionFirstHalfSecond()
    {
        static const PandaRun run = simulate(example(), RunKind::SelfMotion, 0.5);
        return run;
    }

    // The log has one row per control instant, and the summary numbers are the extremes of its columns. A run that
    // stopped early (ok = false) is checked over the rows it has.
    void checkLog(const PandaRun& run, std::size_t rows)
    {
        REQUIRE(run.log.size() == rows);
        // The log was reserved for the whole run before the loop, so it never grew inside it. Grown row by row it
        // would end with a power-of-two capacity (1024 for 1001 rows, 32768 for 20001).
        if (run.ok)
            CHECK(run.log.capacity() == run.log.size());
        double maxXi = 0.0, maxDev = 0.0, maxBeta = 0.0, maxBetaDev = 0.0, maxVel = 0.0, maxAcc = 0.0;
        double minMargin = std::numeric_limits<double>::infinity(), minSigma = minMargin;
        int maxIt = 0;
        const Vec7 dqMax = PandaArm::dqMax(), ddqMax = PandaArm::ddqMax();
        for (std::size_t k = 0; k < run.log.size(); ++k) {
            const PandaLogRow& r = run.log[k];
            REQUIRE(r.t == k * run.dt);
            REQUIRE(r.xiNorm == r.xi.norm());
            REQUIRE(r.eTan == r.etaDot - r.ref);
            REQUIRE(r.jointMargin == PandaArm::jointMargin(r.q));
            maxXi = std::max(maxXi, r.xiNorm);
            maxDev = std::max(maxDev, std::abs(r.xiNorm - r.xiDesign));
            maxBeta = std::max(maxBeta, std::abs(r.beta));
            maxBetaDev = std::max(maxBetaDev, std::abs(r.beta - r.betaDesign));
            maxVel = std::max(maxVel, (r.qd.cwiseAbs().array() / dqMax.array()).maxCoeff());
            maxAcc = std::max(maxAcc, (r.qdd.cwiseAbs().array() / ddqMax.array()).maxCoeff());
            minMargin = std::min(minMargin, r.jointMargin);
            minSigma = std::min(minSigma, r.sigmaMin);
            maxIt = std::max(maxIt, r.iterations);
        }
        CHECK(run.maxXiNorm == maxXi);
        CHECK(run.maxXiDeviation == maxDev);
        CHECK(run.maxBeta == maxBeta);
        CHECK(run.maxBetaDeviation == maxBetaDev);
        CHECK(run.maxVelocityFraction == maxVel);
        CHECK(run.maxAccelerationFraction == maxAcc);
        CHECK(run.minJointMargin == minMargin);
        CHECK(run.minSigma == minSigma);
        CHECK(run.maxIterations == maxIt);
        CHECK(run.startXiNorm == run.log.front().xiNorm);
        CHECK(run.startBeta == run.log.front().beta);
        CHECK(run.finalXiNorm == run.log.back().xiNorm);
        CHECK(run.finalEtaDot == run.log.back().etaDot);
        CHECK(run.finalBeta == run.log.back().beta);

        // Consecutive rows follow the held-acceleration update: q_k+1 = q_k + dt qd_k + dt^2/2 qdd_k, and
        // qd_k+1 = qd_k + dt qdd_k. Relative to the size of the values (at least 1), since an unstable run (see the
        // failure test) reaches thousands of radians, where rounding alone is of order 1e-12.
        double errQ = 0.0, errQd = 0.0;
        for (std::size_t k = 0; k + 1 < run.log.size(); ++k) {
            const PandaLogRow& a = run.log[k];
            const PandaLogRow& b = run.log[k + 1];
            const double dt = run.dt;
            const double scaleQ = std::max(1.0, b.q.cwiseAbs().maxCoeff());
            const double scaleQd = std::max(1.0, b.qd.cwiseAbs().maxCoeff());
            errQ = std::max(errQ, maxDiff(b.q, Vec7(a.q + dt * a.qd + (0.5 * dt * dt) * a.qdd)) / scaleQ);
            errQd = std::max(errQd, maxDiff(b.qd, Vec7(a.qd + dt * a.qdd)) / scaleQd);
        }
        CAPTURE(sci(errQ), sci(errQd));
        CHECK(errQ < 1e-14);  // rounding only; a wrong update (no dt^2/2 term, say) is off by about 1e-6
        CHECK(errQd < 1e-14);
    }

    // Two runs agree row for row on their common rows (the shorter one is the start of the longer one).
    double rowDifference(const PandaRun& a, const PandaRun& b)
    {
        const std::size_t n = std::min(a.log.size(), b.log.size());
        double diff = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            const PandaLogRow& x = a.log[k];
            const PandaLogRow& y = b.log[k];
            diff = std::max({diff, std::abs(x.eta - y.eta), std::abs(x.etaDot - y.etaDot), maxDiff(x.xi, y.xi),
                             std::abs(x.beta - y.beta), maxDiff(x.q, y.q), maxDiff(x.qd, y.qd),
                             maxDiff(x.qdd, y.qdd)});
        }
        return diff;
    }

    // eta, xi and beta of a configuration, as the controller computes them (eq:A)
    struct Outputs {
        double eta = 0.0;
        Vec5 xi = Vec5::Zero();
        double beta = 0.0;
    };

    Outputs outputsAt(const Vec7& q, const RedundancyOutput& out)
    {
        const PandaExample& ex = example();
        Mat67 Jb = Mat67::Zero();
        const Pose g = PandaArm::output(q, &Jb);
        Outputs o;
        REQUIRE(ex.tube().outputs(g, PandaExample::s0, o.eta, o.xi));
        const Mat67 B = ex.tube().JPsi(o.eta, o.xi).partialPivLu().solve(Jb);
        Vec7 dbeta = Vec7::Zero();
        out.evaluate(q, o.eta, Vec7(B.row(0).transpose()), o.beta, dbeta);
        return o;
    }
}

namespace {
    // Whether ex.tube() and ex.redundancy() compile for an expression ex of type T (the deleted &&-overloads make them
    // fail on a temporary, whose parts a controller would keep references to)
    template <class T, class = void>
    struct HasTube : std::false_type {};
    template <class T>
    struct HasTube<T, std::void_t<decltype(std::declval<T>().tube())>> : std::true_type {};
    template <class T, class = void>
    struct HasRedundancy : std::false_type {};
    template <class T>
    struct HasRedundancy<T, std::void_t<decltype(std::declval<T>().redundancy())>> : std::true_type {};
}

TEST_CASE("PandaExample: cannot be copied or moved (its members refer to each other)")
{
    STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<PandaExample>);
    STATIC_REQUIRE_FALSE(std::is_move_constructible_v<PandaExample>);
    STATIC_REQUIRE_FALSE(std::is_copy_assignable_v<PandaExample>);
    CHECK_THROWS_AS(PandaExample(0.0), std::invalid_argument);

    // Its parts are not handed out by a temporary example
    STATIC_REQUIRE(HasTube<const PandaExample&>::value);
    STATIC_REQUIRE(HasRedundancy<const PandaExample&>::value);
    STATIC_REQUIRE_FALSE(HasTube<PandaExample>::value);
    STATIC_REQUIRE_FALSE(HasRedundancy<PandaExample&&>::value);
}

TEST_CASE("PandaExample: the paper's setup (Appendix: Panda)")
{
    const PandaExample& ex = example();

    // Metric l = 0.1 m and the coning circle: L = 0.77008765, gamma(0) at (0.52, 0, 0.45)
    CHECK(maxDiff(ex.metric().matrix(), Mat6(Vec6(0.01, 0.01, 0.01, 1.0, 1.0, 1.0).asDiagonal())) < 1e-16);
    CHECK_THAT(ex.path().length(), WithinAbs(0.77008765, 1e-8));
    CHECK_THAT(ex.tube().length(), WithinAbs(0.77008765, 1e-8));
    CHECK(maxDiff(Vector3d(ex.path().pose(0.0).block<3,1>(0,3)), Vector3d(0.52, 0.0, 0.45)) < 1e-12);
    CHECK_THAT(ex.lift().length(), WithinAbs(0.77008765, 1e-8));

    // R0 = the flange attitude at the rest posture: approach axis down
    const double c = std::sqrt(0.5);
    Matrix3d R0;
    R0 <<  c, -c,  0,
          -c, -c,  0,
           0,  0, -1;
    CHECK(maxDiff(ex.R0(), R0) < 1e-12);

    // The lift is the paper's (the full lift table is in test_lift.cpp)
    CHECK_THAT(ex.lift().diagnostics().minJointMargin, WithinAbs(0.3592, 1e-4));

    // Constants
    CHECK(PandaExample::s0 == 0.1);
    CHECK(maxDiff(PandaExample::xi0(), Vec5(0.02, -0.015, 0.01, -0.005, 0.015)) == 0.0);
    CHECK(PandaExample::nuD == 0.08);
    CHECK(PandaExample::lambda == 4.0);
    CHECK(PandaExample::kt == 5.0);
    CHECK(PandaExample::KpBeta == 16.0);
    CHECK(PandaExample::KdBeta == 8.0);
    CHECK(PandaExample::dt == 1e-3);
    CHECK(PandaExample::T == 20.0);
    CHECK(PandaExample::selfMotionBeta0 == 0.2);
}

TEST_CASE("PandaExample: the paper's start Psi(0.1, xi0) with beta = 0")
{
    const PandaExample& ex = example();
    const Vec7& q0 = ex.paperStart();

    // 15.07 deg and 16.82 mm from gamma(0.1) [15.1 deg, 16.8 mm], |xi0| = 0.03122
    double deg = 0.0, m = 0.0;
    poseError(ex.path().pose(0.1), PandaArm::output(q0, nullptr), deg, m);
    CHECK_THAT(deg, WithinAbs(15.07, 0.005));
    CHECK_THAT(m, WithinAbs(0.01682, 5e-6));

    const Outputs o = outputsAt(q0, ex.redundancy());
    const double errPose = maxDiff(PandaArm::output(q0, nullptr), ex.tube().Psi(0.1, PandaExample::xi0()));
    CAPTURE(sci(errPose), sci(o.beta), sci(o.eta - 0.1));
    CHECK(errPose < 1e-10);
    CHECK(std::abs(o.beta) < 1e-10);
    CHECK(std::abs(o.eta - 0.1) < 1e-10);
    CHECK(maxDiff(o.xi, PandaExample::xi0()) < 1e-10);
    CHECK_THAT(o.xi.norm(), WithinAbs(0.03122, 1e-5));

    // The reference implementation's start configuration (also found by an independent MuJoCo-side run)
    Vec7 expected;
    expected << 0.13435, -0.23612, 0.05396, -2.41230, 0.30751, 2.49181, 0.82647;
    const double errQ = maxDiff(q0, expected);
    CAPTURE(sci(errQ));
    CHECK(errQ < 1e-4);
    CHECK(PandaArm::jointMargin(q0) > 0.1);
}

TEST_CASE("PandaExample: the self-motion start gamma(0.1) with beta = 0.2")
{
    const PandaExample& ex = example();
    const Vec7& q0 = ex.selfMotionStart();
    const Outputs o = outputsAt(q0, ex.redundancy());
    const double errPose = maxDiff(PandaArm::output(q0, nullptr), ex.path().pose(0.1));
    CAPTURE(sci(errPose), sci(o.beta - 0.2), sci(o.xi.norm()));
    CHECK(errPose < 1e-10);
    CHECK(std::abs(o.beta - 0.2) < 1e-10);
    CHECK(o.xi.norm() < 1e-10);
    CHECK(PandaArm::jointMargin(q0) > 0.1);

    // On the self-motion through q*(0.1): beta = X^T (q - q*) = 0.2, so q is about 0.2 rad from q*(0.1)
    Vec7 qs = Vec7::Zero();
    ex.lift().eval(0.1, qs);
    CHECK_THAT((q0 - qs).norm(), WithinAbs(0.200018, 1e-5));
}

TEST_CASE("PandaExample: a start for the coupling output (comparison runs)")
{
    const PandaExample& ex = example();
    const Vec7 q = ex.start(ex.coupling(), PandaExample::xi0(), 0.0);
    CHECK(maxDiff(PandaArm::output(q, nullptr), ex.tube().Psi(0.1, PandaExample::xi0())) < 1e-10);
    CHECK(std::abs(q[2] - q[0]) < 1e-10);
    CHECK(PandaArm::jointMargin(q) > 0.1);
    // the reference implementation's coupled start
    Vec7 expected;
    expected << 0.09368941, -0.23662446, 0.09368941, -2.41144297, 0.32244724, 2.48902589, 0.81586552;
    CHECK(maxDiff(q, expected) < 1e-7);
}

TEST_CASE("PandaExample: construction time (Release)", "[.timing]")
{
    // Release build only: ./build-rel/tests "[.timing]". Tube, lift, kernel frame and the two starts (about
    // 0.03 s; the lift's minimum check is not part of it, see minimumCheck).
#ifndef NDEBUG
    SKIP("timing is meaningful in a Release build only");
#endif
    const auto t0 = std::chrono::steady_clock::now();
    const PandaExample ex;
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    WARN("PandaExample construction: " << seconds << " s");
    CHECK(ex.paperStart().allFinite());
    CHECK(seconds < 0.5);
}

// ===============================================================================================================
// The closed loops: design doc, Panda simulation, and the Panda rows of the Closed loops table

TEST_CASE("PandaExample: gains, reference, designed response and durations of the closed loops")
{
    // lambda = 4: K_p = K_p,beta = 16, K_d = K_d,beta = 8; k_t = 5
    const Gains k = PandaExample::gains();
    CHECK(k.Kp == 16.0);
    CHECK(k.Kd == 8.0);
    CHECK(k.KpBeta == 16.0);
    CHECK(k.KdBeta == 8.0);
    CHECK(k.kt == 5.0);

    // nu_d = 0.08, constant
    const TangentialReference ref = PandaExample::reference();
    REQUIRE(ref.mode == TangentialReference::Mode::Velocity);
    REQUIRE(ref.valid());
    double nu = 0.0, dnuDt = 1.0, dnuDeta = 1.0;
    ref.velocity(3.0, 0.4, nu, dnuDt, dnuDeta);
    CHECK(nu == 0.08);
    CHECK(dnuDt == 0.0);
    CHECK(dnuDeta == 0.0);

    // y0 (1 + 4t) e^{-4t}: y0 at t = 0 with zero slope, 3 y0 e^{-2} at t = 0.5
    CHECK(PandaExample::designedResponse(0.2, 0.0) == 0.2);
    CHECK_THAT(PandaExample::designedResponse(0.2, 0.5), WithinAbs(0.6 * std::exp(-2.0), 1e-15));
    CHECK(std::abs(PandaExample::designedResponse(0.2, 1e-4) - 0.2) < 2e-7);  // flat start: 0.2 * 8e-8
    CHECK(PandaExample::designedResponse(0.0, 1.0) == 0.0);

    // 20 s, and one loop L / nu_d = 9.6261 s
    const PandaExample& ex = example();
    CHECK(ex.duration(RunKind::Paper) == 20.0);
    CHECK_THAT(ex.duration(RunKind::SelfMotion), WithinAbs(0.77008765 / 0.08, 1e-6));
}

TEST_CASE("simulate (Panda): the paper's run, first second (Closed loops: 2.43e-5 at t = 0.41 s, 1.5e-5)")
{
    // The part of the 20 s run that holds its largest deviation from the designed response, its largest beta, its
    // largest joint velocity and acceleration and its smallest sigma_min(A); in Debug it stands in for the 20 s run.
    const PandaRun& run = paperFirstSecond();
    REQUIRE(run.ok);
    CHECK(run.failureTime == -1.0);
    CHECK(run.kind == RunKind::Paper);
    CHECK(run.T == 1.0);
    CHECK(run.dt == 1e-3);
    checkLog(run, 1001);

    // Start: at rest at Psi(0.1, xi0) with beta = 0, 15.07 deg and 16.82 mm from gamma(0.1) [15.1 deg, 16.8 mm]
    const PandaLogRow& first = run.log.front();
    CHECK_THAT(run.startAttitudeDeg, WithinAbs(15.07, 0.005));
    CHECK_THAT(run.startPosition, WithinAbs(0.01682, 5e-6));
    CHECK_THAT(run.startXiNorm, WithinAbs(0.03122, 1e-5));
    CHECK(std::abs(run.startBeta) < 1e-10);
    CHECK(first.q == example().paperStart());
    CHECK(first.qd.isZero(0.0));
    CHECK_THAT(first.eta, WithinAbs(0.1, 1e-10));
    CHECK(maxDiff(first.xi, PandaExample::xi0()) < 1e-10);
    CHECK(first.etaDot == 0.0);
    CHECK(first.ref == 0.08);
    CHECK(first.xiDesign == PandaExample::xi0().norm());
    CHECK(first.betaDesign == 0.0);

    // |xi| follows the designed response |xi0| (1 + 4t) e^{-4t} to 2.43e-5 [paper 2.4e-5], and beta stays at
    // 1.5e-5: both are sampling effects of order dt (see the dt test below). The reference implementation gives
    // 2.43388e-5 and 1.52440e-5, both at t = 0.409 s; a wrong term in A, Adot qdot or mu changes them.
    CAPTURE(sci(run.maxXiDeviation), run.maxXiDeviationTime, sci(run.maxBeta), run.maxBetaTime);
    CHECK(std::abs(run.maxXiDeviation - 2.43388e-5) < 0.01 * 2.43388e-5);
    CHECK_THAT(run.maxXiDeviationTime, WithinAbs(0.409, 0.0015));
    CHECK(std::abs(run.maxBeta - 1.52440e-5) < 0.01 * 1.52440e-5);
    CHECK_THAT(run.maxBetaTime, WithinAbs(0.409, 0.0015));
    CHECK(run.maxXiNorm == run.startXiNorm);  // |xi| only decreases

    // The arm: sigma_min(A) is smallest at the start, 0.0474; the largest joint velocity is 26 % of its limit
    // (joint 5 at t = 0.242 s), the largest commanded acceleration 49 % (joint 5, the first command: no ramp here)
    CAPTURE(run.minSigma, run.maxVelocityFraction, run.maxVelocityTime, run.maxAccelerationFraction);
    CHECK_THAT(run.minSigma, WithinAbs(0.04744, 5e-5));
    CHECK(run.minSigmaTime == 0.0);
    CHECK_THAT(run.maxVelocityFraction, WithinAbs(0.2622, 5e-4));
    CHECK(run.maxVelocityJoint == 4);
    CHECK_THAT(run.maxVelocityTime, WithinAbs(0.242, 0.0015));
    CHECK_THAT(run.maxAccelerationFraction, WithinAbs(0.4938, 5e-4));
    CHECK(run.maxAccelerationJoint == 4);
    CHECK(run.maxAccelerationTime == 0.0);
    CHECK(run.maxIterations == 1);  // the prediction is good enough for one Newton update per cycle
    CHECK(run.minJointMargin > 0.36);  // the smallest margin, 0.3592 rad, comes later (t = 3.54 s)

    // After 1 s: eta_dot near nu_d, |xi| at 9 % of |xi0| as designed (5 e^{-4} = 0.0916)
    CHECK_THAT(run.log.back().etaDot, WithinAbs(0.08, 2e-3));
    CHECK_THAT(run.finalXiNorm / run.startXiNorm, WithinAbs(5.0 * std::exp(-4.0), 1e-3));
}

TEST_CASE("simulate (Panda): the paper's run, 20 s (Closed loops: 2.43e-5, 1.5e-5, 7.1e-7, 0.0800025)",
          SLOW_IN_DEBUG)
{
    const PandaRun& run = paperRun();
    REQUIRE(run.ok);
    CHECK(run.T == 20.0);
    checkLog(run, 20001);
    CHECK(run.log.back().t == 20.0);

    CAPTURE(sci(run.maxXiDeviation), run.maxXiDeviationTime, sci(run.maxBeta), sci(run.finalXiNorm), run.finalEtaDot,
            sci(run.finalBeta));
    CHECK(std::abs(run.maxXiDeviation - 2.43388e-5) < 0.01 * 2.43388e-5);
    CHECK_THAT(run.maxXiDeviationTime, WithinAbs(0.409, 0.0015));
    CHECK(std::abs(run.maxBeta - 1.52440e-5) < 0.01 * 1.52440e-5);
    CHECK_THAT(run.finalXiNorm, WithinAbs(7.086e-7, 0.01e-7));
    CHECK_THAT(run.finalEtaDot, WithinAbs(0.0800025, 1e-7));
    CHECK(std::abs(run.finalBeta) < 2e-6);  // reference 1.2e-6

    // The arm along the run: smallest joint margin 0.3592 rad, joint 4 at its lower limit (the lift's own margin:
    // it is set by the path, not by the transient); sigma_min(A) 0.0474 at the start; 26 % of the velocity limit
    CAPTURE(run.minJointMargin, run.minMarginJoint, run.minMarginTime, run.minSigma, run.maxVelocityFraction);
    CHECK_THAT(run.minJointMargin, WithinAbs(0.35917, 5e-5));
    CHECK(run.minMarginJoint == 3);
    CHECK_THAT(run.minMarginTime, WithinAbs(3.537, 0.0015));
    CHECK_THAT(run.minSigma, WithinAbs(0.04744, 5e-5));
    CHECK(run.minSigmaTime == 0.0);
    CHECK_THAT(run.maxVelocityFraction, WithinAbs(0.2622, 5e-4));
    CHECK(run.maxVelocityJoint == 4);
    CHECK_THAT(run.maxAccelerationFraction, WithinAbs(0.4938, 5e-4));
    CHECK(run.maxIterations == 1);

    // Two loops and a bit: eta wraps around, so the run crosses the seam eta = L -> 0 of the path twice
    int wraps = 0;
    for (std::size_t k = 1; k < run.log.size(); ++k)
        wraps += run.log[k].eta < run.log[k - 1].eta ? 1 : 0;
    CHECK(wraps == 2);

    // Its first second is the 1 s run that stands in for it in Debug, row for row
    const PandaRun& first = paperFirstSecond();
    REQUIRE(run.log.size() > first.log.size());
    CHECK(rowDifference(run, first) == 0.0);
}

TEST_CASE("simulate (Panda): the self-motion run, first half second (Closed loops: 1.51e-4 at t = 0.41 s)")
{
    // Starts on the path with beta0 = 0.2, moving at nu_d: beta decays as designed, and the tool stays on the path.
    // In Debug it stands in for the one-loop run, which is hidden there; it holds that run's largest deviation.
    const PandaRun& run = selfMotionFirstHalfSecond();
    REQUIRE(run.ok);
    CHECK(run.kind == RunKind::SelfMotion);
    checkLog(run, 501);

    // Start: the pose gamma(0.1), beta = 0.2, and A(q0) qd0 = (nu_d, 0, ..., 0): eta_dot = 0.08, xi_dot = 0 and
    // beta_dot = 0, so the designed response beta0 (1 + 4t) e^{-4t} applies from t = 0
    const PandaExample& ex = example();
    const PandaLogRow& first = run.log.front();
    CHECK(first.q == ex.selfMotionStart());
    CHECK(run.startAttitudeDeg < 1e-8);
    CHECK(run.startPosition < 1e-10);
    CHECK(run.startXiNorm < 1e-10);
    CHECK_THAT(run.startBeta, WithinAbs(0.2, 1e-10));
    CHECK_THAT(first.etaDot, WithinAbs(0.08, 1e-12));
    CHECK(std::abs(first.eTan) < 1e-12);
    CHECK(first.xiDesign == 0.0);
    CHECK(first.betaDesign == 0.2);
    Controller<PandaArm> c(ex.tube(), PandaExample::gains(), PandaExample::reference(), &ex.redundancy());
    Controller<PandaArm>::OutputState o;
    REQUIRE(c.outputs(first.q, 0.1, o));
    const Vec7 yDot = o.A * first.qd;
    const double errRate = maxDiff(yDot, Vec7(0.08 * Vec7::Unit(0)));
    CAPTURE(sci(errRate), first.qd.transpose());
    CHECK(errRate < 1e-12);
    CHECK(first.qd.norm() > 0.1);  // the arm moves: the self-motion is not at rest

    // The peak of |beta - beta0 (1 + 4t) e^{-4t}|: 1.5103e-4 at t = 0.405 s in the reference (accept < 2e-4)
    CAPTURE(sci(run.maxBetaDeviation), run.maxBetaDeviationTime, sci(run.maxXiNorm));
    CHECK(run.maxBetaDeviation < 2e-4);
    CHECK(std::abs(run.maxBetaDeviation - 1.51035e-4) < 0.01 * 1.51035e-4);
    CHECK_THAT(run.maxBetaDeviationTime, WithinAbs(0.405, 0.0015));
    CHECK(run.maxXiNorm < 1.3e-6);  // the self-motion does not move the tool off the path
    CHECK(run.maxIterations == 1);
}

TEST_CASE("simulate (Panda): the self-motion run, one loop (Closed loops: 1.2e-6, 1.51e-4, 9.9e-6, 1.3e-5)",
          SLOW_IN_DEBUG)
{
    const PandaRun& run = selfMotionRun();
    REQUIRE(run.ok);
    CHECK_THAT(run.T, WithinAbs(0.77008765 / 0.08, 1e-6));
    checkLog(run, 9627);  // t_k = k dt < L / nu_d: one loop
    CHECK(run.log.back().t == 9626 * 1e-3);

    // The tool stays on the path (max |xi| 1.2e-6) while beta decays as designed (to 1.51e-4, accept < 2e-4)
    CAPTURE(sci(run.maxXiNorm), sci(run.maxBetaDeviation), run.maxBetaDeviationTime);
    CHECK(run.maxXiNorm < 1.3e-6);
    CHECK(std::abs(run.maxXiNorm - 1.24497e-6) < 0.01 * 1.24497e-6);
    CHECK(run.maxBetaDeviation < 2e-4);
    CHECK(std::abs(run.maxBetaDeviation - 1.51035e-4) < 0.01 * 1.51035e-4);
    CHECK_THAT(run.maxBetaDeviationTime, WithinAbs(0.405, 0.0015));
    CHECK(run.maxIterations == 1);

    // After one loop the arm is back at the start of the path and on the closed lift: max_i |q_i - q*_i(eta)|
    // 9.9e-6 and 2-norm 1.3e-5 (accept < 2e-5), from 0.2 rad at the start
    const PandaLogRow& last = run.log.back();
    CAPTURE(sci(run.finalLiftGap), sci(run.finalLiftGapNorm), last.eta, sci(last.beta));
    CHECK(std::abs(last.eta - 0.1) < 1e-3);
    CHECK(run.finalLiftGapNorm < 2e-5);
    CHECK(std::abs(run.finalLiftGap - 9.932e-6) < 0.01 * 9.932e-6);
    CHECK(std::abs(run.finalLiftGapNorm - 1.3024e-5) < 0.01 * 1.3024e-5);
    Vec7 qStar = Vec7::Zero();
    example().lift().eval(last.eta, qStar);
    CHECK(run.finalLiftGap == (last.q - qStar).cwiseAbs().maxCoeff());
    CHECK(std::abs(last.beta) < 5e-6);  // reference -2.6e-6
    CHECK(run.minJointMargin > 0.35);

    // Its first half second is the run that stands in for it in Debug, row for row
    CHECK(rowDifference(run, selfMotionFirstHalfSecond()) == 0.0);
}

TEST_CASE("simulate (Panda): the deviations are sample-and-hold effects: halving dt halves them", SLOW_IN_DEBUG)
{
    // Design doc, Three kinds of test: if halving dt does not halve the deviation from the designed response, the
    // controller has a bug that the sampling does not explain. Both peaks are in the first second (t = 0.41 s).
    const PandaRun& coarse = paperFirstSecond();
    const PandaRun fine = simulate(example(), RunKind::Paper, 1.0, 0.5e-3);
    REQUIRE(fine.ok);
    checkLog(fine, 2001);
    const double ratioXi = coarse.maxXiDeviation / fine.maxXiDeviation;
    const double ratioBeta = coarse.maxBeta / fine.maxBeta;
    CAPTURE(sci(fine.maxXiDeviation), sci(fine.maxBeta), ratioXi, ratioBeta);
    CHECK(ratioXi > 1.9);
    CHECK(ratioXi < 2.1);
    CHECK(ratioBeta > 1.9);
    CHECK(ratioBeta < 2.1);
    CHECK(fine.log[0].qdd == coarse.log[0].qdd);  // the first command does not depend on dt
}

TEST_CASE("simulate (Panda): argument checks; T = 0 gives the first command only")
{
    const PandaExample& ex = example();
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(simulate(ex, RunKind::Paper, -1.0), std::invalid_argument);
    CHECK_THROWS_AS(simulate(ex, RunKind::Paper, nan), std::invalid_argument);
    CHECK_THROWS_AS(simulate(ex, RunKind::Paper, inf), std::invalid_argument);
    CHECK_THROWS_AS(simulate(ex, RunKind::Paper, 1.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(simulate(ex, RunKind::SelfMotion, 1.0, -1e-3), std::invalid_argument);
    CHECK_THROWS_AS(simulate(ex, RunKind::SelfMotion, 1.0, nan), std::invalid_argument);

    // Finite but absurd lengths: T / dt above PandaExample::maxSteps (1e7). Without that check, 1e18 steps would
    // make reserve throw std::length_error, 1e303 steps would overflow the cast to long (undefined behaviour), and
    // T / dt = 1 / 1e-320 overflows to infinity. (Values just above the bound are not tried here: if the check were
    // missing, such a run would take hours instead of failing.)
    STATIC_REQUIRE(PandaExample::maxSteps == 10'000'000);
    CHECK_THROWS_AS(simulate(ex, RunKind::Paper, 1e15), std::invalid_argument);
    CHECK_THROWS_AS(simulate(ex, RunKind::Paper, 1e300), std::invalid_argument);
    CHECK_THROWS_AS(simulate(ex, RunKind::SelfMotion, 1.0, 1e-320), std::invalid_argument);

    for (RunKind kind : {RunKind::Paper, RunKind::SelfMotion}) {
        const PandaRun one = simulate(ex, kind, 0.0);
        REQUIRE(one.ok);
        REQUIRE(one.log.size() == 1);
        const PandaRun& longer = kind == RunKind::Paper ? paperFirstSecond() : selfMotionFirstHalfSecond();
        CHECK(rowDifference(one, longer) == 0.0);  // the same first row (and command) as the longer run
        CHECK(one.finalXiNorm == one.startXiNorm);
    }

    // The paper's first command at rest: 49 % of joint 5's acceleration limit (the hardware ramp brings it down)
    const PandaRun one = simulate(ex, RunKind::Paper, 0.0);
    CHECK_THAT(one.maxAccelerationFraction, WithinAbs(0.4938, 5e-4));
}

TEST_CASE("simulate (Panda): a failed step stops the run; the summary covers the rows before it")
{
    // With dt = 0.3 s the sampled loop is unstable (the poles at -4 need 4 dt well below 1): the arm swings ever
    // further until the step at t = 7 dt = 2.1 s fails. There |q - q*(eta)| is about 4e10 rad, so the X' term of
    // d beta (eq:dbeta) makes the last row of A of order 1e8, and sigma_min(A) = 1e-10 is below the step's 1e-9.
    const PandaRun run = simulate(example(), RunKind::Paper, 20.0, 0.3);
    REQUIRE_FALSE(run.ok);
    CAPTURE(run.failureTime, run.log.size(), run.minJointMargin);
    CHECK_THAT(run.failureTime, WithinAbs(2.1, 1e-12));
    CHECK(run.T == 20.0);
    CHECK(run.log.capacity() == 67);  // reserved for t = 0, 0.3, ..., 19.8 before the loop, filled up to 1.8 s
    checkLog(run, 7);                 // t = 0, ..., 1.8 s: the extremes and the update over these rows
    CHECK_THAT(run.log.back().t, WithinAbs(1.8, 1e-12));

    // ok covers the controller only: before its step failed, the arm left its joint range by thousands of radians
    // (a model without joint limits), which only minJointMargin < 0 reports.
    CHECK(run.minJointMargin < -1000.0);
    CHECK(run.log.front().q == example().paperStart());  // the same start as the 1 ms run
}
