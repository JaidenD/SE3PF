// Tests of the controller (Section V, Algorithm 1): generic parts, the rigid body (n = 6), and the Panda (n = 7),
// first with the joint coupling and a test-only redundancy output, then with the lift's beta of Lemma 3 in the
// paper's example.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "se3pf/controller.hpp"
#include "se3pf/experiments.hpp"
#include "se3pf/lift.hpp"
#include "se3pf/panda.hpp"
#include "se3pf/panda_example.hpp"
#include "se3pf/rigid_body.hpp"
#include "se3pf/rigid_body_example.hpp"
#include "se3pf/safety.hpp"
#include "se3pf/se3.hpp"
#include "panda_fixture.hpp"
#include "test_helpers.hpp"
#include "tube_fixture.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace Eigen;
using namespace se3pf;
using Catch::Matchers::WithinAbs;
using testing::arcDiff;
using testing::maxDiff;
using testing::pandaExample;
using testing::pandaTube;
using testing::randn;
using testing::rigidTube;
using testing::sci;
using testing::uniform;

// ---------------------------------------------------------------------------------------------------------------
// Counting global operator new, for the no-allocation tests (design doc, Module: controller, Tests). It replaces
// the global operator new of the whole test program, so allocations of every test and every thread are counted; only
// the difference across a window of control steps matters, and that window runs on one thread. Eigen allocates with
// malloc instead, which EIGEN_RUNTIME_NO_MALLOC checks (Debug only).

namespace {
    // Atomic: a test may allocate on a thread of its own (started with std::async, say), and unsynchronized
    // increments of a plain long from two threads would be a data race. Relaxed order is enough for a counter.
    std::atomic<long> gNewCalls{0};
}

void* operator new(std::size_t size)
{
    gNewCalls.fetch_add(1, std::memory_order_relaxed);
    if (void* p = std::malloc(size == 0 ? 1 : size))
        return p;
    throw std::bad_alloc();
}

void* operator new(std::size_t size, std::align_val_t alignment)
{
    gNewCalls.fetch_add(1, std::memory_order_relaxed);
    const std::size_t a = static_cast<std::size_t>(alignment);
    if (void* p = std::aligned_alloc(a, (size + a - 1) / a * a))  // the size must be a multiple of the alignment
        return p;
    throw std::bad_alloc();
}

void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete(void* p, std::align_val_t) noexcept { std::free(p); }
void operator delete(void* p, std::size_t, std::align_val_t) noexcept { std::free(p); }

namespace {

    // -----------------------------------------------------------------------------------------------------------
    // Fixtures

    using RigidController = Controller<RigidBody>;
    using PandaController = Controller<PandaArm>;

    // The tubes are those of tube_fixture.hpp: rigidTube() (tumbling circle, r = 1 m, W = I) and pandaTube()
    // (coning circle, l = 0.1 m).

    Vec5 rigidXi0() { return RigidBodyExample::xi0(); }

    // The rigid body's start pose Psi(0.1, xi0): 25.45 deg and 0.4389 m off the path, so (eta, xi) and A are far
    // from their values on the path.
    Pose rigidStart()
    {
        return rigidTube().Psi(0.1, rigidXi0());
    }

    // Gains that differ from each other, so that a gain used in the wrong place shows.
    Gains oddGains()
    {
        Gains k;
        k.Kp = 11.0;
        k.Kd = 7.0;
        k.kt = 3.0;
        k.kpEta = 5.0;
        k.kdEta = 2.0;
        k.KpBeta = 13.0;
        k.KdBeta = 3.0;
        return k;
    }

    // References with every term of eq:tan-vel and eq:tan-pos active.
    // nu_d(t, eta) = 0.4 + 0.1 sin t + 0.05 cos(2 pi eta / L), with its partial derivatives
    double nuTest(double t, double eta, double L)
    {
        return 0.4 + 0.1 * std::sin(t) + 0.05 * std::cos(2.0 * EIGEN_PI * eta / L);
    }
    double dnuDtTest(double t) { return 0.1 * std::cos(t); }
    double dnuDetaTest(double eta, double L)
    {
        return -0.05 * (2.0 * EIGEN_PI / L) * std::sin(2.0 * EIGEN_PI * eta / L);
    }

    TangentialReference velocityTest(double L)
    {
        return TangentialReference::fromVelocity([L](double t, double eta, double& nu, double& dnuDt, double& dnuDeta) {
            nu = nuTest(t, eta, L);
            dnuDt = dnuDtTest(t);
            dnuDeta = dnuDetaTest(eta, L);
        });
    }

    // eta_d(t) = 1 + 0.3 t + 0.2 sin t, real-valued: at t = 40 it is one loop past L on the rigid tube.
    double etaDTest(double t) { return 1.0 + 0.3 * t + 0.2 * std::sin(t); }
    double detaDTest(double t) { return 0.3 + 0.2 * std::cos(t); }
    double ddetaDTest(double t) { return -0.2 * std::sin(t); }

    TangentialReference positionTest()
    {
        return TangentialReference::fromPosition([](double t, double& e, double& de, double& dde) {
            e = etaDTest(t);
            de = detaDTest(t);
            dde = ddetaDTest(t);
        });
    }

    // eta_d = constant
    TangentialReference positionConstant(double etaD)
    {
        return TangentialReference::fromPosition([etaD](double, double& e, double& de, double& dde) {
            e = etaD;
            de = 0.0;
            dde = 0.0;
        });
    }

    // A curve through g with velocity zeta that differs from RigidBody::retract at second order: the screw motion
    // g exp(h zeta^) = g (exp(h omega), J_l(h omega) h v), with J_l(w) = J_r(w)^T. A central difference along any such
    // curve gives the same derivative, so checks that use it do not just repeat the controller's computation.
    Pose screw(const Pose& g, const Vec6& zeta, double h)
    {
        const Vector3d w = h * zeta.head<3>();
        const Vector3d v = h * zeta.tail<3>();
        return g * makePose(expSO3(w), Jr(w).transpose() * v);
    }

    Pose curve(const Pose& g, const Vec6& zeta, double h) { return screw(g, zeta, h); }
    Vec7 curve(const Vec7& q, const Vec7& dq, double h) { return q + h * dq; }

    // ybar = (eta, xi[, beta]) at q, with eta relative to etaRef on the circle (so it can be differenced)
    template <class System>
    Vec<System::kDof> ybar(const Controller<System>& c, const typename System::Config& q, double etaRef)
    {
        typename Controller<System>::OutputState o;
        Vec<System::kDof> y = Vec<System::kDof>::Constant(std::numeric_limits<double>::quiet_NaN());
        if (!c.outputs(q, etaRef, o))
            return y;
        y[0] = arcDiff(o.eta, etaRef, c.tube().length());
        y.template segment<5>(1) = o.xi;
        if constexpr (System::kDof == 7)
            y[6] = o.beta;
        return y;
    }

    // Adot qdot by a fourth-order central difference along curve(), not along System::retract:
    // (-A(2h) + 8 A(h) - 8 A(-h) + A(-2h)) qdot / 12h, independent of the controller's own second-order difference.
    // The frame spline N(s) is only C^1: N'' jumps at the grid nodes, and so does the derivative of A along the
    // motion. A stencil that spans a node mixes the two one-sided values (an error up to about 1e-6), so h is small
    // (the stencil moves eta by about 4h eta_dot) and the tests put eta well inside a grid interval.
    template <class System>
    Vec<System::kDof> adotQdot(const Controller<System>& c, const typename System::Config& q,
                               const Vec<System::kDof>& dq, double eta)
    {
        constexpr int n = System::kDof;
        const double h = 1e-4;
        Mat<n, n> sum = Mat<n, n>::Zero();
        const double steps[4] = {-2.0, -1.0, 1.0, 2.0}, weights[4] = {1.0, -8.0, 8.0, -1.0};
        for (int i = 0; i < 4; ++i) {
            typename Controller<System>::OutputState o;
            REQUIRE(c.outputs(curve(q, dq, steps[i] * h), eta, o));
            sum += weights[i] * o.A;
        }
        return (sum * dq) / (12.0 * h);
    }

    // A is the Jacobian of ybar (eq:A): (ybar(retract(q, w, h)) - ybar(retract(q, w, -h))) / 2h = A w for random w.
    template <class System>
    void checkJacobian(const Controller<System>& c, const typename System::Config& q, double guess, double scale)
    {
        constexpr int n = System::kDof;
        typename Controller<System>::OutputState o;
        REQUIRE(c.outputs(q, guess, o));
        for (int i = 0; i < 6; ++i) {
            const Vec<n> w = randn<Vec<n>>(scale);
            const double h = 1e-6;
            const Vec<n> fd = (ybar(c, System::retract(q, w, h), o.eta) - ybar(c, System::retract(q, w, -h), o.eta))
                            / (2.0 * h);
            const Vec<n> Aw = o.A * w;
            const double err = maxDiff(fd, Aw);
            CAPTURE(i, sci(err), fd.transpose(), Aw.transpose());
            CHECK(err < 1e-6);
        }
    }

    // mu of Section V computed again from what the step reported: eq:tan-vel or eq:tan-pos with the given
    // reference values, eq:transverse and (n = 7) eq:beta-law.
    template <class System>
    Vec<System::kDof> expectedMu(const Gains& k, const typename Controller<System>::Status& st, bool velocity,
                                 double ref, double dref, double ddref, double L)
    {
        Vec<System::kDof> mu = Vec<System::kDof>::Zero();
        if (velocity)  // ref = nu_d, dref = d nu_d / dt, ddref = d nu_d / d eta
            mu[0] = dref + ddref * st.etaDot - k.kt * (st.etaDot - ref);
        else  // ref = eta_d, dref = eta_d', ddref = eta_d''
            mu[0] = ddref - k.kdEta * (st.etaDot - dref) - k.kpEta * arcDiff(st.eta, ref, L);
        for (int i = 0; i < 5; ++i)
            mu[i + 1] = -k.Kp * st.xi[i] - k.Kd * st.xiDot[i];
        if constexpr (System::kDof == 7)
            mu[6] = -k.KpBeta * st.beta - k.KdBeta * st.betaDot;
        return mu;
    }

    // -----------------------------------------------------------------------------------------------------------
    // Panda fixtures (n = 7)

    // The joint coupling beta = q3 - q1 (CouplingRedundancy, design doc, Module: lift and redundancy output). Linear
    // in q, so d beta = e3 - e1 everywhere, and {beta = value} is the hyperplane (e3 - e1)^T q = value (tested in
    // test_lift.cpp).
    const CouplingRedundancy& coupling()
    {
        static const CouplingRedundancy c;
        return c;
    }

    // A second test-only output that depends on eta, so that the d eta term of eq:dbeta is active:
    // beta = q3 - q1 + a sin(k eta), d beta = e3 - e1 + a k cos(k eta) d eta, with a = 0.2 and k = 8 (one wave per
    // 0.79 of arc length, about one loop of the Panda path). The coupling ignores eta and d eta, so with it alone a
    // controller that passed the wrong row of J_Psi^-1 J_b as d eta, or evaluated beta at the guess instead of at
    // the projected eta, would still pass.
    class WavyCouplingOutput final : public RedundancyOutput {
        public:
            static constexpr double kAmplitude = 0.2, kWavenumber = 8.0;  // a and k

            // The part a sin(k eta) and its derivative a k cos(k eta), for the tests' own expected values
            static double wave(double eta) { return kAmplitude * std::sin(kWavenumber * eta); }
            static double waveSlope(double eta) { return kAmplitude * kWavenumber * std::cos(kWavenumber * eta); }

            void evaluate(const Vec7& q, double eta, const Vec7& deta, double& beta, Vec7& dbeta) const override
            {
                beta = q[2] - q[0] + wave(eta);
                dbeta = Vec7::Unit(2) - Vec7::Unit(0) + waveSlope(eta) * deta;
            }

            // At a fixed eta, {beta = value} is the hyperplane q3 - q1 = value - a sin(k eta).
            void levelSet(double eta, double value, Vec7& p, double& b) const override
            {
                p = Vec7::Unit(2) - Vec7::Unit(0);
                b = value - wave(eta);
            }
    };

    const WavyCouplingOutput& wavyCoupling()
    {
        static const WavyCouplingOutput c;
        return c;
    }

    // A configuration with flange pose g and q3 - q1 = beta (design doc, Start configurations, with the coupling
    // in place of the lift): Newton's method on [J_b; p^T] dq = (e, beta - p^T q), e = Log_I(Pi(q)^-1 g), from the
    // rest posture, with steps capped at 0.2 rad.
    Vec7 pandaConfiguration(const Pose& g, double beta)
    {
        Vec7 p = Vec7::Zero();
        double b = 0.0;
        coupling().levelSet(0.0, beta, p, b);
        Vec7 q = PandaArm::restPosture();
        for (int k = 0; k < 100; ++k) {
            Mat67 Jb = Mat67::Zero();
            const Pose gq = PandaArm::output(q, &Jb);
            const Pose gErr = inverse(gq) * g;
            Vec6 e = Vec6::Zero();
            e << logSO3(Matrix3d(gErr.block<3,3>(0,0))), gErr.block<3,1>(0,3);
            const double c = b - p.dot(q);
            if (e.cwiseAbs().maxCoeff() < 1e-13 && std::abs(c) < 1e-13)
                return q;
            Mat<7, 7> M = Mat<7, 7>::Zero();
            M.topRows<6>() = Jb;
            M.row(6) = p.transpose();
            Vec7 rhs = Vec7::Zero();
            rhs << e, c;
            Vec7 dq = M.partialPivLu().solve(rhs);
            if (dq.norm() > 0.2)
                dq *= 0.2 / dq.norm();
            q += dq;
        }
        FAIL("pandaConfiguration: Newton's method did not converge");
        return q;
    }

    Vec5 pandaXi0() { return PandaExample::xi0(); }

    // The Panda's start pose Psi(0.1, xi0) with beta = q3 - q1 = 0.03: off the path and off the coupled lift.
    const Vec7& pandaStart()
    {
        static const Vec7 q = pandaConfiguration(pandaTube().Psi(0.1, pandaXi0()), 0.03);
        return q;
    }

    // A test-only System: a rigid body that can hardly translate along its z axis (J_b = I with its last column
    // scaled by 1e-12). A = J_Psi^-1 J_b is then singular to working precision but not exactly, so an LU solve
    // returns huge finite values instead of failing; the step must refuse it by its sigma_min check.
    struct StuckBody {
        static constexpr int kDof = 6;
        using Config = Pose;
        static Pose output(const Pose& g, Mat6* Jb)
        {
            *Jb = Mat6::Identity();
            (*Jb)(5, 5) = 1e-12;
            return g;
        }
        static Pose retract(const Pose& g, const Vec6& zeta, double h) { return RigidBody::retract(g, zeta, h); }
    };

}  // namespace

// ===============================================================================================================
// Generic parts: gains, references, the ramp, construction

TEST_CASE("Gains: the paper's defaults; fromLambda puts the transverse and redundancy poles at -lambda")
{
    const Gains k;
    CHECK(k.Kp == 16.0);
    CHECK(k.Kd == 8.0);
    CHECK(k.kt == 5.0);
    CHECK(k.kpEta == 9.0);
    CHECK(k.kdEta == 6.0);
    CHECK(k.KpBeta == 16.0);
    CHECK(k.KdBeta == 8.0);

    const Gains g = Gains::fromLambda(3.0);  // s^2 + Kd s + Kp = (s + 3)^2
    CHECK(g.Kp == 9.0);
    CHECK(g.Kd == 6.0);
    CHECK(g.KpBeta == 9.0);
    CHECK(g.KdBeta == 6.0);
    CHECK(g.kt == 5.0);
}

TEST_CASE("quinticRamp: 0 and 1 at the ends with flat first and second derivatives, clamped outside")
{
    CHECK(quinticRamp(0.0) == 0.0);
    CHECK(quinticRamp(1.0) == 1.0);
    CHECK_THAT(quinticRamp(0.5), WithinAbs(0.5, 1e-15));  // symmetric: r(x) + r(1 - x) = 1
    CHECK_THAT(quinticRamp(0.3) + quinticRamp(0.7), WithinAbs(1.0, 1e-15));
    CHECK(quinticRamp(-0.2) == 0.0);
    CHECK(quinticRamp(1.7) == 1.0);
    CHECK(quinticRamp(std::numeric_limits<double>::quiet_NaN()) == 0.0);

    // r(h) = 10 h^3 + O(h^4) and 1 - r(1 - h) = 10 h^3 + O(h^4): value, slope and curvature vanish at both ends
    const double h = 1e-3;
    CHECK_THAT(quinticRamp(h) / (h * h * h), WithinAbs(10.0, 0.02));
    CHECK_THAT((1.0 - quinticRamp(1.0 - h)) / (h * h * h), WithinAbs(10.0, 0.02));
    // and it increases
    for (int i = 0; i < 100; ++i)
        CHECK(quinticRamp(0.01 * (i + 1)) > quinticRamp(0.01 * i));
}

TEST_CASE("wrapSymmetric: the representative in (-L/2, L/2]")
{
    const double L = 2.0;
    CHECK_THAT(wrapSymmetric(0.3, L), WithinAbs(0.3, 1e-15));
    CHECK_THAT(wrapSymmetric(1.7, L), WithinAbs(-0.3, 1e-15));
    CHECK_THAT(wrapSymmetric(-1.7, L), WithinAbs(0.3, 1e-15));
    CHECK_THAT(wrapSymmetric(6.3, L), WithinAbs(0.3, 1e-14));
    CHECK_THAT(wrapSymmetric(-5.9, L), WithinAbs(0.1, 1e-14));
    CHECK(wrapSymmetric(1.0, L) == 1.0);   // L/2 is in the range
    CHECK(wrapSymmetric(-1.0, L) == 1.0);  // -L/2 is not
    for (int i = 0; i < 100; ++i) {
        const double x = uniform(-20.0, 20.0), w = wrapSymmetric(x, L);
        CHECK(w > -0.5 * L);
        CHECK(w <= 0.5 * L);
        CHECK_THAT(std::remainder(x - w, L), WithinAbs(0.0, 1e-13));  // differs from x by a multiple of L
    }
}

TEST_CASE("stepCount: the instants k dt <= T of a simulated run, and its argument checks")
{
    CHECK(stepCount(20.0, 1e-3, 10'000'000) == 20000);  // t = 20 s included despite rounding: 20001 rows
    CHECK(stepCount(20.0, 0.002, 10'000'000) == 10000);
    CHECK(stepCount(20.0, 0.3, 10'000'000) == 66);       // 66 * 0.3 = 19.8 <= 20 < 67 * 0.3
    CHECK(stepCount(0.77008765 / 0.08, 1e-3, 10'000'000) == 9626);  // one loop of the Panda path: 9627 rows
    CHECK(stepCount(0.0, 1e-3, 10'000'000) == 0);
    CHECK(stepCount(1.0, 1e-3, 1000) == 1000);           // exactly at the limit
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(stepCount(1.0, 1e-3, 999), std::invalid_argument);  // one step too many
    CHECK_THROWS_AS(stepCount(-1.0, 1e-3, 1000), std::invalid_argument);
    CHECK_THROWS_AS(stepCount(1.0, 0.0, 1000), std::invalid_argument);
    CHECK_THROWS_AS(stepCount(1.0, -1e-3, 1000), std::invalid_argument);
    CHECK_THROWS_AS(stepCount(nan, 1e-3, 1000), std::invalid_argument);
    CHECK_THROWS_AS(stepCount(inf, 1e-3, 1000), std::invalid_argument);
    CHECK_THROWS_AS(stepCount(1.0, nan, 1000), std::invalid_argument);
    CHECK_THROWS_AS(stepCount(1.0, 1e-320, 10'000'000), std::invalid_argument);  // T / dt overflows to inf
}

TEST_CASE("TangentialReference: factories set the mode and the function; valid() checks the selected one")
{
    const TangentialReference v = TangentialReference::constantVelocity(0.3);
    CHECK(v.mode == TangentialReference::Mode::Velocity);
    REQUIRE(v.valid());
    double nu = 0.0, dt = 1.0, de = 1.0;
    v.velocity(2.0, 1.0, nu, dt, de);
    CHECK(nu == 0.3);
    CHECK(dt == 0.0);
    CHECK(de == 0.0);

    const TangentialReference p = positionConstant(1.5);
    CHECK(p.mode == TangentialReference::Mode::Position);
    CHECK(p.valid());

    TangentialReference empty;
    CHECK_FALSE(empty.valid());
    empty.position = p.position;  // set, but not for the selected mode (velocity)
    CHECK_FALSE(empty.valid());
}

TEST_CASE("Controller: construction checks; a temporary tube does not compile")
{
    const Tube& T = rigidTube();
    CHECK_NOTHROW(RigidController(T, Gains(), TangentialReference::constantVelocity(0.5)));
    CHECK_THROWS_AS(RigidController(T, Gains(), TangentialReference()), std::invalid_argument);  // no function
    CHECK_THROWS_AS(RigidController(T, Gains(), TangentialReference::constantVelocity(0.5), &coupling()),
                    std::invalid_argument);  // n = 6 has no redundancy output

    RigidController c(T, Gains(), TangentialReference::constantVelocity(0.5));
    CHECK_THROWS_AS(c.setReference(TangentialReference()), std::invalid_argument);
    CHECK_THROWS_AS(c.setRampTime(-1.0), std::invalid_argument);
    CHECK_THROWS_AS(c.setRampTime(std::numeric_limits<double>::infinity()), std::invalid_argument);
    CHECK_THROWS_AS(c.setRampTime(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    CHECK(c.rampTime() == 0.0);  // off by default
    c.setRampTime(1.5);
    CHECK(c.rampTime() == 1.5);
    c.setReference(positionConstant(2.0));
    CHECK(c.reference().mode == TangentialReference::Mode::Position);

    STATIC_REQUIRE(std::is_constructible_v<RigidController, const Tube&, const Gains&, TangentialReference>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<RigidController, Tube&&, const Gains&, TangentialReference>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<PandaController, Tube&&, const Gains&, TangentialReference,
                                                 const RedundancyOutput*>);
}

TEST_CASE("Controller: reset wraps the prediction into [0, L); the constructor sets it to 0")
{
    const Tube& T = rigidTube();
    const double L = T.length();
    RigidController c(T, Gains(), TangentialReference::constantVelocity(0.5));
    CHECK(c.prediction() == 0.0);
    c.reset(0.1);
    CHECK_THAT(c.prediction(), WithinAbs(0.1, 1e-15));
    c.reset(L + 0.3);
    CHECK_THAT(c.prediction(), WithinAbs(0.3, 1e-14));
    c.reset(-0.2);
    CHECK_THAT(c.prediction(), WithinAbs(L - 0.2, 1e-14));
    c.reset(-3.0 * L + 1.0);
    CHECK_THAT(c.prediction(), WithinAbs(1.0, 1e-13));
    c.reset(L);
    CHECK(c.prediction() == 0.0);
}

// ===============================================================================================================
// Rigid body (n = 6)

TEST_CASE("Controller<RigidBody>: A is the Jacobian of ybar = (eta, xi) off the path (eq:A)")
{
    const Tube& T = rigidTube();
    const RigidController c(T, Gains(), TangentialReference::constantVelocity(0.5));

    // For the rigid body J_b = I, so A = J_Psi^-1
    RigidController::OutputState o;
    REQUIRE(c.outputs(rigidStart(), 0.1, o));
    CHECK_THAT(o.eta, WithinAbs(0.1, 1e-9));
    CHECK(maxDiff(o.xi, rigidXi0()) < 1e-9);
    CHECK(o.iterations >= 1);
    const double errInverse = maxDiff(Mat6(o.A * T.JPsi(o.eta, o.xi)), Mat6::Identity());
    CAPTURE(sci(errInverse));
    CHECK(errInverse < 1e-12);

    // Central differences of (eta, xi) along random body twists, at the start pose and at random poses off the path
    checkJacobian(c, rigidStart(), 0.1, 1.0);
    for (int i = 0; i < 4; ++i) {
        const double s = uniform(0.0, T.length());
        const Vec5 xi = randn<Vec5>(0.25);
        CAPTURE(s, xi.norm());
        checkJacobian(c, T.Psi(s, xi), s, 1.0);
    }
}

TEST_CASE("Controller<RigidBody>: on the path xi = 0 and A is orthogonal (W = I), so sigma_min(A) = 1")
{
    // On the path J_Psi = [zeta_gamma | N], a W-orthonormal basis; with W = I and J_b = I, A = J_Psi^-1 is
    // orthogonal (Proposition 2: A is invertible on the path).
    const Tube& T = rigidTube();
    RigidController c(T, Gains(), TangentialReference::constantVelocity(0.5));
    for (int i = 0; i < 10; ++i) {
        const double s = uniform(0.0, T.length());
        RigidController::OutputState o;
        REQUIRE(c.outputs(T.path().pose(s), s + 0.01, o));
        const double sigma = JacobiSVD<Mat6>(o.A).singularValues().minCoeff();
        const double errOrtho = maxDiff(Mat6(o.A.transpose() * o.A), Mat6::Identity());
        CAPTURE(s, sci(o.xi.norm()), sci(errOrtho));
        CHECK(std::abs(arcDiff(o.eta, s, T.length())) < 1e-9);
        CHECK(o.xi.norm() < 1e-9);
        CHECK(errOrtho < 1e-9);
        CHECK_THAT(sigma, WithinAbs(1.0, 1e-9));

        // and the step reports the same sigma_min
        c.reset(s);
        Vec6 qdd = Vec6::Zero();
        RigidController::Status st;
        REQUIRE(c.step(T.path().pose(s), Vec6::Zero(), 0.0, 0.002, qdd, st));
        CHECK_THAT(st.sigmaMin, WithinAbs(1.0, 1e-9));
    }
}

TEST_CASE("Controller<RigidBody>: at rest Adot qdot = 0, so A qdd* = mu exactly")
{
    const Tube& T = rigidTube();
    const Gains k = oddGains();
    RigidController c(T, k, velocityTest(T.length()));
    const double t = 0.7;

    // Garbage in the output and the status first: the step must overwrite all of it.
    Vec6 qdd = Vec6::Constant(1e300);
    RigidController::Status st;
    st.ok = true;
    st.etaDot = st.beta = 1e300;
    c.reset(0.1);
    REQUIRE(c.step(rigidStart(), Vec6::Zero(), t, 0.002, qdd, st));
    CHECK(st.ok);
    CHECK(st.t == t);
    CHECK(st.etaDot == 0.0);
    CHECK(st.xiDot.isZero(0.0));
    CHECK(st.beta == 0.0);  // n = 6: no redundancy output
    CHECK(st.ramp == 1.0);

    RigidController::OutputState o;
    REQUIRE(c.outputs(rigidStart(), 0.1, o));
    const Vec6 mu = expectedMu<RigidBody>(k, st, true, nuTest(t, st.eta, T.length()), dnuDtTest(t),
                                         dnuDetaTest(st.eta, T.length()), T.length());
    const double err = maxDiff(Vec6(o.A * qdd), mu);
    CAPTURE(sci(err), mu.transpose());
    CHECK(err < 1e-12);
    CHECK(mu.cwiseAbs().maxCoeff() > 1.0);  // a real test: mu is not small here
}

TEST_CASE("Controller<RigidBody>: step returns qdd* with A qdd* + Adot qdot = mu (eq:fb, Algorithm 1)")
{
    const Tube& T = rigidTube();
    const double L = T.length();
    const Gains k = oddGains();

    for (bool velocity : {true, false}) {
        RigidController c(T, k, velocity ? velocityTest(L) : positionTest());
        for (int i = 0; i < 4; ++i) {
            // Off the path, moving: every term of A, Adot qdot and mu is active. t = 40 puts eta_d past L. The foot
            // point s is the middle of a grid interval of the frame spline (see adotQdot).
            const double ds = L / T.gridIntervals();
            const double s = (std::floor(uniform(0.0, T.gridIntervals())) + 0.5) * ds;
            const double t = i == 3 ? 40.0 : uniform(0.0, 10.0), dt = 0.002;
            const Pose g = T.Psi(s, randn<Vec5>(0.2));
            const Vec6 dq = randn<Vec6>(0.7);
            c.reset(s + 0.01);
            Vec6 qdd = Vec6::Zero();
            RigidController::Status st;
            REQUIRE(c.step(g, dq, t, dt, qdd, st));

            // Independent: A at q, ybar_dot = A qdot, Adot qdot along another curve, mu from the status
            RigidController::OutputState o;
            REQUIRE(c.outputs(g, s, o));
            const Vec6 yDot = o.A * dq;
            const Vec6 Adq = adotQdot(c, g, dq, o.eta);
            const Vec6 mu = velocity
                ? expectedMu<RigidBody>(k, st, true, nuTest(t, st.eta, L), dnuDtTest(t), dnuDetaTest(st.eta, L), L)
                : expectedMu<RigidBody>(k, st, false, etaDTest(t), detaDTest(t), ddetaDTest(t), L);

            const double errEta = std::abs(st.eta - o.eta);
            const double errYDot =
                std::max(std::abs(st.etaDot - yDot[0]), maxDiff(st.xiDot, Vec5(yDot.segment<5>(1))));
            const double errXi = maxDiff(st.xi, o.xi);
            const double errMu = maxDiff(Vec6(o.A * qdd + Adq), mu);
            CAPTURE(velocity, s, t, sci(errEta), sci(errYDot), sci(errXi), sci(errMu), Adq.transpose());
            CHECK(errEta < 1e-12);
            CHECK(errYDot < 1e-12);
            CHECK(errXi < 1e-12);
            // The controller's Adot qdot is a second-order difference with h = 1e-5, whose error grows as |qdot|^3:
            // about 1e-9 typically, up to 2.4e-8 over 100 random draws of qdot ~ N(0, 0.7^2) (the test's own
            // fourth-order difference adds less). The bound leaves room for that, and is still 1e5 times below
            // |Adot qdot| (next check), which a missing or wrong Adot qdot would leave in errMu.
            CHECK(errMu < 1e-7);
            CHECK(Adq.cwiseAbs().maxCoeff() > 0.01);  // Adot qdot matters: dropping it would fail errMu

            // Status: the reference and the tangential error
            if (velocity) {
                CHECK_THAT(st.ref, WithinAbs(nuTest(t, st.eta, L), 1e-15));
                CHECK_THAT(st.eTan, WithinAbs(st.etaDot - st.ref, 1e-15));
            } else {
                CHECK_THAT(st.ref, WithinAbs(etaDTest(t), 1e-15));
                CHECK_THAT(st.eTan, WithinAbs(arcDiff(st.eta, st.ref, L), 1e-12));
                CHECK(std::abs(st.eTan) <= 0.5 * L);
            }
            CHECK(st.iterations >= 1);

            // The prediction for the next cycle: eta + eta_dot dt, wrapped
            CHECK(std::abs(arcDiff(c.prediction(), st.eta + st.etaDot * dt, L)) < 1e-12);
            CHECK(c.prediction() >= 0.0);
            CHECK(c.prediction() < L);
        }
    }
}

TEST_CASE("Controller<RigidBody>: signs of the tangential and transverse laws (eq:tan-vel, eq:tan-pos, eq:transverse)")
{
    const Tube& T = rigidTube();
    const double L = T.length(), s = 2.0;
    const Pose g = T.path().pose(s);
    Vec6 qdd = Vec6::Zero();
    RigidController::Status st;
    RigidController::OutputState o;
    REQUIRE(RigidController(T, Gains(), TangentialReference::constantVelocity(1.0)).outputs(g, s, o));

    // At rest on the path, behind the reference: eta must accelerate forward, eta_dd = mu_par > 0. At rest
    // Adot qdot = 0, so eta_dd = (A qdd)[0].
    struct Case {
        TangentialReference ref;
        double muPar;  // expected eta_dd
    };
    const Case cases[] = {
        {TangentialReference::constantVelocity(1.0), +5.0},   // e_v = -1: mu_par = -k_t e_v = +5
        {TangentialReference::constantVelocity(-1.0), -5.0},  // and backwards
        {positionConstant(s + 0.5), +4.5},                    // e_eta = -0.5: mu_par = -k_p e_eta = +4.5
        {positionConstant(s - 0.5), -4.5},
        {positionConstant(s + 0.5 - L), +4.5},                // eta_d one loop behind: e_eta wraps to -0.5
        {positionConstant(s - 0.5 + 2.0 * L), -4.5},          // two loops ahead: e_eta wraps to +0.5
    };
    for (const Case& cs : cases) {
        RigidController c(T, Gains(), cs.ref);
        c.reset(s);
        REQUIRE(c.step(g, Vec6::Zero(), 0.0, 0.002, qdd, st));
        const Vec6 yDd = o.A * qdd;
        CAPTURE(cs.muPar, st.eTan, yDd.transpose());
        CHECK_THAT(yDd[0], WithinAbs(cs.muPar, 1e-9));
        CHECK(yDd.tail<5>().cwiseAbs().maxCoeff() < 1e-9);  // on the path xi = 0: no transverse command
    }

    // Moving forward along the path at the reference position: the damping term brakes, mu_par = -k_d eta_dot.
    {
        RigidController c(T, Gains(), positionConstant(s));
        const Vec6 dq = 0.3 * T.path().velocity(s);  // eta_dot = 0.3, xi_dot = 0
        c.reset(s);
        REQUIRE(c.step(g, dq, 0.0, 0.002, qdd, st));
        const Vec6 yDd = o.A * qdd + adotQdot(c, g, dq, s);
        CAPTURE(st.etaDot, yDd.transpose());
        CHECK_THAT(st.etaDot, WithinAbs(0.3, 1e-9));
        CHECK_THAT(yDd[0], WithinAbs(-6.0 * 0.3, 1e-8));
        CHECK(yDd.tail<5>().cwiseAbs().maxCoeff() < 1e-8);
    }

    // Off the path at rest: the transverse command pulls back, xi_dd = -K_p xi.
    {
        RigidController c(T, Gains(), TangentialReference::constantVelocity(0.0));
        c.reset(0.1);
        REQUIRE(c.step(rigidStart(), Vec6::Zero(), 0.0, 0.002, qdd, st));
        RigidController::OutputState os;
        REQUIRE(c.outputs(rigidStart(), 0.1, os));
        const Vec6 yDd = os.A * qdd;
        const double err = maxDiff(Vec5(yDd.tail<5>()), Vec5(-16.0 * rigidXi0()));
        CAPTURE(sci(err));
        CHECK(err < 1e-9);
        CHECK_THAT(yDd[0], WithinAbs(0.0, 1e-9));  // nu_d = 0 and eta_dot = 0: no tangential command
    }
}

TEST_CASE("Controller<RigidBody>: the mu ramp multiplies mu, not Adot qdot; r(0) = 0 and r(T) = 1")
{
    const Tube& T = rigidTube();
    const Gains k = oddGains();
    const double rampT = 2.0;
    RigidController plain(T, k, velocityTest(T.length()));
    RigidController ramped(T, k, velocityTest(T.length()));
    ramped.setRampTime(rampT);
    const Pose g = rigidStart();

    auto run = [&](RigidController& c, const Vec6& dq, double t, Vec6& qdd, RigidController::Status& st) {
        c.reset(0.1);
        REQUIRE(c.step(g, dq, t, 0.002, qdd, st));
    };
    Vec6 q0 = Vec6::Zero(), q1 = Vec6::Zero();
    RigidController::Status s0, s1;

    // At rest at t = 0: the first command is exactly zero, so the input starts continuously from rest.
    run(ramped, Vec6::Zero(), 0.0, q0, s0);
    CHECK(s0.ramp == 0.0);
    CHECK((q0.array() == 0.0).all());

    // At t = T the ramp is 1: the same command as without a ramp. Halfway, r = 1/2 and at rest qdd* is halved.
    run(ramped, Vec6::Zero(), rampT, q0, s0);
    run(plain, Vec6::Zero(), rampT, q1, s1);
    CHECK(s0.ramp == 1.0);
    CHECK(maxDiff(q0, q1) < 1e-15);
    run(ramped, Vec6::Zero(), 0.5 * rampT, q0, s0);
    run(plain, Vec6::Zero(), 0.5 * rampT, q1, s1);
    CHECK_THAT(s0.ramp, WithinAbs(0.5, 1e-15));
    CHECK(maxDiff(q0, Vec6(0.5 * q1)) < 1e-12);
    CHECK(s1.ramp == 1.0);  // no ramp: r = 1 always

    // Moving at t = 0: mu is switched off, Adot qdot is not, so A qdd* = -Adot qdot (not zero).
    const Vec6 dq = randn<Vec6>(0.7);
    run(ramped, dq, 0.0, q0, s0);
    RigidController::OutputState o;
    REQUIRE(ramped.outputs(g, 0.1, o));
    const Vec6 Adq = adotQdot(ramped, g, dq, o.eta);
    const double err = maxDiff(Vec6(o.A * q0), Vec6(-Adq));
    CAPTURE(sci(err), Adq.transpose());
    CHECK(err < 1e-8);
    CHECK(Adq.cwiseAbs().maxCoeff() > 0.01);
}

TEST_CASE("Controller<RigidBody>: failure is false, not an exception, and leaves qdd unchanged")
{
    const Tube& T = rigidTube();
    RigidController c(T, Gains(), TangentialReference::constantVelocity(0.5));
    const Vec6 marker = Vec6::Constant(7.0);
    Vec6 qdd = marker;
    RigidController::Status st;

    Pose bad = rigidStart();
    bad(0, 3) = std::numeric_limits<double>::quiet_NaN();
    c.reset(0.1);
    CHECK_FALSE(c.step(bad, Vec6::Zero(), 0.0, 0.002, qdd, st));
    CHECK_FALSE(st.ok);
    CHECK(qdd == marker);

    Vec6 dq = Vec6::Zero();
    dq[2] = std::numeric_limits<double>::infinity();
    c.reset(0.1);
    CHECK_FALSE(c.step(rigidStart(), dq, 0.0, 0.002, qdd, st));
    CHECK(qdd == marker);

    RigidController::OutputState o;
    CHECK_FALSE(c.outputs(bad, 0.1, o));

    // A reference that returns NaN makes mu, and so qdd*, non-finite: refused after the solve
    RigidController nanRef(T, Gains(), TangentialReference::constantVelocity(std::numeric_limits<double>::quiet_NaN()));
    nanRef.reset(0.1);
    CHECK_FALSE(nanRef.step(rigidStart(), Vec6::Zero(), 0.0, 0.002, qdd, st));
    CHECK_FALSE(st.ok);
    CHECK(qdd == marker);
}

TEST_CASE("Controller: a nearly singular A (sigma_min < 1e-9) is refused, for any System")
{
    // StuckBody can hardly move along its z axis, so A = J_Psi^-1 J_b is nearly singular. The LU solve alone would
    // return a huge but finite qdd*; the sigma_min check must refuse it.
    const Tube& T = rigidTube();
    Controller<StuckBody> c(T, Gains(), TangentialReference::constantVelocity(0.5));
    Controller<StuckBody>::OutputState o;
    REQUIRE(c.outputs(rigidStart(), 0.1, o));
    const double sigma = JacobiSVD<Mat6>(o.A).singularValues().minCoeff();
    CAPTURE(sci(sigma));
    CHECK(sigma < 1e-11);
    CHECK(sigma > 0.0);
    CHECK(Vec6(o.A.partialPivLu().solve(Vec6::Ones())).allFinite());  // LU does not notice

    Vec6 qdd = Vec6::Constant(7.0);
    Controller<StuckBody>::Status st;
    c.reset(0.1);
    CHECK_FALSE(c.step(rigidStart(), Vec6::Zero(), 0.0, 0.002, qdd, st));
    CHECK_FALSE(st.ok);
    CHECK(qdd == Vec6::Constant(7.0));
}

TEST_CASE("Controller<RigidBody>: step does not allocate")
{
    // States along a short closed loop, computed before the window (the plant and the log are not under test).
    const Tube& T = rigidTube();
    const RigidBody body = RigidBody::example();
    RigidController c(T, Gains(), TangentialReference::constantVelocity(0.5));
    std::vector<Pose> poses;
    std::vector<Vec6> twists;
    poses.reserve(10);
    twists.reserve(10);
    {
        Pose g = rigidStart();
        Vec6 zeta = Vec6::Zero();
        c.reset(0.1);
        for (int k = 0; k < 10; ++k) {
            Vec6 qdd = Vec6::Zero();
            RigidController::Status st;
            REQUIRE(c.step(g, zeta, 0.002 * k, 0.002, qdd, st));
            const Vec6 u = body.mass() * qdd + body.coriolis(zeta);
            for (int i = 0; i < 4; ++i)
                body.integrate(g, zeta, u, 0.0005);
            poses.push_back(g);
            twists.push_back(zeta);
        }
    }
    REQUIRE(twists.back().norm() > 0.01);  // moving: the two extra projections for Adot qdot run

    c.reset(0.1);
    bool allOk = true;
    double sink = 0.0;
    Vec6 qdd = Vec6::Zero();
    RigidController::Status st;
    // No Catch2 macros inside the window: they allocate.
    const long before = gNewCalls.load(std::memory_order_relaxed);
#ifdef EIGEN_RUNTIME_NO_MALLOC
    Eigen::internal::set_is_malloc_allowed(false);
#endif
    for (int k = 0; k < 1000; ++k) {
        allOk = c.step(poses[k % 10], twists[k % 10], 0.002 * k, 0.002, qdd, st) && allOk;
        sink += qdd[0];
    }
#ifdef EIGEN_RUNTIME_NO_MALLOC
    Eigen::internal::set_is_malloc_allowed(true);
#endif
    const long after = gNewCalls.load(std::memory_order_relaxed);
    CAPTURE(before, after, sink);
    CHECK(allOk);
    REQUIRE(after == before);
}

TEST_CASE("Timing: Controller<RigidBody>::step", "[.timing]")
{
    // Release build only: ./build-rel/tests "[.timing]"
#ifndef NDEBUG
    SKIP("timing is meaningful in a Release build only");
#endif
    const Tube& T = rigidTube();
    RigidController c(T, Gains(), TangentialReference::constantVelocity(0.5));
    const Pose g = rigidStart();
    const Vec6 dq = randn<Vec6>(0.5);
    Vec6 qdd = Vec6::Zero();
    RigidController::Status st;
    const int calls = 20000;
    double sink = 0.0;
    int failures = 0;

    const auto start = std::chrono::steady_clock::now();
    for (int k = 0; k < calls; ++k) {
        c.reset(0.1);  // the same state every time: one Newton update per projection, as in the closed loop
        if (!c.step(g, dq, 0.002 * k, 0.002, qdd, st))
            ++failures;
        sink += qdd[k % 6];
    }
    const auto stop = std::chrono::steady_clock::now();

    const double us = std::chrono::duration<double, std::micro>(stop - start).count() / calls;
    WARN("Controller<RigidBody>::step (moving, three projections): " << us << " us per call (checksum " << sink
         << ")");
    CHECK(failures == 0);
    CHECK(us < 50.0);  // a few microseconds; the 1 kHz budget is 1000 us
}

// ===============================================================================================================
// Panda (n = 7), with two simple redundancy outputs: the coupling beta = q3 - q1, and the test-only output that adds
// a sin(k eta), which checks that eta and d eta reach eq:dbeta. The tests with the lift's beta of Lemma 3
// (eq:beta-flat) follow below; the Panda closed loops are in test_panda_example.cpp.

TEST_CASE("Controller<PandaArm>: compiles, and needs a redundancy output")
{
    const Tube& T = pandaTube();
    CHECK_THROWS_AS(PandaController(T, Gains(), TangentialReference::constantVelocity(0.08), nullptr),
                    std::invalid_argument);
    CHECK_THROWS_AS(PandaController(T, Gains(), TangentialReference::constantVelocity(0.08)), std::invalid_argument);
    CHECK_NOTHROW(PandaController(T, Gains(), TangentialReference::constantVelocity(0.08), &coupling()));
    const PandaController c(T, Gains(), TangentialReference::constantVelocity(0.08), &coupling());
    CHECK(c.redundancy() == &coupling());
    STATIC_REQUIRE(PandaController::n == 7);
    STATIC_REQUIRE(std::is_same_v<PandaController::VecN, Vec7>);
    STATIC_REQUIRE(std::is_same_v<RigidController::Config, Pose>);
}

TEST_CASE("Controller<PandaArm>: A is the Jacobian of ybar = (eta, xi, beta) off the path (eq:A)")
{
    const Tube& T = pandaTube();
    const PandaController c(T, Gains(), TangentialReference::constantVelocity(0.08), &coupling());
    const Vec7& q = pandaStart();

    // The start configuration has the requested outputs: eta = 0.1, xi = xi0, beta = 0.03
    PandaController::OutputState o;
    REQUIRE(c.outputs(q, 0.1, o));
    CAPTURE(o.eta, o.xi.transpose(), o.beta);
    CHECK_THAT(o.eta, WithinAbs(0.1, 1e-9));
    CHECK(maxDiff(o.xi, pandaXi0()) < 1e-9);
    CHECK_THAT(o.beta, WithinAbs(0.03, 1e-12));
    CHECK(o.A.row(6) == (Vec7::Unit(2) - Vec7::Unit(0)).transpose());  // d beta of the coupling, as given

    // Top rows: J_Psi^-1 J_b
    Mat67 Jb = Mat67::Zero();
    PandaArm::output(q, &Jb);
    const double errTop = maxDiff(Mat67(T.JPsi(o.eta, o.xi) * o.A.topRows<6>()), Jb);
    CAPTURE(sci(errTop));
    CHECK(errTop < 1e-12);

    checkJacobian(c, q, 0.1, 1.0);
}

TEST_CASE("Controller<PandaArm>: step returns qdd* with A qdd* + Adot qdot = mu, including mu_beta (eq:beta-law)")
{
    const Tube& T = pandaTube();
    const double L = T.length();
    const Gains k = oddGains();
    PandaController c(T, k, velocityTest(L), &coupling());
    const Vec7& q = pandaStart();

    for (int i = 0; i < 3; ++i) {
        Vec7 dq = randn<Vec7>(0.2);
        dq[2] = dq[0] + 0.1;  // beta_dot = 0.1, so mu_beta = -13 (0.03) - 3 (0.1) = -0.69
        const double t = uniform(0.0, 10.0);
        c.reset(0.1);
        Vec7 qdd = Vec7::Zero();
        PandaController::Status st;
        REQUIRE(c.step(q, dq, t, 0.001, qdd, st));
        CHECK_THAT(st.beta, WithinAbs(0.03, 1e-12));
        CHECK_THAT(st.betaDot, WithinAbs(dq[2] - dq[0], 1e-12));
        CHECK(st.sigmaMin > 1e-3);

        PandaController::OutputState o;
        REQUIRE(c.outputs(q, 0.1, o));
        const Vec7 Adq = adotQdot(c, q, dq, o.eta);
        const Vec7 mu =
            expectedMu<PandaArm>(k, st, true, nuTest(t, st.eta, L), dnuDtTest(t), dnuDetaTest(st.eta, L), L);
        const double errMu = maxDiff(Vec7(o.A * qdd + Adq), mu);
        CAPTURE(t, sci(errMu), mu.transpose(), Adq.transpose());
        CHECK(errMu < 1e-8);
        CHECK_THAT(mu[6], WithinAbs(-0.69, 1e-12));  // the beta channel is active
    }

    // At rest with the ramp at t = 0: zero command, for n = 7 too
    c.setRampTime(1.0);
    c.reset(0.1);
    Vec7 qdd = Vec7::Constant(1.0);
    PandaController::Status st;
    REQUIRE(c.step(q, Vec7::Zero(), 0.0, 0.001, qdd, st));
    CHECK((qdd.array() == 0.0).all());
}

TEST_CASE("Controller<PandaArm>: a beta that depends on eta gets the projected eta and d eta (eq:beta-flat, eq:dbeta)")
{
    // With beta = q3 - q1 + a sin(k eta), the last row of A is d beta = e3 - e1 + a k cos(k eta) d eta, where d eta
    // is row 0 of J_Psi^-1 J_b. Like the lift's beta of Lemma 3, it needs the eta and d eta that outputs computes.
    const Tube& T = pandaTube();
    const WavyCouplingOutput& wavy = wavyCoupling();
    const PandaController c(T, Gains(), TangentialReference::constantVelocity(0.08), &wavy);
    const Vec7& q = pandaStart();

    // The guess is off by 0.02, so a beta evaluated at the guess instead of at eta would show.
    PandaController::OutputState o;
    REQUIRE(c.outputs(q, 0.12, o));
    CHECK_THAT(o.eta, WithinAbs(0.1, 1e-9));
    const double betaExpected = q[2] - q[0] + WavyCouplingOutput::wave(o.eta);
    const Vec7 deta = o.A.row(0).transpose();
    const Vec7 dbeta = o.A.row(6).transpose();
    const Vec7 dbetaExpected = Vec7::Unit(2) - Vec7::Unit(0) + WavyCouplingOutput::waveSlope(o.eta) * deta;
    const double errDbeta = maxDiff(dbeta, dbetaExpected);
    CAPTURE(o.beta, betaExpected, sci(errDbeta), deta.transpose());
    CHECK_THAT(o.beta, WithinAbs(betaExpected, 1e-12));
    CHECK(errDbeta < 1e-12);
    // The d eta term is not small: a k cos(0.8) = 1.11 times d eta, whose largest entry is about 0.4
    CHECK(maxDiff(dbeta, Vec7(Vec7::Unit(2) - Vec7::Unit(0))) > 0.3);

    // and A is the Jacobian of ybar = (eta, xi, beta), where beta changes with eta too
    checkJacobian(c, q, 0.1, 1.0);
}

TEST_CASE("Controller<PandaArm>: with a beta that depends on eta, step still gives A qdd* + Adot qdot = mu")
{
    const Tube& T = pandaTube();
    const double L = T.length();
    const Gains k = oddGains();
    const WavyCouplingOutput& wavy = wavyCoupling();
    PandaController c(T, k, velocityTest(L), &wavy);
    const Vec7& q = pandaStart();

    for (int i = 0; i < 3; ++i) {
        const Vec7 dq = randn<Vec7>(0.2);
        const double t = uniform(0.0, 10.0);
        c.reset(0.11);  // the prediction differs from eta = 0.1, as in the closed loop
        Vec7 qdd = Vec7::Zero();
        PandaController::Status st;
        REQUIRE(c.step(q, dq, t, 0.001, qdd, st));

        // beta at the projected eta, and beta_dot = (e3 - e1)^T qdot + a k cos(k eta) eta_dot
        CHECK_THAT(st.eta, WithinAbs(0.1, 1e-9));
        CHECK_THAT(st.beta, WithinAbs(q[2] - q[0] + WavyCouplingOutput::wave(st.eta), 1e-12));
        CHECK_THAT(st.betaDot, WithinAbs(dq[2] - dq[0] + WavyCouplingOutput::waveSlope(st.eta) * st.etaDot, 1e-12));

        PandaController::OutputState o;
        REQUIRE(c.outputs(q, 0.1, o));
        const Vec7 Adq = adotQdot(c, q, dq, o.eta);
        const Vec7 mu =
            expectedMu<PandaArm>(k, st, true, nuTest(t, st.eta, L), dnuDtTest(t), dnuDetaTest(st.eta, L), L);
        const double errMu = maxDiff(Vec7(o.A * qdd + Adq), mu);
        CAPTURE(t, sci(errMu), mu.transpose(), Adq.transpose());
        CHECK(errMu < 1e-8);
        // The beta row of Adot qdot is zero for the coupling; here the eta terms make it active: it is
        // a k (-k sin(k eta) eta_dot^2 + cos(k eta) (Adot qdot)_0), whose two terms can nearly cancel: 0.28, 0.43 and
        // 0.006 for these three velocities (as low as 7e-4 for other draws). The bound 1e-4 is still 1e4 times the
        // bound on errMu, so a controller that dropped this row would fail the check above.
        CHECK(std::abs(Adq[6]) > 1e-4);
    }
}

// ===============================================================================================================
// Panda (n = 7) with the lift's beta = X(eta)^T (q - q*(eta)) of Lemma 3 (eq:beta-flat, eq:dbeta): the paper's
// example (testing::pandaExample(), panda_fixture.hpp), with its tube, rest-posture lift, kernel frame and
// LiftRedundancy. The closed loops themselves are in test_panda_example.cpp; this file checks only that they
// allocate before their loop, with its counting operator new.

TEST_CASE("Controller<PandaArm> with the lift's beta: A is the Jacobian of ybar off the lift, at both starts (eq:A)")
{
    // Off the lift q - q*(eta) != 0, so the X' term of eq:dbeta is active; on the lift it drops out, and a sign
    // error in it would go unseen. At the paper start (|q - q*| = 0.42 rad) the term changes the last row of A by
    // up to 0.31; at the self-motion start (q - q* = 0.2 X + O(0.2^2), and X'^T X = 0 since |X| = 1) by 1.3e-3.
    const PandaExample& ex = pandaExample();
    const PandaController c(ex.tube(), PandaExample::gains(), PandaExample::reference(), &ex.redundancy());
    struct Start {
        Vec7 q;
        Vec5 xi;
        double beta;
        double xPrimeTerm;  // a lower bound of the X' term's largest entry
    };
    const Start starts[] = {
        {ex.paperStart(), PandaExample::xi0(), 0.0, 0.3},
        {ex.selfMotionStart(), Vec5::Zero(), PandaExample::selfMotionBeta0, 1e-3},
    };
    for (const Start& start : starts) {
        CAPTURE(start.beta);
        PandaController::OutputState o;
        REQUIRE(c.outputs(start.q, 0.1, o));
        CHECK_THAT(o.eta, WithinAbs(0.1, 1e-10));
        CHECK(maxDiff(o.xi, start.xi) < 1e-10);
        CHECK_THAT(o.beta, WithinAbs(start.beta, 1e-10));

        // The last row of A is eq:dbeta at the projected eta, with d eta = row 0 of J_Psi^-1 J_b
        Vec7 qs = Vec7::Zero(), dqs = Vec7::Zero(), X = Vec7::Zero(), dX = Vec7::Zero();
        ex.lift().eval(o.eta, qs, &dqs);
        ex.kernelFrame().eval(o.eta, X, &dX);
        const Vec7 deta = o.A.row(0).transpose();
        const Vec7 dbeta = o.A.row(6).transpose();
        const Vec7 xPrimeTerm = dX.dot(start.q - qs) * deta;
        const double errDbeta = maxDiff(dbeta, Vec7(X + xPrimeTerm - X.dot(dqs) * deta));
        CAPTURE(sci(errDbeta), xPrimeTerm.cwiseAbs().maxCoeff());
        CHECK(errDbeta < 1e-13);
        CHECK(xPrimeTerm.cwiseAbs().maxCoeff() > start.xPrimeTerm);

        // Central differences of (eta, xi, beta) along random joint velocities, h = 1e-6, eta wrapped
        checkJacobian(c, start.q, 0.1, 1.0);
    }
}

TEST_CASE("Controller<PandaArm> with the lift's beta: on the lift xi = 0, beta = 0 and A is invertible (Proposition 2)")
{
    // On the lift J_b(q*) X = 0 and d beta X = X^T X = 1 (the d eta term vanishes with J_b X), so A X = e_7: the
    // self-motion moves beta alone. The smallest singular value of A along the lift is 0.054 (largest 0.084).
    const PandaExample& ex = pandaExample();
    const PandaController c(ex.tube(), PandaExample::gains(), PandaExample::reference(), &ex.redundancy());
    const double L = ex.tube().length();
    std::vector<double> points = {0.0, 0.1, 0.5 * L, L - 1e-7};  // the seam s = 0 = L too
    for (int i = 0; i < 16; ++i)
        points.push_back(uniform(0.0, L));
    double sigmaLow = std::numeric_limits<double>::infinity();
    for (double s : points) {
        Vec7 q = Vec7::Zero(), X = Vec7::Zero();
        ex.lift().eval(s, q);
        ex.kernelFrame().eval(s, X);
        PandaController::OutputState o;
        REQUIRE(c.outputs(q, s + 0.01, o));
        const double sigma = JacobiSVD<Mat<7, 7>>(o.A).singularValues().minCoeff();
        const double errAX = maxDiff(Vec7(o.A * X), Vec7(Vec7::Unit(6)));
        CAPTURE(s, sci(o.xi.norm()), sci(o.beta), sci(errAX), sigma);
        CHECK(std::abs(arcDiff(o.eta, s, L)) < 1e-10);
        CHECK(o.xi.norm() < 1e-10);
        CHECK(std::abs(o.beta) < 1e-10);
        CHECK(errAX < 1e-10);
        CHECK(sigma > 0.05);
        sigmaLow = std::min(sigmaLow, sigma);
    }
    CHECK(sigmaLow < 0.07);  // a real test: sigma_min(A) is not of order one here
}

TEST_CASE("Controller<PandaArm> with the lift's beta: step does not allocate")
{
    // States along a short closed loop from the self-motion start, which moves from the first step (A(q0) qd0 =
    // (nu_d, 0, ..., 0)), computed before the window. Every step runs the two extra projections for Adot qdot and
    // evaluates the lift and kernel-frame splines for beta (LiftRedundancy::evaluate).
    const PandaExample& ex = pandaExample();
    PandaController c(ex.tube(), PandaExample::gains(), PandaExample::reference(), &ex.redundancy());
    std::vector<Vec7> qs, qds;
    qs.reserve(10);
    qds.reserve(10);
    {
        Vec7 q = ex.selfMotionStart();
        PandaController::OutputState o;
        REQUIRE(c.outputs(q, 0.1, o));
        Vec7 qd = o.A.partialPivLu().solve(Vec7(0.08 * Vec7::Unit(0)));
        c.reset(0.1);
        for (int k = 0; k < 10; ++k) {
            Vec7 qdd = Vec7::Zero();
            PandaController::Status st;
            REQUIRE(c.step(q, qd, 0.001 * k, 0.001, qdd, st));
            q += 0.001 * qd + (0.5e-6) * qdd;
            qd += 0.001 * qdd;
            qs.push_back(q);
            qds.push_back(qd);
        }
    }
    REQUIRE(qds.back().norm() > 0.1);  // moving

    c.reset(0.1);
    bool allOk = true;
    double sink = 0.0;
    Vec7 qdd = Vec7::Zero();
    PandaController::Status st;
    // No Catch2 macros inside the window: they allocate.
    const long before = gNewCalls.load(std::memory_order_relaxed);
#ifdef EIGEN_RUNTIME_NO_MALLOC
    Eigen::internal::set_is_malloc_allowed(false);
#endif
    for (int k = 0; k < 1000; ++k) {
        allOk = c.step(qs[k % 10], qds[k % 10], 0.001 * k, 0.001, qdd, st) && allOk;
        sink += qdd[0];
    }
#ifdef EIGEN_RUNTIME_NO_MALLOC
    Eigen::internal::set_is_malloc_allowed(true);
#endif
    const long after = gNewCalls.load(std::memory_order_relaxed);
    CAPTURE(before, after, sink);
    CHECK(allOk);
    REQUIRE(after == before);
}

TEST_CASE("simulate (Panda): the closed loop allocates only before its loop, so a longer run allocates no more")
{
    // The closed loops of test_panda_example.cpp seen through the counting operator new: a run's one allocation is
    // the log's reserve before the loop, so a 51-step run allocates as often as a one-step run. A log that grew row
    // by row inside the loop would allocate 7 times for 51 rows (capacity 1, 2, 4, ..., 64). Eigen's own heap use
    // in the loop would be the step's, which the test above checks with EIGEN_RUNTIME_NO_MALLOC.
    const PandaExample& ex = pandaExample();
    for (PandaExample::RunKind kind : {PandaExample::RunKind::Paper, PandaExample::RunKind::SelfMotion}) {
        const long before = gNewCalls.load(std::memory_order_relaxed);
        const PandaRun one = simulate(ex, kind, 0.0);
        const long middle = gNewCalls.load(std::memory_order_relaxed);
        const PandaRun longer = simulate(ex, kind, 0.05);
        const long after = gNewCalls.load(std::memory_order_relaxed);
        CAPTURE(middle - before, after - middle);
        REQUIRE(one.ok);
        REQUIRE(longer.ok);
        REQUIRE(one.log.size() == 1);
        REQUIRE(longer.log.size() == 51);
        CHECK(after - middle == middle - before);
        CHECK(middle - before == 1);  // the reserve: the controller and the setup allocate nothing
    }
}

TEST_CASE("Timing: Controller<PandaArm>::step with the lift's beta", "[.timing]")
{
    // Release build only: ./build-rel/tests "[.timing]". The paper reports 27 us per cycle (2.1 GHz Xeon), the
    // reference implementation about 24 us; the 1 kHz budget is 1000 us, and the design doc asks for < 100 us.
#ifndef NDEBUG
    SKIP("timing is meaningful in a Release build only");
#endif
    const PandaExample& ex = pandaExample();
    PandaController c(ex.tube(), PandaExample::gains(), PandaExample::reference(), &ex.redundancy());
    const Vec7& q = ex.paperStart();
    const Vec7 dq = randn<Vec7>(0.2);
    Vec7 qdd = Vec7::Zero();
    PandaController::Status st;
    const int calls = 20000;
    double sink = 0.0;
    int failures = 0, maxIterations = 0;

    const auto start = std::chrono::steady_clock::now();
    for (int k = 0; k < calls; ++k) {
        c.reset(0.1);  // the same state every time: one Newton update per projection, as in the closed loop
        if (!c.step(q, dq, 0.001 * k, 0.001, qdd, st))
            ++failures;
        maxIterations = std::max(maxIterations, st.iterations);
        sink += qdd[k % 7];
    }
    const auto stop = std::chrono::steady_clock::now();

    const double us = std::chrono::duration<double, std::micro>(stop - start).count() / calls;
    WARN("Controller<PandaArm>::step (moving, three projections, lift's beta): " << us
         << " us per call (paper 27 us; checksum " << sink << ")");
    CHECK(failures == 0);
    CHECK(maxIterations == 1);
    CHECK(us < 100.0);
}

namespace {
    // The Panda closed loop at the joint-acceleration level with a constant disturbance d added to every commanded
    // qdd* (a joint torque the model does not know, divided by the inertia: friction, a payload), from rest on the
    // path at s0, with the reference given (hold at s0, or a constant nu_d). Returns the last status and the largest
    // |xi| on the way; stops early if a step fails.
    bool disturbedLoop(const PandaExample& ex, const Gains& gains, const TangentialReference& ref, const Vec7& d,
                       double T, Controller<PandaArm>::Status& st, double* maxXi = nullptr)
    {
        Controller<PandaArm> c(ex.tube(), gains, ref, &ex.redundancy());
        c.reset(PandaExample::s0);
        Vec7 q = Vec7::Zero(), dq = Vec7::Zero();
        ex.lift().eval(PandaExample::s0, q);
        const double dt = 1e-3;
        double mx = 0.0;
        for (long k = 0; k * dt <= T; ++k) {
            Vec7 qdd = Vec7::Zero();
            if (!c.step(q, dq, k * dt, dt, qdd, st))
                return false;
            mx = std::max(mx, st.xi.norm());
            qdd += d;
            q += dt * dq + (0.5 * dt * dt) * qdd;
            dq += dt * qdd;
        }
        if (maxXi)
            *maxXi = mx;
        return true;
    }

    TangentialReference holdAtS0()
    {
        return TangentialReference::fromPosition([](double, double& e, double& de, double& dde) {
            e = PandaExample::s0;
            de = 0.0;
            dde = 0.0;
        });
    }

    // A d in output coordinates at the configuration of the status: what the integral terms must carry
    Vec7 outputDisturbance(const PandaExample& ex, const Gains& gains, const Vec7& q, double eta, const Vec7& d)
    {
        Controller<PandaArm>::OutputState out;
        REQUIRE(Controller<PandaArm>(ex.tube(), gains, holdAtS0(), &ex.redundancy()).outputs(q, eta, out));
        return out.A * d;
    }
}  // namespace

TEST_CASE("Integral action: a constant disturbance leaves a steady error without it and none with it (CTRL-1)")
{
    const PandaExample& ex = pandaExample();
    const Vec7 d = (Vec7() << 0.1, 0.0, 0.1, 0.0, 0.3, 0.0, 0.3).finished();  // rad/s^2
    Gains pd = PandaExample::gains();
    Gains pid = pd;
    pid.kiEta = 25.0;
    pid.Ki = 36.0;
    pid.KiBeta = 36.0;
    Controller<PandaArm>::Status st;

    SECTION("holding a point of the path (A constant): the PID removes the error, its integrals carry A d")
    {
        REQUIRE(disturbedLoop(ex, pd, holdAtS0(), d, 12.0, st));
        const double xiPd = st.xi.norm(), betaPd = std::abs(st.beta), ePd = std::abs(st.eTan);
        CAPTURE(sci(xiPd), sci(betaPd), sci(ePd));
        CHECK(xiPd > 1e-3);  // the paper's PD: a steady error, A d / K_p
        CHECK(st.iXi.isZero());  // and no integral terms
        CHECK(st.iBeta == 0.0);
        CHECK(st.iEta == 0.0);
        // The error is the disturbance over the gains (the hold's k_p = 9 on eta, K_p = 16 on xi and beta)
        Vec7 q = Vec7::Zero();
        ex.lift().eval(PandaExample::s0, q);
        const Vec7 Ad = outputDisturbance(ex, pd, q, PandaExample::s0, d);
        CHECK_THAT(st.eTan, WithinAbs(Ad[0] / 9.0, 0.1 * std::abs(Ad[0] / 9.0)));
        CHECK_THAT(st.xi.norm(), WithinAbs(Ad.segment<5>(1).norm() / 16.0, 0.1 * Ad.segment<5>(1).norm() / 16.0));

        // 20 s: with the hold's k_p = 9, k_d = 6 the eta channel's integral mode (s^3 + 6 s^2 + 9 s + 25) is lightly
        // damped at -0.4 +- 2.2i, so it takes that long to settle to 1e-5
        double maxXi = 0.0;
        REQUIRE(disturbedLoop(ex, pid, holdAtS0(), d, 20.0, st, &maxXi));
        CAPTURE(sci(st.xi.norm()), sci(st.beta), sci(st.eTan), maxXi);
        CHECK(st.xi.norm() < 1e-5);
        CHECK(std::abs(st.beta) < 1e-5);
        CHECK(std::abs(st.eTan) < 1e-5);
        CHECK(maxXi < 0.02);  // no large transient on the way
        // The integral terms are the disturbance in output coordinates
        const Vec7 AdEnd = outputDisturbance(ex, pid, q, PandaExample::s0, d);
        CHECK_THAT(st.iEta, WithinAbs(AdEnd[0], 1e-4));
        CHECK(maxDiff(st.iXi, Vec5(AdEnd.segment<5>(1))) < 1e-4);
        CHECK_THAT(st.iBeta, WithinAbs(AdEnd[6], 1e-4));
    }
    SECTION("moving along the path (A d varies with eta): the PID keeps the errors well below the PD's")
    {
        REQUIRE(disturbedLoop(ex, pd, PandaExample::reference(), d, 12.0, st));
        const double xiPd = st.xi.norm(), ePd = std::abs(st.eTan);
        REQUIRE(disturbedLoop(ex, pid, PandaExample::reference(), d, 12.0, st));
        CAPTURE(sci(xiPd), sci(ePd), sci(st.xi.norm()), sci(st.beta), sci(st.eTan));
        CHECK(xiPd > 3e-3);
        CHECK(st.xi.norm() < 0.2 * xiPd);
        CHECK(std::abs(st.eTan) < 0.5 * ePd);
        CHECK(std::abs(st.beta) < 2e-3);
        CHECK(st.iXi.norm() > 1e-3);  // the integrals carry the (slowly varying) disturbance
        CHECK(std::abs(st.iEta) > 1e-3);
    }
}

TEST_CASE("Integral action: the clamps hold, the integrals reset, bad gains are refused (CTRL-1)")
{
    const PandaExample& ex = pandaExample();
    Gains pid = PandaExample::gains();
    pid.kiEta = 25.0;
    pid.Ki = 36.0;
    pid.KiBeta = 36.0;
    pid.integralMaxEta = 0.05;  // far below what the disturbance needs: every integral runs into its clamp
    pid.integralMaxXi = 0.05;
    pid.integralMaxBeta = 0.005;
    const Vec7 d = (Vec7() << 0.2, 0.0, 0.2, 0.0, 0.6, 0.0, 0.6).finished();
    Controller<PandaArm>::Status st;
    double maxXi = 0.0;
    REQUIRE(disturbedLoop(ex, pid, holdAtS0(), d, 10.0, st, &maxXi));
    CAPTURE(st.iEta, st.iXi.transpose(), st.iBeta, maxXi);
    CHECK_THAT(std::abs(st.iEta), WithinAbs(pid.integralMaxEta, 1e-12));
    CHECK_THAT(std::abs(st.iBeta), WithinAbs(pid.integralMaxBeta, 1e-12));
    CHECK(st.iXi.cwiseAbs().maxCoeff() <= pid.integralMaxXi + 1e-12);
    CHECK_THAT(st.iXi.cwiseAbs().maxCoeff(), WithinAbs(pid.integralMaxXi, 1e-12));
    CHECK(maxXi < 0.04);  // the PD still holds the arm in the tube with the integrals saturated
    CHECK(st.xi.norm() > 1e-4);  // ... with a steady error left, since the integrals are clamped
    // Reset zeroes the integrals: the next step at rest on the lift carries at most one period of error
    Controller<PandaArm> c(ex.tube(), pid, PandaExample::reference(), &ex.redundancy());
    Vec7 q = Vec7::Zero(), qdd = Vec7::Zero();
    ex.lift().eval(PandaExample::s0, q);
    for (int k = 0; k < 3000; ++k)  // wind the integrals up: at rest, e_v = -nu_d every period
        REQUIRE(c.step(q, Vec7::Zero(), 1e-3 * k, 1e-3, qdd, st));
    CHECK_THAT(std::abs(st.iEta), WithinAbs(pid.integralMaxEta, 1e-12));
    c.reset(PandaExample::s0);
    REQUIRE(c.step(q, Vec7::Zero(), 3.0, 1e-3, qdd, st));
    CHECK_THAT(st.iEta, WithinAbs(-pid.kiEta * PandaExample::nuD * 1e-3, 1e-12));  // one period of e_v = -nu_d
    CHECK(st.iXi.norm() < 1e-9);
    CHECK(std::abs(st.iBeta) < 1e-9);
    // Refused: negative or non-finite gains, bounds that are not positive
    Gains bad = pid;
    bad.Ki = -1.0;
    CHECK_THROWS_AS(Controller<PandaArm>(ex.tube(), bad, PandaExample::reference(), &ex.redundancy()),
                    std::invalid_argument);
    bad = pid;
    bad.KiBeta = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(Controller<PandaArm>(ex.tube(), bad, PandaExample::reference(), &ex.redundancy()),
                    std::invalid_argument);
    bad = pid;
    bad.integralMaxXi = 0.0;
    CHECK_THROWS_AS(Controller<PandaArm>(ex.tube(), bad, PandaExample::reference(), &ex.redundancy()),
                    std::invalid_argument);
    bad = pid;
    bad.integralMaxEta = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(Controller<PandaArm>(ex.tube(), bad, PandaExample::reference(), &ex.redundancy()),
                    std::invalid_argument);
}

TEST_CASE("ExperimentLoop::follow, observe, the safety stop and the safety checks do not allocate")
{
    // What the robot driver's 1 kHz callback calls (apps/robot/robot_driver.cpp; Rules for the 1 kHz loop): one
    // control cycle of an experiment with the model terms, the safety stop's torque and the outputs of the stop,
    // and each check of the safety monitor directly. States along a short run of the paper experiment, moving.
    const PandaExample& ex = pandaExample();
    const ExperimentSpec* spec = findExperiment("paper");
    REQUIRE(spec != nullptr);
    const ExperimentSetup setup(ex, *spec);
    const SafetyMonitor monitor = setup.monitor();
    std::vector<Vec7> qs, qds;
    qs.reserve(10);
    qds.reserve(10);
    {
        ExperimentLoop warm(setup, monitor);
        Vec7 q = setup.start(), qd = Vec7::Zero();
        ExperimentLoop::Cycle c;
        for (int k = 0; k < 1500; ++k) {
            warm.follow(0.001 * k, q, qd, nullptr, c);
            REQUIRE(c.violation == Violation::None);
            q += 0.001 * qd + 0.5e-6 * c.qdd;
            qd += 0.001 * c.qdd;
            if (k >= 1490) {
                qs.push_back(q);
                qds.push_back(qd);
            }
        }
    }
    REQUIRE(qds.back().norm() > 0.01);  // moving

    ExperimentLoop loop(setup, monitor);
    // ... and the same cycle with the friction feedforward on (the lift's spline and tanh at every cycle: CTRL-1)
    ExperimentSpec withFriction = *spec;
    withFriction.friction = FrictionFeedforward::gaz2019();
    withFriction.gains.Ki = withFriction.gains.KiBeta = withFriction.gains.kiEta = 36.0;  // the integral path too
    const ExperimentSetup setupF(ex, withFriction);
    ExperimentLoop loopF(setupF, monitor);
    ModelTerms terms;
    terms.M = 2.0 * Mat<7, 7>::Identity();
    terms.c = Vec7::Constant(0.01);
    terms.G = Vec7::Constant(1.0);
    ExperimentLoop::Cycle c;
    ExperimentLoop::Status st;
    int violations = 0;
    double sink = 0.0;
    const Vec5 xi = Vec5::Constant(0.001);
    SafetyStop stop;
    stop.begin(1.49, qs[0]);
    bool stopDone = false;
    // No Catch2 macros inside the window: they allocate.
    const long before = gNewCalls.load(std::memory_order_relaxed);
#ifdef EIGEN_RUNTIME_NO_MALLOC
    Eigen::internal::set_is_malloc_allowed(false);
#endif
    for (int k = 0; k < 1000; ++k) {
        loop.follow(1.49 + 0.001 * k, qs[k % 10], qds[k % 10], &terms, c);
        violations += c.violation == Violation::None ? 0 : 1;
        sink += c.tau[0] + SafetyStop::torque(qs[k % 10], qds[k % 10], qs[0])[1];
        loopF.follow(1.49 + 0.001 * k, qs[k % 10], qds[k % 10], &terms, c);
        violations += c.violation == Violation::None ? 0 : 1;
        sink += c.tauF[4];
        sink += stop.step(1.49 + 0.001 * k, qs[k % 10], qds[k % 10], stopDone)[2];
        sink += loop.observe(qs[k % 10], qds[k % 10], st) ? st.xi[0] : 1.0;
        violations += monitor.checkOutputs(xi, 0.001, 0.05) == Violation::None ? 0 : 1;
        violations += monitor.checkState(qs[k % 10], qds[k % 10]) == Violation::None ? 0 : 1;
        violations += monitor.checkAcceleration(c.qdd) == Violation::None ? 0 : 1;
        violations += monitor.checkTorque(c.tauG, c.tau, c.tau, 1e-3) == Violation::None ? 0 : 1;
    }
#ifdef EIGEN_RUNTIME_NO_MALLOC
    Eigen::internal::set_is_malloc_allowed(true);
#endif
    const long after = gNewCalls.load(std::memory_order_relaxed);
    CAPTURE(before, after, sink, violations);
    REQUIRE(after == before);
}
