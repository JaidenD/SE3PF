#include <catch2/catch_test_macros.hpp>
#include "se3pf/metric.hpp"
#include "se3pf/rigid_body.hpp"
#include "se3pf/se3.hpp"
#include "test_helpers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

using namespace Eigen;
using namespace se3pf;
using testing::maxDiff;
using testing::randn;
using testing::randomPose;
using testing::sci;
using testing::twist;

namespace {
    Matrix3d diag3(double a, double b, double c)
    {
        return Matrix3d(Vector3d(a, b, c).asDiagonal());
    }

    // Fixtures for these tests.

    // A fixed start pose away from the identity, so that frame mix-ups (R v versus v) show.
    Pose startPose()
    {
        return makePose(expSO3(Vector3d(0.3, -0.5, 0.8)), Vector3d(0.1, 0.2, -0.3));
    }

    // A fixed rotation Q from principal axes to a tilted body frame. In the tilted frame the inertia Q D Q^T is a
    // full matrix, which catches code that only works on principal axes (that uses just the diagonal of J).
    Matrix3d tilt()
    {
        return expSO3(Vector3d(0.4, -0.2, 0.7));
    }

    // The paper's example body, D = diag(1.0, 1.2, 0.8) and m = 2, described in the tilted body frame.
    RigidBody tiltedExample()
    {
        const Matrix3d Q = tilt();
        return RigidBody(Matrix3d(Q * diag3(1.0, 1.2, 0.8) * Q.transpose()), 2.0);
    }

    // Free motion (u = 0) from startPose() over the given number of 0.5 ms steps: the energy and the spatial
    // angular and linear momenta (R J omega, m R v) stay constant, and R stays in SO(3). The energy alone cannot
    // see a sign error in C zeta, or C computed from the diagonal of a full J; the angular momentum can.
    void checkFreeMotion(const RigidBody& body, int steps)
    {
        Vec6 zeta;
        zeta << 1.0, -0.7, 0.5, 0.4, 0.1, -0.3;
        Pose g = startPose();
        const double h = 0.5e-3;

        const Mat6 M = body.mass();
        const double E0 = 0.5 * zeta.dot(M * zeta);
        const Vec6 mu = M * zeta;
        const Vector3d angular0 = rotation(g) * mu.head<3>();
        const Vector3d linear0 = rotation(g) * mu.tail<3>();

        double energyDrift = 0.0, orthoErr = 0.0, momentumErr = 0.0;
        for (int k = 0; k < steps; ++k) {
            body.integrate(g, zeta, Vec6::Zero(), h);
            const Matrix3d R = rotation(g);
            const Vec6 muk = M * zeta;
            energyDrift = std::max(energyDrift, std::abs(0.5 * zeta.dot(muk) - E0) / E0);
            orthoErr = std::max(orthoErr, maxDiff(Matrix3d(R.transpose() * R), Matrix3d::Identity()));
            momentumErr = std::max(momentumErr, maxDiff(Vector3d(R * muk.head<3>()), angular0));
            momentumErr = std::max(momentumErr, maxDiff(Vector3d(R * muk.tail<3>()), linear0));
        }
        CAPTURE(sci(energyDrift), sci(orthoErr), sci(momentumErr));
        CHECK(energyDrift < 1e-10);  // the reference reaches 6.5e-15
        CHECK(orthoErr < 1e-12);
        CHECK(momentumErr < 1e-10);
    }
}

TEST_CASE("RigidBody::example and mass(): M = diag(J, m I)")
{
    const RigidBody body = RigidBody::example();
    Mat6 M = Mat6::Zero();
    M.diagonal() << 1.0, 1.2, 0.8, 2.0, 2.0, 2.0;
    CHECK(maxDiff(body.mass(), M) == 0.0);

    // a full (non-diagonal) inertia goes into the top-left block unchanged
    const Matrix3d Q = tilt();
    const Matrix3d J = Q * diag3(1.0, 2.0, 3.0) * Q.transpose();
    REQUIRE(maxDiff(J, Matrix3d(J.transpose())) == 0.0);  // so that its symmetric part is J itself
    const RigidBody tilted(J, 0.5);
    CHECK(maxDiff(Matrix3d(tilted.mass().block<3,3>(0,0)), J) == 0.0);
    CHECK(maxDiff(Matrix3d(tilted.mass().block<3,3>(3,3)), Matrix3d(0.5 * Matrix3d::Identity())) == 0.0);
    CHECK(tilted.mass().block<3,3>(0,3).isZero(0.0));
    CHECK(tilted.mass().block<3,3>(3,0).isZero(0.0));
}

TEST_CASE("RigidBody rejects an inertia that is not symmetric positive definite and a mass that is not positive")
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    Matrix3d asym = diag3(1.0, 1.2, 0.8);
    asym(0, 1) = 0.1;

    CHECK_THROWS_AS(RigidBody(asym, 2.0), std::invalid_argument);                  // not symmetric
    CHECK_THROWS_AS(RigidBody(diag3(1.0, -1.0, 1.0), 2.0), std::invalid_argument);  // indefinite
    CHECK_THROWS_AS(RigidBody(diag3(1.0, 1.0, 0.0), 2.0), std::invalid_argument);   // singular
    CHECK_THROWS_AS(RigidBody(Matrix3d::Zero(), 2.0), std::invalid_argument);
    CHECK_THROWS_AS(RigidBody(diag3(1.0, nan, 1.0), 2.0), std::invalid_argument);
    CHECK_THROWS_AS(RigidBody(diag3(1.0, inf, 1.0), 2.0), std::invalid_argument);
    CHECK_THROWS_AS(RigidBody(diag3(1.0, 1.2, 0.8), 0.0), std::invalid_argument);
    CHECK_THROWS_AS(RigidBody(diag3(1.0, 1.2, 0.8), -2.0), std::invalid_argument);
    CHECK_THROWS_AS(RigidBody(diag3(1.0, 1.2, 0.8), nan), std::invalid_argument);
    CHECK_THROWS_AS(RigidBody(diag3(1.0, 1.2, 0.8), inf), std::invalid_argument);
    CHECK_NOTHROW(RigidBody(diag3(1.0, 1.2, 0.8), 2.0));
}

TEST_CASE("RigidBody accepts an inertia that is symmetric up to rounding, and uses its exactly symmetric part")
{
    // A computed inertia differs from its transpose in the last digits: Q D Q^T here is already off by one ulp
    // (5.6e-17), and a J read from another program can be off by more.
    const Matrix3d Q = tilt();
    const Matrix3d J0 = Q * diag3(1.0, 1.2, 0.8) * Q.transpose();
    Matrix3d J = J0;
    J(0, 1) += 1e-15;
    REQUIRE(J(0, 1) != J(1, 0));
    CHECK_NOTHROW(RigidBody(J, 2.0));
    CHECK_NOTHROW(RigidBody(Matrix3d(1e7 * J), 2.0));  // the tolerance is relative: J in g cm^2 instead of kg m^2

    // mass(), coriolis() and acceleration() must describe one matrix: the stored inertia is (J + J^T) / 2, and
    // acceleration() inverts that matrix, not the one in J's lower triangle. An asymmetry near the tolerance
    // (relative 1e-12) is large enough for the second to show.
    J(0, 1) = J0(0, 1) + 5e-13;
    const RigidBody body(J, 2.0);
    const Matrix3d stored = body.mass().block<3,3>(0,0);
    CHECK(maxDiff(stored, Matrix3d(stored.transpose())) == 0.0);
    CHECK(maxDiff(stored, Matrix3d(0.5 * (J + J.transpose()))) == 0.0);
    for (int i = 0; i < 20; ++i) {
        const Vec6 z = randn<Vec6>(), u = randn<Vec6>();
        const double err = maxDiff(Vec6(body.mass() * body.acceleration(z, u) + body.coriolis(z)), u);
        CAPTURE(sci(err));
        CHECK(err < 1e-14);
    }
}

TEST_CASE("RigidBody::coriolis: C zeta = -(J omega x omega, m v x omega), and it does no work")
{
    // By hand, for J = diag(1, 1, 2), m = 2, omega = (1, 0, 1), v = (0, 1, 0):
    //   J omega x omega = (1, 0, 2) x (1, 0, 1) = (0, 1, 0),   m v x omega = 2 (0, 1, 0) x (1, 0, 1) = (2, 0, -2).
    const RigidBody body(diag3(1.0, 1.0, 2.0), 2.0);
    Vec6 expected;
    expected << 0.0, -1.0, 0.0, -2.0, 0.0, 2.0;
    const Vec6 c = body.coriolis(twist(Vector3d(1.0, 0.0, 1.0), Vector3d(0.0, 1.0, 0.0)));
    CAPTURE(c.transpose());
    CHECK(maxDiff(c, expected) < 1e-15);

    // gyroscopic terms are orthogonal to the velocity: zeta^T C zeta = 0, so free motion keeps its energy
    const RigidBody bodies[] = {RigidBody::example(), tiltedExample()};
    for (const RigidBody& b : bodies) {
        for (int i = 0; i < 20; ++i) {
            const Vec6 z = randn<Vec6>();
            const double power = std::abs(z.dot(b.coriolis(z)));
            CAPTURE(sci(power));
            CHECK(power < 1e-13);
        }
    }

    // A full J: C is covariant under a change of body frame Q. With D = diag(1.0, 1.2, 0.8),
    //     C_{Q D Q^T}(Q omega, Q v) = diag(Q, Q) C_D(omega, v),   since (Q D omega) x (Q omega) = Q (D omega x omega).
    // Code that uses only the diagonal of J (the principal-axis form of Euler's equations) breaks this. The checks
    // above cannot see that: zeta^T C zeta = 0 still holds, and the acceleration test compares with coriolis() itself.
    const Matrix3d Q = tilt();
    const RigidBody principal = RigidBody::example();
    const RigidBody tilted = tiltedExample();
    for (int i = 0; i < 20; ++i) {
        const Vector3d w = randn<Vector3d>(), v = randn<Vector3d>();
        const Vec6 cp = principal.coriolis(twist(w, v));
        const Vec6 expected = twist(Vector3d(Q * cp.head<3>()), Vector3d(Q * cp.tail<3>()));
        const double err = maxDiff(tilted.coriolis(twist(Vector3d(Q * w), Vector3d(Q * v))), expected);
        CAPTURE(sci(err));
        CHECK(err < 1e-13);
    }
}

TEST_CASE("RigidBody::acceleration solves M zeta' + C zeta = u (eq:EL with G = 0)")
{
    // the example, and a full J with a different mass
    const Matrix3d Q = tilt();
    const RigidBody bodies[] = {RigidBody::example(),
                                RigidBody(Matrix3d(Q * diag3(1.0, 2.0, 3.0) * Q.transpose()), 0.5)};
    for (const RigidBody& body : bodies) {
        for (int i = 0; i < 20; ++i) {
            const Vec6 z = randn<Vec6>(), u = randn<Vec6>();
            const Vec6 a = body.acceleration(z, u);
            const double err = maxDiff(Vec6(body.mass() * a + body.coriolis(z)), u);
            CAPTURE(sci(err));
            CHECK(err < 1e-13);
        }
    }
}

TEST_CASE("RigidBody as a System: Pi = id and J_b = I")
{
    STATIC_REQUIRE(RigidBody::kDof == 6);
    const Pose g = randomPose();
    Mat6 Jb = Mat6::Zero();
    CHECK(maxDiff(RigidBody::output(g, &Jb), g) == 0.0);
    CHECK(maxDiff(Jb, Mat6::Identity()) == 0.0);
    CHECK(maxDiff(RigidBody::output(g, nullptr), g) == 0.0);  // Jb is optional
}

TEST_CASE("RigidBody::retract: a curve through g with velocity zeta, equal to g Exp_I(h zeta)")
{
    const Metric W = Metric::product(0.01, 1.0);  // the product-type Exp_I does not depend on w_R, w_p
    const double h = 1e-6;
    for (int i = 0; i < 20; ++i) {
        const Pose g = randomPose();
        const Vec6 z = randn<Vec6>();

        const double zeroErr = maxDiff(RigidBody::retract(g, z, 0.0), g);
        CAPTURE(sci(zeroErr));
        CHECK(zeroErr < 1e-15);

        // derivative at h = 0 in body coordinates: (g^-1 (c(h) - c(-h)) / 2h)^vee = zeta
        const Matrix4d dc = (RigidBody::retract(g, z, h) - RigidBody::retract(g, z, -h)) / (2 * h);
        const double err = maxDiff(vee(Matrix4d(inverse(g) * dc)), z);
        CAPTURE(sci(err));
        CHECK(err < 1e-8);

        const double t = testing::uniform(-1.0, 1.0);
        const double expErr = maxDiff(RigidBody::retract(g, z, t), Pose(g * W.Exp(Vec6(t * z))));
        CAPTURE(sci(expErr));
        CHECK(expErr < 1e-13);
    }
}

TEST_CASE("RigidBody::integrate, free motion: energy and momenta conserved, R stays in SO(3)")
{
    SECTION("principal axes, J = diag(1.0, 1.2, 0.8), 10 s")
    {
        checkFreeMotion(RigidBody::example(), 20000);
    }

    SECTION("tilted body frame, full J = Q diag(1.0, 1.2, 0.8) Q^T, 1 s")
    {
        checkFreeMotion(tiltedExample(), 2000);  // 1 s is enough here and keeps the Debug run short
    }
}

TEST_CASE("RigidBody::integrate, closed-form motion of a symmetric body (catches a sign error in C zeta)")
{
    // J = diag(1, 1, 2), omega0 = (1, 0, 1): Euler's equations give omega1' = -omega2, omega2' = omega1,
    // omega3' = 0, so omega(t) = (cos t, sin t, 1). A flipped sign in C gives (cos t, -sin t, 1).
    // The attitude is a precession about the fixed spatial angular momentum a = R0 J omega0 = R0 (1, 0, 2) at
    // rate |a| / J_1 combined with a spin about e3 at rate (1 - J_3 / J_1) omega3 = -1:
    //     R(t) = exp(t a^) R0 exp(-t e3^).
    // With no force, the centre of mass moves on a straight line: R v = R0 v0 and p(t) = p0 + t R0 v0.
    const RigidBody body(diag3(1.0, 1.0, 2.0), 2.0);
    const Pose g0 = startPose();
    const Matrix3d R0 = rotation(g0);
    const Vector3d p0 = translation(g0);
    const Vector3d v0(0.3, -0.2, 0.5);
    const Vector3d a = R0 * Vector3d(1.0, 0.0, 2.0);  // J_1 = 1, so a is also the precession rate

    Pose g = g0;
    Vec6 zeta = twist(Vector3d(1.0, 0.0, 1.0), v0);
    const double h = 0.5e-3;
    const int steps = 20000;  // 10 s

    double omegaErr = 0.0, attitudeErr = 0.0, velocityErr = 0.0, positionErr = 0.0;
    for (int k = 1; k <= steps; ++k) {
        body.integrate(g, zeta, Vec6::Zero(), h);
        const double t = k * h;
        const Matrix3d R = rotation(g);
        const Matrix3d Rexact = expSO3(Vector3d(t * a)) * R0 * expSO3(Vector3d(0.0, 0.0, -t));
        omegaErr = std::max(omegaErr, maxDiff(Vector3d(zeta.head<3>()), Vector3d(std::cos(t), std::sin(t), 1.0)));
        attitudeErr = std::max(attitudeErr, maxDiff(R, Rexact));
        velocityErr = std::max(velocityErr, maxDiff(Vector3d(R * zeta.tail<3>()), Vector3d(R0 * v0)));
        positionErr = std::max(positionErr, maxDiff(translation(g), Vector3d(p0 + t * R0 * v0)));
    }

    CAPTURE(sci(omegaErr), sci(attitudeErr), sci(velocityErr), sci(positionErr));
    CHECK(omegaErr < 1e-12);
    CHECK(attitudeErr < 1e-10);
    CHECK(velocityErr < 1e-10);
    CHECK(positionErr < 1e-10);
}

TEST_CASE("RigidBody::integrate with a constant input: force and torque enter with the right sign and scale")
{
    const RigidBody body = RigidBody::example();  // J = diag(1.0, 1.2, 0.8), m = 2
    const Pose g0 = startPose();
    const Matrix3d R0 = rotation(g0);
    const Vector3d p0 = translation(g0);
    const double h = 1e-3, T = 1.0;
    const int steps = 1000;

    SECTION("a body force without rotation: v = v0 + (f / m) t, p = p0 + R0 (v0 t + f t^2 / 2m)")
    {
        const Vector3d v0(0.1, 0.0, -0.2), f(1.0, -2.0, 0.5);
        Pose g = g0;
        Vec6 zeta = twist(Vector3d::Zero(), v0);
        for (int k = 0; k < steps; ++k)
            body.integrate(g, zeta, twist(Vector3d::Zero(), f), h);

        const double attitudeErr = maxDiff(rotation(g), R0);
        const double velocityErr = maxDiff(Vector3d(zeta.tail<3>()), Vector3d(v0 + (T / 2.0) * f));
        const double positionErr = maxDiff(translation(g), Vector3d(p0 + R0 * (T * v0 + (T * T / 4.0) * f)));
        CAPTURE(sci(attitudeErr), sci(velocityErr), sci(positionErr));
        CHECK(attitudeErr < 1e-15);
        CHECK(velocityErr < 1e-12);
        CHECK(positionErr < 1e-12);
    }

    SECTION("a torque about the principal axis e3: omega3 = w0 + (tau / J3) t, R = R0 Rz(w0 t + tau t^2 / 2 J3)")
    {
        const double w0 = 0.5, tau = 0.4;
        Pose g = g0;
        Vec6 zeta = twist(Vector3d(0.0, 0.0, w0), Vector3d::Zero());
        for (int k = 0; k < steps; ++k)
            body.integrate(g, zeta, twist(Vector3d(0.0, 0.0, tau), Vector3d::Zero()), h);

        const double angle = w0 * T + 0.5 * (tau / 0.8) * T * T;
        const double twistErr = maxDiff(zeta, twist(Vector3d(0.0, 0.0, w0 + (tau / 0.8) * T), Vector3d::Zero()));
        const double attitudeErr = maxDiff(rotation(g), Matrix3d(R0 * expSO3(Vector3d(0.0, 0.0, angle))));
        CAPTURE(sci(twistErr), sci(attitudeErr));
        CHECK(twistErr < 1e-12);
        CHECK(attitudeErr < 1e-12);
        CHECK(maxDiff(translation(g), p0) == 0.0);
    }
}
