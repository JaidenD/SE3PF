#include <catch2/catch_test_macros.hpp>
#include "se3pf/metric.hpp"
#include "se3pf/se3.hpp"
#include "test_helpers.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

using namespace Eigen;
using namespace se3pf;
using testing::maxDiff;
using testing::randn;
using testing::uniform;

TEST_CASE("Metric::product builds diag(wR I, wp I); inner product, norm and flat map")
{
    const Metric W = Metric::product(0.01, 1.0);
    Mat6 M = Mat6::Zero();
    M.diagonal() << 0.01, 0.01, 0.01, 1.0, 1.0, 1.0;
    CHECK(maxDiff(W.matrix(), M) == 0.0);

    Vec6 v;
    v << 1, 0, 0, 0, 0, 2;
    CHECK(std::abs(W(v, v) - 4.01) < 1e-15);  // 0.01 * 1 + 1 * 4
    for (int i = 0; i < 10; ++i) {
        const Vec6 a = randn<Vec6>(), b = randn<Vec6>();
        CHECK(std::abs(W(a, b) - a.dot(M * b)) < 1e-14);
        CHECK(std::abs(W(a, b) - W(b, a)) < 1e-14);
        CHECK(std::abs(W.norm(a) - std::sqrt(a.dot(M * a))) < 1e-14);
        CHECK(maxDiff(W(a), Vec6(M * a)) < 1e-15);
        // the named forms and the shorthands are the same functions
        CHECK(W.inner(a, b) == W(a, b));
        CHECK(W.flat(a) == W(a));
    }
}

TEST_CASE("Metric accepts only the product type W = diag(wR I, wp I), wR, wp > 0 (Example 1)")
{
    // Exp, Log and Gamma are the closed forms of Example 1: any other W would silently give the wrong geodesics.
    Mat6 P = Mat6::Zero();
    P.diagonal() << 0.01, 0.01, 0.01, 1.0, 1.0, 1.0;
    CHECK_NOTHROW(Metric(P));
    CHECK(Metric(P).matrix() == P);
    CHECK_NOTHROW(Metric(Mat6::Identity()));
    CHECK_NOTHROW(Metric::product(0.1 * 0.1, 1.0));  // l = 0.1 m, as PandaExample builds it

    Mat6 W = Mat6::Identity();
    W(2, 2) = 2.0;  // diagonal, but not of product type: the speed of the example paths stays constant in it
    CHECK_THROWS_AS(Metric(W), std::invalid_argument);
    W = Mat6::Identity();
    W(0, 3) = W(3, 0) = 0.3;  // couples rotation and translation
    CHECK_THROWS_AS(Metric(W), std::invalid_argument);
    W = P;
    W(4, 4) += 1e-9;  // a relative difference of 1e-9 is not rounding
    CHECK_THROWS_AS(Metric(W), std::invalid_argument);
    W = P;
    W(4, 4) += 1e-15;  // this is
    CHECK_NOTHROW(Metric(W));
    CHECK_THROWS_AS(Metric::product(0.0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(Metric::product(1.0, -1.0), std::invalid_argument);
    CHECK_THROWS_AS(Metric::product(std::numeric_limits<double>::quiet_NaN(), 1.0), std::invalid_argument);
    CHECK_THROWS_AS(Metric::product(std::numeric_limits<double>::infinity(), 1.0), std::invalid_argument);
}

TEST_CASE("Exp_I of Example 1: rotation exp(omega), translation exactly v (not the screw motion)")
{
    const Metric W = Metric::product(0.01, 1.0);
    const Vec6 z = randn<Vec6>();
    const Pose g = W.Exp(z);
    CHECK(maxDiff(Matrix3d(g.block<3,3>(0,0)), expSO3(Vector3d(z.head<3>()))) < 1e-15);
    CHECK(maxDiff(Vector3d(g.block<3,1>(0,3)), Vector3d(z.tail<3>())) == 0.0);
}

TEST_CASE("Log_I inverts Exp_I below the injectivity radius, and fails at a rotation by pi")
{
    const Metric W = Metric::product(0.01, 1.0);
    for (int i = 0; i < 30; ++i) {
        Vec6 z = randn<Vec6>();
        z.head<3>() *= uniform(1e-9, 3.0) / z.head<3>().norm();  // rotation angle in (0, 3) < pi
        Vec6 l;
        REQUIRE(W.Log(W.Exp(z), l));
        CHECK(maxDiff(l, z) < 1e-12);
    }
    Vec6 l;
    const Pose flip = makePose(Matrix3d(Vector3d(1.0, -1.0, -1.0).asDiagonal()), Vector3d(0.1, 0.2, 0.3));
    CHECK_FALSE(W.Log(flip, l));  // rotation by pi about x: trace = -1
}

TEST_CASE("Gamma is the derivative of Exp_I (eq:Gamma); GammaInv inverts it")
{
    const Metric W = Metric::product(0.01, 1.0);
    const double h = 1e-6;
    for (int i = 0; i < 20; ++i) {
        const Vec6 z = randn<Vec6>(0.8), d = randn<Vec6>();
        // body-coordinate derivative of Exp_I along d: (Exp(z)^-1 (Exp(z + h d) - Exp(z - h d)) / 2h)^vee
        const Matrix4d dE = (W.Exp(Vec6(z + h * d)) - W.Exp(Vec6(z - h * d))) / (2 * h);
        const Vec6 fd = vee(Matrix4d(inverse(W.Exp(z)) * dE));
        CHECK(maxDiff(fd, Vec6(W.Gamma(z) * d)) < 1e-8);
        CHECK(maxDiff(Mat6(W.GammaInv(z) * W.Gamma(z)), Mat6::Identity()) < 1e-12);
    }
}

TEST_CASE("distance: dist(I, Exp(z)) = |z|_W, left-invariant, symmetric")
{
    const Metric W = Metric::product(0.01, 1.0);
    for (int i = 0; i < 20; ++i) {
        const Vec6 z = randn<Vec6>(0.5);
        CHECK(std::abs(W.distance(Pose::Identity(), W.Exp(z)) - W.norm(z)) < 1e-12);

        const Pose a = W.Exp(randn<Vec6>(0.5)), b = W.Exp(randn<Vec6>(0.5)), g = W.Exp(randn<Vec6>());
        CHECK(std::abs(W.distance(Pose(g * a), Pose(g * b)) - W.distance(a, b)) < 1e-12);
        CHECK(std::abs(W.distance(a, b) - W.distance(b, a)) < 1e-12);
    }
}
