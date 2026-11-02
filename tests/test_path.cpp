#include <catch2/catch_test_macros.hpp>
#include "se3pf/path.hpp"
#include "se3pf/se3.hpp"
#include "test_helpers.hpp"

#include <cmath>
#include <stdexcept>

using namespace Eigen;
using namespace se3pf;
using testing::maxDiff;
using testing::restAttitude;
using testing::uniform;

namespace {
    // Properties every arc-length path must have, checked at many s.
    void checkArcLengthPath(const Path& P, const Metric& W)
    {
        const double h = 1e-6, L = P.length();
        for (int i = 0; i < 50; ++i) {
            const double s = uniform(0.0, L);
            CHECK(std::abs(W.norm(P.velocity(s)) - 1.0) < 1e-12);  // unit speed in the metric

            // zeta_gamma(s) = (gamma(s)^-1 gamma'(s))^vee, gamma' by a central difference
            const Matrix4d dg = (P.pose(s + h) - P.pose(s - h)) / (2 * h);
            CHECK(maxDiff(vee(Matrix4d(inverse(P.pose(s)) * dg)), P.velocity(s)) < 1e-8);

            // zeta_gamma'(s) by a central difference of zeta_gamma
            const Vec6 dz = (P.velocity(s + h) - P.velocity(s - h)) / (2 * h);
            CHECK(maxDiff(dz, P.acceleration(s)) < 1e-6);
        }
        // closed and periodic
        CHECK(maxDiff(P.pose(L), P.pose(0.0)) < 1e-12);
        CHECK(maxDiff(P.velocity(L), P.velocity(0.0)) < 1e-12);
        CHECK(maxDiff(P.pose(-0.1), P.pose(L - 0.1)) < 1e-12);
    }
}

TEST_CASE("coning circle (Panda example, l = 0.1 m): reference values")
{
    const Metric W = Metric::product(0.01, 1.0);
    const ConingCircle path(restAttitude(), Vector3d(0.40, 0.0, 0.45), 0.12, 0.25, W);
    CHECK(std::abs(path.length() - 0.77008765) < 1e-8);

    const Pose g0 = path.pose(0.0);
    CHECK(maxDiff(Vector3d(g0.block<3,1>(0,3)), Vector3d(0.52, 0.0, 0.45)) < 1e-15);

    Vec6 z0;  // zeta_gamma(0) of the reference implementation
    z0 << 0.0, 2.01858182142, -0.253645174458, -0.692318537798, -0.670796031053, 0.171282347316;
    CHECK(maxDiff(path.velocity(0.0), z0) < 1e-10);
}

TEST_CASE("coning circle: unit speed, consistent derivatives, closed")
{
    const Metric W = Metric::product(0.01, 1.0);
    const ConingCircle path(restAttitude(), Vector3d(0.40, 0.0, 0.45), 0.12, 0.25, W);
    checkArcLengthPath(path, W);

    const Metric W1 = Metric::product(1.0, 1.0);  // l = 1 m
    const ConingCircle path1(restAttitude(), Vector3d(0.40, 0.0, 0.45), 0.12, 0.25, W1);
    CHECK(std::abs(path1.length() - 1.73869664) < 1e-8);
    checkArcLengthPath(path1, W1);
}

TEST_CASE("tumbling circle (rigid-body example, W = I): L = 2 sqrt(3) pi, unit speed, closed")
{
    const Metric W = Metric::product(1.0, 1.0);
    const TumblingCircle path(1.0, W);
    CHECK(std::abs(path.length() - 2.0 * std::sqrt(3.0) * EIGEN_PI) < 1e-12);

    Vec6 z0;
    z0 << 1, 0, 1, 0, 1, 0;
    CHECK(maxDiff(path.velocity(0.0), Vec6(z0 / std::sqrt(3.0))) < 1e-15);
    checkArcLengthPath(path, W);
}

TEST_CASE("paths work through a const Path& (virtual functions)")
{
    const Metric W = Metric::product(1.0, 1.0);
    const TumblingCircle tumbling(1.0, W);
    const Path& p = tumbling;
    CHECK(std::abs(p.length() - tumbling.length()) == 0.0);
    CHECK(maxDiff(p.pose(1.0), tumbling.pose(1.0)) == 0.0);
}

TEST_CASE("a metric in which the speed would not be constant is rejected (by Metric itself)")
{
    // The paths have constant speed in every product-type W, the only kind Metric accepts; a W that couples
    // rotation and translation does not get as far as the path's own constant-speed check.
    Mat6 M = Mat6::Identity();
    M(0, 3) = M(3, 0) = 0.3;  // couples rotation and translation: not of product type
    CHECK_THROWS_AS(Metric(M), std::invalid_argument);
}

// ---------------------------------------------------------------------------------------------------------------
// SampledPath and the paper paths (2026-11-02)

namespace {
    // As checkArcLengthPath, with tolerances for the numerical derivatives of SampledPath (larger step h)
    void checkSampledPath(const SampledPath& P, const Metric& W)
    {
        const double h = 1e-4, L = P.length();
        for (int i = 0; i < 50; ++i) {
            const double s = uniform(0.0, L);
            CHECK(std::abs(W.norm(P.velocity(s)) - 1.0) < 1e-12);
            Vec6 zp, zm;  // a smaller step for the pose difference: its truncation error grows with the curvature
            REQUIRE(W.Log(inverse(P.pose(s)) * P.pose(s + 0.1 * h), zp));
            REQUIRE(W.Log(inverse(P.pose(s)) * P.pose(s - 0.1 * h), zm));
            CHECK(maxDiff(Vec6((zp - zm) / (0.2 * h)), P.velocity(s)) < 1e-6);
            const Vec6 dz = (P.velocity(s + h) - P.velocity(s - h)) / (2 * h);
            CHECK(maxDiff(dz, P.acceleration(s)) < 1e-4 * std::max(1.0, P.acceleration(s).norm()));
        }
        CHECK(maxDiff(P.pose(L), P.pose(0.0)) < 1e-12);
        CHECK(maxDiff(P.pose(-0.1), P.pose(L - 0.1)) < 1e-9);
        CHECK(maxDiff(P.velocity(L), P.velocity(0.0)) < 1e-9);
    }
}

TEST_CASE("SampledPath: the coning circle as a curve gives the coning circle")
{
    const Metric W = Metric::product(0.01, 1.0);
    const ConingCircle c(restAttitude(), Vector3d(0.40, 0.0, 0.45), 0.12, 0.25, W);
    const SampledPath sp([&c](double th) { return c.pose(th / (2.0 * EIGEN_PI) * c.length()); }, W);
    CHECK(std::abs(sp.length() - c.length()) < 1e-10);
    for (int i = 0; i < 50; ++i) {
        const double s = uniform(0.0, c.length());
        CHECK(maxDiff(sp.pose(s), c.pose(s)) < 1e-10);
        CHECK(maxDiff(sp.velocity(s), c.velocity(s)) < 1e-8);
        CHECK(maxDiff(sp.acceleration(s), c.acceleration(s)) < 1e-6);
    }
    checkSampledPath(sp, W);
}

TEST_CASE("SampledPath: arc length is the integral of the W-speed of the curve")
{
    const Metric W = Metric::product(0.01, 1.0);
    const SampledPath p(inspectionOrbit(restAttitude(), InspectionOrbit()), W);
    // theta(s(theta)) = theta, s by the trapezoid rule on a fine grid
    double s = 0.0, thPrev = 0.0, vPrev = W.norm(p.curveVelocity(0.0));
    for (int k = 1; k <= 20000; ++k) {
        const double th = 2.0 * EIGEN_PI * k / 20000;
        const double v = W.norm(p.curveVelocity(th));
        s += 0.5 * (v + vPrev) * (th - thPrev);
        thPrev = th;
        vPrev = v;
        if (k % 2000 == 0 && k < 20000)
            CHECK(std::abs(p.theta(s) - th) < 1e-6);
    }
    CHECK(std::abs(s - p.length()) < 1e-6);
}

TEST_CASE("SampledPath: rejects too few nodes, no curve, or a curve that stops")
{
    const Metric W = Metric::product(0.01, 1.0);
    CHECK_THROWS_AS(SampledPath(inspectionOrbit(restAttitude(), InspectionOrbit()), W, 8), std::invalid_argument);
    CHECK_THROWS_AS(SampledPath(SampledPath::Curve(), W), std::invalid_argument);
    const Pose g = Pose::Identity();
    CHECK_THROWS_AS(SampledPath([g](double) { return g; }, W), std::invalid_argument);
}

TEST_CASE("paper paths: the inspection orbit and the dome figure-eight")
{
    const Metric W = Metric::product(0.01, 1.0);
    const Matrix3d R0 = restAttitude();

    const InspectionOrbit o;
    const SampledPath orbit(inspectionOrbit(R0, o), W);
    CHECK(std::abs(orbit.length() - 0.78868) < 1e-4);  // the experiments' L (src/experiments.cpp)
    checkSampledPath(orbit, W);
    for (int i = 0; i < 20; ++i) {  // the tool z-axis points at the target
        const Pose g = orbit.pose(uniform(0.0, orbit.length()));
        const Vector3d toTarget = (o.center - Vector3d(0.0, 0.0, o.depth) - g.block<3, 1>(0, 3)).normalized();
        CHECK(maxDiff(Vector3d(g.block<3, 1>(0, 2)), toTarget) < 1e-12);
    }

    const DomeFigureEight d;
    const SampledPath dome(domeFigureEight(R0, d), W);
    CHECK(std::abs(dome.length() - 0.98243) < 1e-4);
    checkSampledPath(dome, W);
    for (int i = 0; i < 20; ++i) {  // on the dome, the tool z-axis along the inward normal
        const Pose g = dome.pose(uniform(0.0, dome.length()));
        const Vector3d r = g.block<3, 1>(0, 3) - d.center;
        CHECK(std::abs(r.norm() - d.radius) < 1e-12);
        CHECK(maxDiff(Vector3d(g.block<3, 1>(0, 2)), Vector3d(-r.normalized())) < 1e-12);
    }
    // The two passes over the top have the same position but attitudes 2 psi0 apart: the path does not cross itself
    const Pose top1 = domeFigureEight(R0, d)(0.0), top2 = domeFigureEight(R0, d)(EIGEN_PI);
    CHECK(maxDiff(Vector3d(top1.block<3, 1>(0, 3)), Vector3d(top2.block<3, 1>(0, 3))) < 1e-12);
    CHECK(std::abs(logSO3(Matrix3d(top1.block<3, 3>(0, 0).transpose() * top2.block<3, 3>(0, 0))).norm()
                   - 2.0 * d.twist) < 1e-12);
}

TEST_CASE("paper paths: the straight rail is nearly straight around s0, at constant attitude")
{
    const Metric W = Metric::product(0.01, 1.0);
    const StraightRail r;
    const SampledPath rail(straightRail(restAttitude(), r), W);
    checkSampledPath(rail, W);
    // s0 = 0.1 is the middle of the near side
    const Pose g0 = rail.pose(0.1);
    CHECK(std::abs(g0(1, 3) - r.center.y()) < 2e-3);
    CHECK(std::abs(g0(0, 3) - (r.center.x() - r.halfWidth)) < 1e-9);
    // +-10 cm along the rail the side is within 1.5 mm of the line x = c_x - w
    for (double ds : {-0.1, -0.05, 0.05, 0.1}) {
        const Pose g = rail.pose(0.1 + ds);
        CHECK(std::abs(g(0, 3) - (r.center.x() - r.halfWidth)) < 1.5e-3);
        CHECK(maxDiff(Matrix3d(g.block<3, 3>(0, 0)), restAttitude()) < 1e-12);
    }
}
