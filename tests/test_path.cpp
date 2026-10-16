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
