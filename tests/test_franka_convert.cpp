// Tests of the robot driver's conversions between libfranka's arrays and se3pf's types (apps/robot/franka_convert.hpp):
// libfranka's matrices are column-major, and its 6 x 7 body Jacobian has the rows (v, omega) where se3pf has (omega,
// v). The libfranka-style arrays here are built from first principles (finite differences of the flange pose), not
// from PandaArm's Jacobian, so a wrong swap or a transposed layout fails.
#include <catch2/catch_test_macros.hpp>
#include "franka_convert.hpp"
#include "se3pf/panda.hpp"
#include "se3pf/se3.hpp"
#include "test_helpers.hpp"

#include <array>
#include <cstddef>

using namespace se3pf;
using namespace se3pf::robot;
using testing::maxDiff;

TEST_CASE("franka_convert: column-major poses, 7 x 7 and 3 x 3 matrices, vectors")
{
    const Pose g = PandaArm::output(testing::randomConfig(), nullptr);
    std::array<double, 16> a{};
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
            a[static_cast<std::size_t>(c * 4 + r)] = g(r, c);
    CHECK(maxDiff(poseFromColumnMajor(a), g) == 0.0);
    // libfranka's O_T_EE layout: the translation is in elements 12, 13, 14
    CHECK(a[12] == g(0, 3));
    CHECK(a[14] == g(2, 3));

    Mat<7, 7> M = testing::randn<Mat<7, 7>>();
    std::array<double, 49> m{};
    for (int r = 0; r < 7; ++r)
        for (int c = 0; c < 7; ++c)
            m[static_cast<std::size_t>(c * 7 + r)] = M(r, c);
    CHECK(maxDiff(mat7FromColumnMajor(m), M) == 0.0);

    const Matrix3d I = testing::randn<Matrix3d>();
    std::array<double, 9> i9{};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            i9[static_cast<std::size_t>(c * 3 + r)] = I(r, c);
    CHECK(maxDiff(mat3FromColumnMajor(i9), I) == 0.0);

    const Vec7 v = testing::randn<Vec7>();
    CHECK(maxDiff(vec7(array7(v)), v) == 0.0);
    CHECK(maxDiff(vec3({{1.0, 2.0, 3.0}}), Vector3d(1.0, 2.0, 3.0)) == 0.0);
}

TEST_CASE("franka_convert: a libfranka body Jacobian (rows v, omega; column-major) becomes PandaArm's J_b")
{
    // libfranka's bodyJacobian(Frame::kFlange) of joint c: the flange twist for qd = e_c in the flange frame,
    // translational part first: v = R^T dp/dq_c, omega = (R^T dR/dq_c)^vee. Built by central differences.
    for (int trial = 0; trial < 5; ++trial) {
        const Vec7 q = testing::randomConfig(0.1);
        std::array<double, 42> franka{};
        const double h = 1e-6;
        const Pose g = PandaArm::output(q, nullptr);
        for (int c = 0; c < 7; ++c) {
            Vec7 dq = Vec7::Zero();
            dq[c] = h;
            const Pose gp = PandaArm::output(Vec7(q + dq), nullptr), gm = PandaArm::output(Vec7(q - dq), nullptr);
            const Vector3d v = rotation(g).transpose() * (translation(gp) - translation(gm)) / (2.0 * h);
            const Vector3d w = vee(Matrix3d(rotation(g).transpose() * (rotation(gp) - rotation(gm)) / (2.0 * h)));
            for (int r = 0; r < 3; ++r) {
                franka[static_cast<std::size_t>(c * 6 + r)] = v[r];
                franka[static_cast<std::size_t>(c * 6 + 3 + r)] = w[r];
            }
        }
        Mat67 J = Mat67::Zero();
        PandaArm::output(q, &J);
        CHECK(maxDiff(bodyJacobianFromFranka(franka), J) < 1e-8);
    }
}
