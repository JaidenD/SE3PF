#include <catch2/catch_test_macros.hpp>
#include "se3pf/se3.hpp"
#include "test_helpers.hpp"

using namespace Eigen;
using namespace se3pf;
using testing::maxDiff;
using testing::randn;
using testing::randomPose;

TEST_CASE("hat: skew-symmetric, hat(w) x = w cross x, vee inverts it")
{
    for (int i = 0; i < 20; ++i) {
        const Vector3d w = randn<Vector3d>(), x = randn<Vector3d>();
        const Matrix3d W = hat(w);
        CHECK(maxDiff(W, Matrix3d(-W.transpose())) < 1e-15);
        CHECK(maxDiff(Vector3d(W * x), Vector3d(w.cross(x))) < 1e-14);
        CHECK(maxDiff(vee(W), w) == 0.0);
    }
}

TEST_CASE("hat on se(3): [hat(omega) v; 0 0], and vee inverts it")
{
    const Vec6 z = randn<Vec6>();
    const Matrix4d X = hat(z);
    CHECK(maxDiff(Vector3d(X.block<3,1>(0,3)), Vector3d(z.tail<3>())) == 0.0);
    CHECK(X.row(3).cwiseAbs().maxCoeff() == 0.0);  // an element of the Lie algebra: bottom row all zero
    CHECK(maxDiff(vee(X), z) == 0.0);
}

TEST_CASE("expSO3: a rotation, agrees with Eigen, turns the right way, exp(0) = I")
{
    for (int i = 0; i < 20; ++i) {
        const Vector3d w = randn<Vector3d>();
        const Matrix3d R = expSO3(w);
        CHECK(maxDiff(Matrix3d(R.transpose() * R), Matrix3d::Identity()) < 1e-14);
        CHECK(std::abs(R.determinant() - 1.0) < 1e-14);
        const Matrix3d Reigen = AngleAxisd(w.norm(), w.normalized()).toRotationMatrix();
        CHECK(maxDiff(R, Reigen) < 1e-14);
    }
    CHECK(maxDiff(expSO3(Vector3d::Zero()), Matrix3d::Identity()) == 0.0);
    // a quarter turn about z takes x to y
    const Vector3d y = expSO3(Vector3d(0.0, 0.0, EIGEN_PI / 2)) * Vector3d::UnitX();
    CHECK(maxDiff(y, Vector3d(Vector3d::UnitY())) < 1e-15);
}

TEST_CASE("logSO3 inverts expSO3, including tiny angles and angles near pi")
{
    const double pi = EIGEN_PI;  // EIGEN_PI is a long double; mixing it into the list below would not compile
    for (int i = 0; i < 20; ++i) {
        const Vector3d dir = randn<Vector3d>().normalized();
        for (double angle : {1e-12, 1e-9, 1e-4, 0.5, 2.0, pi - 1e-3, pi - 1e-6}) {
            const Vector3d w = angle * dir;
            CHECK(maxDiff(logSO3(expSO3(w)), w) < 1e-9);
        }
    }
}

TEST_CASE("Jr is the derivative of expSO3; JrInv inverts it")
{
    const double h = 1e-6;
    for (int i = 0; i < 20; ++i) {
        const Vector3d w = randn<Vector3d>(0.8), d = randn<Vector3d>();
        // d/dt exp(w + t d) at t = 0, in body coordinates: (exp(w)^T (exp(w + h d) - exp(w - h d)) / 2h)^vee
        const Matrix3d dR = (expSO3(Vector3d(w + h * d)) - expSO3(Vector3d(w - h * d))) / (2 * h);
        const Vector3d fd = vee(Matrix3d(expSO3(w).transpose() * dR));
        CHECK(maxDiff(fd, Vector3d(Jr(w) * d)) < 1e-8);
        CHECK(maxDiff(Matrix3d(JrInv(w) * Jr(w)), Matrix3d::Identity()) < 1e-12);
    }
    // the small-angle branches
    const Vector3d tiny(1e-9, -2e-9, 3e-9);
    CHECK(maxDiff(Matrix3d(JrInv(tiny) * Jr(tiny)), Matrix3d::Identity()) < 1e-14);
}

TEST_CASE("makePose puts R and p in the right blocks, rotation and translation read them; inverse is the group inverse")
{
    const Matrix3d R = expSO3(Vector3d(0.3, -0.2, 0.5));
    const Vector3d p(1.0, 2.0, 3.0);
    const Pose g = makePose(R, p);
    CHECK(maxDiff(Matrix3d(g.block<3,3>(0,0)), R) == 0.0);
    CHECK(maxDiff(Vector3d(g.block<3,1>(0,3)), p) == 0.0);
    CHECK(maxDiff(rotation(g), R) == 0.0);
    CHECK(maxDiff(translation(g), p) == 0.0);
    CHECK(maxDiff(Vector4d(g.row(3).transpose()), Vector4d(0, 0, 0, 1)) == 0.0);
    for (int i = 0; i < 10; ++i) {
        const Pose h = randomPose();
        CHECK(maxDiff(Matrix4d(h * inverse(h)), Matrix4d::Identity()) < 1e-14);
        CHECK(maxDiff(inverse(h), Matrix4d(h.inverse())) < 1e-13);
    }
}

TEST_CASE("Ad_g z = (g hat(z) g^-1)^vee; Ad is a homomorphism")
{
    for (int i = 0; i < 20; ++i) {
        const Pose g = randomPose(), h = randomPose();
        const Vec6 z = randn<Vec6>();
        const Vec6 byDefinition = vee(Matrix4d(g * hat(z) * inverse(g)));
        CHECK(maxDiff(Vec6(Ad(g) * z), byDefinition) < 1e-12);
        CHECK(maxDiff(Ad(Pose(g * h)), Mat6(Ad(g) * Ad(h))) < 1e-12);
    }
}

TEST_CASE("ad_a b = (hat(a) hat(b) - hat(b) hat(a))^vee")
{
    for (int i = 0; i < 20; ++i) {
        const Vec6 a = randn<Vec6>(), b = randn<Vec6>();
        const Vec6 commutator = vee(Matrix4d(hat(a) * hat(b) - hat(b) * hat(a)));
        CHECK(maxDiff(Vec6(ad(a) * b), commutator) < 1e-13);
    }
}
