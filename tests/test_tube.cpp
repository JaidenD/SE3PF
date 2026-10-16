#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "se3pf/panda_example.hpp"
#include "se3pf/rigid_body_example.hpp"
#include "se3pf/se3.hpp"
#include "se3pf/tube.hpp"
#include "test_helpers.hpp"
#include "tube_fixture.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace Eigen;
using namespace se3pf;
using Catch::Matchers::WithinAbs;
using testing::arcDiff;
using testing::maxDiff;
using testing::pandaTube;
using testing::pandaTube1;
using testing::poseError;
using testing::randn;
using testing::rigidTube;
using testing::sci;
using testing::uniform;

namespace {
    // Start offsets of the two examples (Appendix)
    Vec5 pandaXi0() { return PandaExample::xi0(); }
    Vec5 rigidXi0() { return RigidBodyExample::xi0(); }

    // Random transverse offset with |xi| = r
    Vec5 randomXi(double r)
    {
        const Vec5 x = randn<Vec5>();
        return (r / x.norm()) * x;
    }

    // dist(g, C) by brute force: the smallest distance to the samples of the path
    double sampledDistance(const Tube& T, const std::vector<Pose>& samples, const Pose& g)
    {
        double best = std::numeric_limits<double>::infinity();
        for (const Pose& gp : samples)
            best = std::min(best, T.metric().distance(gp, g));
        return best;
    }

    // Rotation by th in the (i, j) coordinate plane of R^5: e_i -> cos e_i + sin e_j
    Mat5 planeRotation(int i, int j, double th)
    {
        Mat5 R = Mat5::Identity();
        R(i, i) = R(j, j) = std::cos(th);
        R(j, i) = std::sin(th);
        R(i, j) = -std::sin(th);
        return R;
    }
}

// ---------------------------------------------------------------------------------------------------------------
// Rotations in SO(5) by invariant planes (the logarithm of the holonomy)

TEST_CASE("rotationPlanes: logarithm and exponential of a rotation in SO(5)")
{
    const HouseholderQR<Mat5> qr(randn<Mat5>());
    const Mat5 Q = qr.householderQ();  // random orthogonal matrix

    // One plane with angle 2.5, three fixed directions
    const Mat5 H = Q * planeRotation(0, 1, 2.5) * Q.transpose();
    const RotationPlanes P = rotationPlanes(H);
    REQUIRE(P.count == 1);
    CHECK(std::abs(P.theta[0] - 2.5) < 1e-12);
    const Mat5 Phi = P.log();
    const double errSkew = maxDiff(Phi, Mat5(-Phi.transpose()));
    const double errExp = maxDiff(P.exp(1.0), H);
    const double errGroup = maxDiff(Mat5(P.exp(0.3) * P.exp(0.7)), P.exp(1.0));
    const double errInverse = maxDiff(P.exp(-0.4), Mat5(P.exp(0.4).transpose()));
    CAPTURE(sci(errSkew), sci(errExp), sci(errGroup), sci(errInverse));
    CHECK(errSkew < 1e-15);
    CHECK(errExp < 1e-13);
    CHECK(errGroup < 1e-14);  // a one-parameter group
    CHECK(errInverse < 1e-15);

    // Two planes, the second with angle 1e-5 (1 - cos = 5e-11: a plane, not a fixed direction). Its eigenvalue
    // cos(1e-5) is only 5e-11 from the fixed direction's 1, so its eigenvectors, and exp(log H), are accurate to
    // about 2 eps / theta = 5e-11 (8e-14 for the tubes' holonomies). The constant in front depends on Q: over 100
    // random Q the error reached 1.03e-10, hence the bound of 1e-9 (a dropped plane would give 1e-5).
    const Mat5 H2 = Q * planeRotation(0, 1, 2.5) * planeRotation(2, 3, 1e-5) * Q.transpose();
    const RotationPlanes P2 = rotationPlanes(H2);
    REQUIRE(P2.count == 2);
    CHECK(std::abs(P2.theta[0] - 2.5) < 1e-12);
    CHECK(std::abs(P2.theta[1] - 1e-5) < 1e-14);
    const double err2 = maxDiff(P2.exp(1.0), H2);
    CAPTURE(sci(err2));
    CHECK(err2 < 1e-9);

    // 1 - cos(theta) < 1e-12 counts as fixed: an angle of 1e-7 is dropped
    const RotationPlanes P1 = rotationPlanes(Mat5(Q * planeRotation(0, 1, 1e-7) * Q.transpose()));
    CHECK(P1.count == 0);
    CHECK(rotationPlanes(Mat5::Identity()).count == 0);
    CHECK(maxDiff(rotationPlanes(Mat5::Identity()).exp(0.5), Mat5::Identity()) == 0.0);

    // not a rotation
    Mat5 F = Mat5::Identity();
    F(4, 4) = -1.0;
    CHECK_THROWS_AS(rotationPlanes(F), std::invalid_argument);
}

TEST_CASE("rotationPlanes: a rotation by pi has a logarithm too")
{
    // H u = -u on the whole plane, so H does not tell which v belongs to u. Any partner from the (-1)-eigenspace
    // gives a logarithm, and any logarithm closes the frame (eq:holonomy). Check that the one returned is one:
    // exp(log H) = H, and exp(log H / 2) is a rotation whose square is H (it fails if v is not a unit vector of
    // the (-1)-eigenspace orthogonal to u).
    const HouseholderQR<Mat5> qr(randn<Mat5>());
    const Mat5 Q = qr.householderQ();
    const Mat5 onePlane = Q * planeRotation(0, 1, EIGEN_PI) * Q.transpose();
    const Mat5 twoPlanes = Q * planeRotation(0, 1, EIGEN_PI) * planeRotation(2, 3, EIGEN_PI) * Q.transpose();
    const Mat5 mixed = Q * planeRotation(0, 1, EIGEN_PI) * planeRotation(2, 3, 0.7) * Q.transpose();

    struct Case {
        Mat5 H;
        int planes;
    };
    for (const Case& c : {Case{onePlane, 1}, Case{twoPlanes, 2}, Case{mixed, 2}}) {
        const RotationPlanes P = rotationPlanes(c.H);
        REQUIRE(P.count == c.planes);
        CHECK(std::abs(P.theta[0] - EIGEN_PI) < 1e-9);
        const Mat5 Phi = P.log();
        const Mat5 half = P.exp(0.5);
        const double errSkew = maxDiff(Phi, Mat5(-Phi.transpose()));
        const double errExp = maxDiff(P.exp(1.0), c.H);
        const double errOrtho = maxDiff(Mat5(half.transpose() * half), Mat5::Identity());
        const double errSquare = maxDiff(Mat5(half * half), c.H);
        CAPTURE(c.planes, sci(errSkew), sci(errExp), sci(errOrtho), sci(errSquare));
        CHECK(errSkew < 1e-15);
        CHECK(errExp < 1e-12);
        CHECK(errOrtho < 1e-12);
        CHECK(errSquare < 1e-12);
    }
    CHECK(std::abs(rotationPlanes(mixed).theta[1] - 0.7) < 1e-12);

    // Not orthogonal (det = 0.5 > 0): H e_1 = -e_1, but no other direction is flipped, so e_1 has no partner.
    // (A rotation in SO(5) always has an even-dimensional (-1)-eigenspace.)
    Mat5 lone = Mat5::Identity();
    lone(0, 0) = -1.0;
    lone(4, 4) = -0.5;
    CHECK_THROWS_AS(rotationPlanes(lone), std::runtime_error);
}

// ---------------------------------------------------------------------------------------------------------------
// TUBE-1: the normal frame (Lemma 1)

TEST_CASE("Tube: a temporary path does not compile (the tube keeps a reference)")
{
    STATIC_REQUIRE(std::is_constructible_v<Tube, const ConingCircle&, const Metric&>);
    STATIC_REQUIRE(std::is_constructible_v<Tube, const TumblingCircle&, const Metric&, int>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<Tube, ConingCircle, const Metric&>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<Tube, TumblingCircle&&, const Metric&, int>);
}

TEST_CASE("Tube: initial frame by Gram-Schmidt on zeta_gamma(0), e_1, ..., e_6, in that order")
{
    // Rigid body: zeta_gamma(0) = (e_1 + e_3 + e_5)/sqrt(3), so e_5 drops out and N(0) holds the residuals of
    // e_1, e_2, e_3, e_4, e_6. (exp(-(0/L) Phi) = I, so N(0) = N0(0).)
    Mat65 expected = Mat65::Zero();
    expected.col(0) << 2, 0, -1, 0, -1, 0;
    expected.col(0) /= std::sqrt(6.0);
    expected(1, 1) = 1.0;
    expected.col(2) << 0, 0, 1, 0, -1, 0;
    expected.col(2) /= std::sqrt(2.0);
    expected(3, 3) = 1.0;
    expected(5, 4) = 1.0;
    Mat65 N;
    rigidTube().frame(0.0, N);
    const double errRigid = maxDiff(N, expected);
    CAPTURE(sci(errRigid));
    CHECK(errRigid < 1e-12);

    // Panda: zeta_gamma(0) has no omega_x, so the first column is e_1 / |e_1|_W = 10 e_1; e_6 drops out.
    pandaTube().frame(0.0, N);
    const double errPanda = maxDiff(Vec6(N.col(0)), Vec6(10.0 * Vec6::Unit(0)));
    CAPTURE(sci(errPanda));
    CHECK(errPanda < 1e-12);
}

TEST_CASE("Tube: N^T W N = I and N^T W zeta_gamma = 0 (eq:frame-props)")
{
    for (const Tube* T : {&pandaTube(), &pandaTube1(), &rigidTube()}) {
        const Mat6& W = T->metric().matrix();
        for (int i = 0; i < 50; ++i) {
            const double s = uniform(-T->length(), 2.0 * T->length());
            Mat65 N;
            T->frame(s, N);
            const double errOrtho = maxDiff(Mat5(N.transpose() * W * N), Mat5::Identity());
            const double errNormal = Vec5(N.transpose() * W * T->path().velocity(s)).cwiseAbs().maxCoeff();
            CAPTURE(T->length(), s, sci(errOrtho), sci(errNormal));
            CHECK(errOrtho < 1e-9);
            CHECK(errNormal < 1e-9);
        }
    }
}

TEST_CASE("Tube: holonomy H is a rotation; Phi = log H is skew and has two planes")
{
    for (const Tube* T : {&pandaTube(), &pandaTube1(), &rigidTube()}) {
        const TubeDiagnostics& d = T->diagnostics();
        CAPTURE(T->length(), sci(d.detH - 1.0), sci(d.orthogonalityError), sci(d.skewError), sci(d.logError));
        CHECK(std::abs(d.detH - 1.0) < 1e-10);
        CHECK(d.orthogonalityError < 1e-10);
        CHECK(d.skewError < 1e-12);
        // exp(Phi) = H_r: no plane lost (8e-5 for the Panda), no sign flipped. Not to 1e-15: the eigenvalue
        // cos(8e-5) of the Panda's small plane is 3.2e-9 from the fixed direction's 1, which limits the plane's
        // eigenvectors to about 2 eps / theta = 3e-12.
        CHECK(d.logError < 1e-10);
        CHECK(d.planes == 2);
    }
}

TEST_CASE("Tube: rotation angles of the holonomy (Testing anchors)")
{
    // They do not depend on the initial frame (H becomes Q^T H Q) and hardly on the grid.
    const TubeDiagnostics& a = pandaTube().diagnostics();
    CHECK_THAT(a.angles[0], WithinAbs(0.161182, 1e-5));
    CHECK_THAT(a.angles[1], WithinAbs(8.0e-5, 1e-5));
    const TubeDiagnostics& b = pandaTube1().diagnostics();
    CHECK_THAT(b.angles[0], WithinAbs(0.707417, 1e-5));
    CHECK_THAT(b.angles[1], WithinAbs(3.09e-4, 1e-5));
    const TubeDiagnostics& c = rigidTube().diagnostics();
    CHECK_THAT(c.angles[0], WithinAbs(2.166040, 1e-5));
    CHECK_THAT(c.angles[1], WithinAbs(0.492067, 1e-5));
}

TEST_CASE("Tube: the frame closes, N(s_n) = N(s_0) on the node values (eq:holonomy)")
{
    // What eq:holonomy is for. The frame properties above hold for the open frame N0 too, so they cannot see a
    // wrong sign in exp(-(s/L) Phi) or a dropped plane (the Panda's second plane alone leaves a gap of about 1e-4).
    // Measured on the node values: frame() wraps s = L to 0 and returns node 0 exactly.
    for (const Tube* T : {&pandaTube(), &pandaTube1(), &rigidTube()}) {
        const double gap = T->diagnostics().closureGap;
        CAPTURE(T->length(), sci(gap));
        CHECK(gap < 1e-9);
    }
}

TEST_CASE("Tube: N' agrees with central differences of N (eq:frame-ode)")
{
    for (const Tube* T : {&pandaTube(), &pandaTube1(), &rigidTube()}) {
        // At the nodes, against a fourth-order difference of the node values (step = node spacing). This tests the
        // node derivatives from eq:frame-ode themselves: a difference with a step below the node spacing only sees
        // the cubic interpolant, whose derivative is consistent with its values whatever the node derivatives are.
        const int n = T->gridIntervals();
        const double h = T->length() / n;
        for (int i = 0; i < 50; ++i) {
            const int k = (i * 997) % n;  // spread over the grid, including k = 0 (the differences wrap around)
            Mat65 N, dN, a, b, c, d;
            T->frame(k * h, N, &dN);
            T->frame((k - 2) * h, a);
            T->frame((k - 1) * h, b);
            T->frame((k + 1) * h, c);
            T->frame((k + 2) * h, d);
            const Mat65 fd = (a - 8.0 * b + 8.0 * c - d) / (12.0 * h);
            const double err = maxDiff(dN, fd);
            CAPTURE(T->length(), k, sci(err));
            CHECK(err < 1e-6);
        }
        // Between the nodes: the derivative of the interpolant
        for (int i = 0; i < 20; ++i) {
            const double s = uniform(0.0, T->length()), e = 1e-6;
            Mat65 N, dN, Np, Nm;
            T->frame(s, N, &dN);
            T->frame(s + e, Np);
            T->frame(s - e, Nm);
            const double err = maxDiff(dN, Mat65((Np - Nm) / (2.0 * e)));
            CAPTURE(T->length(), s, sci(err));
            CHECK(err < 1e-6);
        }
    }
}

TEST_CASE("Tube: construction rejects too coarse a grid")
{
    const Metric W = Metric::product(1.0, 1.0);
    const TumblingCircle path(1.0, W);
    CHECK_THROWS_AS(Tube(path, W, 2), std::invalid_argument);
}

// ---------------------------------------------------------------------------------------------------------------
// TUBE-2: the tubular map and the outputs (Proposition 1)

TEST_CASE("Tube::wrap maps s into [0, L)")
{
    const Tube& T = pandaTube();
    const double L = T.length();
    CHECK(T.wrap(0.0) == 0.0);
    CHECK(T.wrap(L) == 0.0);
    CHECK(T.wrap(-1e-17) == 0.0);  // -1e-17 + L rounds to L
    CHECK_THAT(T.wrap(-0.1), WithinAbs(L - 0.1, 1e-15));
    CHECK_THAT(T.wrap(3.0 * L + 0.2), WithinAbs(0.2, 1e-14));
    for (int i = 0; i < 100; ++i) {
        const double w = T.wrap(uniform(-10.0, 10.0));
        CHECK(w >= 0.0);
        CHECK(w < L);
    }
}

TEST_CASE("Tube: Psi(s, 0) = gamma(s), and the start poses of the examples")
{
    for (const Tube* T : {&pandaTube(), &rigidTube()}) {
        const double s = uniform(0.0, T->length());
        CHECK(maxDiff(T->Psi(s, Vec5::Zero()), T->path().pose(s)) < 1e-15);
    }

    // These numbers hold only for the Gram-Schmidt order of Lemma 1 step 1 (reversed columns give the rigid body
    // 15.50 deg and 0.563 m instead).
    double deg = 0.0, m = 0.0;
    poseError(rigidTube().path().pose(0.1), rigidTube().Psi(0.1, rigidXi0()), deg, m);
    CHECK_THAT(deg, WithinAbs(25.45, 0.005));
    CHECK_THAT(m, WithinAbs(0.4389, 5e-5));

    poseError(pandaTube().path().pose(0.1), pandaTube().Psi(0.1, pandaXi0()), deg, m);
    CHECK_THAT(deg, WithinAbs(15.07, 0.005));
    CHECK_THAT(m, WithinAbs(0.01682, 5e-6));

    // and the Panda start round-trips: |xi0| = 0.031225
    double eta = -1.0;
    Vec5 xi = Vec5::Zero();
    REQUIRE(pandaTube().outputs(pandaTube().Psi(0.1, pandaXi0()), 0.1, eta, xi));
    CHECK_THAT(eta, WithinAbs(0.1, 1e-9));
    CHECK(maxDiff(xi, pandaXi0()) < 1e-9);
    CHECK_THAT(xi.norm(), WithinAbs(0.031225, 1e-6));
}

TEST_CASE("Tube: Psi^-1 round trip and |xi| = dist(g, C) on the Panda tube (Proposition 1)")
{
    const Tube& T = pandaTube();
    std::vector<Pose> samples(20000);  // for the brute-force distance
    for (int k = 0; k < 20000; ++k)
        samples[k] = T.path().pose(T.length() * k / 20000);

    for (int i = 0; i < 20; ++i) {
        const double s = uniform(0.0, T.length());
        const Vec5 xi = randomXi(uniform(0.005, 0.06));
        const Pose g = T.Psi(s, xi);

        double eta = -1.0;
        Vec5 xr = Vec5::Zero();
        int iterations = 0;
        REQUIRE(T.outputs(g, s + 0.02, eta, xr, &iterations));
        const double errEta = std::abs(arcDiff(eta, s, T.length()));
        const double errXi = maxDiff(xr, xi);
        CAPTURE(s, xi.norm(), eta, iterations, sci(errEta), sci(errXi));
        CHECK(eta >= 0.0);
        CHECK(eta < T.length());
        CHECK(errEta < 1e-9);
        CHECK(errXi < 1e-9);
        CHECK(iterations >= 1);
        CHECK(iterations <= 20);

        // Proposition 1(c): gamma(eta) is the closest point and dist(g, C) = |xi|. For every fourth pose only,
        // since 20000 distances take a while in Debug.
        if (i % 4 == 0) {
            const double errDist = std::abs(sampledDistance(T, samples, g) - xi.norm());
            CAPTURE(sci(errDist));
            CHECK(errDist < 2e-6);
        }
    }
}

TEST_CASE("Tube: round trip on the rigid-body tube, and its start pose")
{
    const Tube& T = rigidTube();
    double eta = -1.0;
    Vec5 xi = Vec5::Zero();
    REQUIRE(T.outputs(T.Psi(0.1, rigidXi0()), 0.1, eta, xi));
    CHECK_THAT(eta, WithinAbs(0.1, 1e-9));
    CHECK(maxDiff(xi, rigidXi0()) < 1e-9);
    CHECK_THAT(xi.norm(), WithinAbs(0.6245, 5e-5));

    for (int i = 0; i < 10; ++i) {
        const double s = uniform(0.0, T.length());
        const Vec5 x = randomXi(uniform(0.05, 0.6));
        REQUIRE(T.outputs(T.Psi(s, x), s - 0.05, eta, xi));
        const double errEta = std::abs(arcDiff(eta, s, T.length()));
        const double errXi = maxDiff(xi, x);
        CAPTURE(s, x.norm(), eta, sci(errEta), sci(errXi));
        CHECK(eta >= 0.0);
        CHECK(eta < T.length());
        CHECK(errEta < 1e-9);
        CHECK(errXi < 1e-9);
    }
}

TEST_CASE("Tube: phi_g' against a central difference of phi_g, and phi_g'(eta) = -1 on the path")
{
    for (const Tube* T : {&pandaTube(), &rigidTube()}) {
        const double scale = T == &pandaTube() ? 0.03 : 0.4;
        for (int i = 0; i < 10; ++i) {
            const double s = uniform(0.0, T->length());
            const Pose g = T->Psi(s, randomXi(scale));
            const double t = s + uniform(-0.05, 0.05), h = 1e-6;  // not at the root
            double phi, dphi, phiP, phiM, unused;
            Vec6 e, eP, eM;
            REQUIRE(T->projectionCondition(g, t, phi, dphi, e));
            REQUIRE(T->projectionCondition(g, t + h, phiP, unused, eP));
            REQUIRE(T->projectionCondition(g, t - h, phiM, unused, eM));
            const double err = std::abs(dphi - (phiP - phiM) / (2.0 * h));
            CAPTURE(T->length(), t, sci(err));
            CHECK(err < 1e-6);
            // e = e_g(t) is the logarithm of gamma(t)^-1 g
            CHECK(maxDiff(T->metric().Exp(e), Pose(inverse(T->path().pose(t)) * g)) < 1e-12);
        }
        // On the path: e_g(s) = 0, phi_g(s) = 0 and phi_g'(s) = -|zeta_gamma|_W^2 = -1
        for (int i = 0; i < 10; ++i) {
            const double s = uniform(0.0, T->length());
            double phi, dphi;
            Vec6 e;
            REQUIRE(T->projectionCondition(T->path().pose(s), s, phi, dphi, e));
            CAPTURE(sci(phi), sci(dphi + 1.0));
            CHECK(std::abs(phi) < 1e-12);
            CHECK(std::abs(dphi + 1.0) < 1e-9);
        }
    }
}

TEST_CASE("Tube: Newton converges from guesses offset by +-0.05; at least one update")
{
    for (const Tube* T : {&pandaTube(), &rigidTube()}) {
        const double scale = T == &pandaTube() ? 0.04 : 0.5;
        for (int i = 0; i < 20; ++i) {
            const double s = uniform(0.0, T->length());
            const Pose g = T->Psi(s, randomXi(uniform(0.0, scale)));
            for (double offset : {-0.05, 0.05}) {
                double eta = -1.0;
                Vec5 xi = Vec5::Zero();
                REQUIRE(T->outputs(g, s + offset, eta, xi));
                const double err = std::abs(arcDiff(eta, s, T->length()));
                CAPTURE(T->length(), s, offset, eta, sci(err));
                CHECK(eta >= 0.0);
                CHECK(eta < T->length());
                CHECK(err < 1e-9);
            }
        }
    }

    // Started at the root, Newton still makes one update (the count is the number of updates)
    const Tube& T = pandaTube();
    double eta = -1.0;
    Vec5 xi = Vec5::Zero();
    int iterations = 0;
    REQUIRE(T.outputs(T.Psi(0.3, pandaXi0()), 0.3, eta, xi, &iterations));
    CHECK(iterations == 1);
}

TEST_CASE("Tube: Newton steps across the seam s = 0 = L and still returns eta in [0, L)")
{
    // Random s rarely lands this close to the seam. The guess s - 0.05 for s = 0.01 wraps to the far end, near L,
    // and Newton's update must wrap back past L to 0.01 (the other way round for s = L - 0.01). Without the wrap
    // in the update, eta would come out as s + L or s - L.
    for (const Tube* T : {&pandaTube(), &rigidTube()}) {
        const double L = T->length();
        const Vec5 xi0 = T == &pandaTube() ? pandaXi0() : rigidXi0();
        for (const double s : {0.01, L - 0.01}) {
            const double guess = s < 0.5 * L ? s - 0.05 : s + 0.05;
            double eta = -1.0;
            Vec5 xi = Vec5::Zero();
            REQUIRE(T->outputs(T->Psi(s, xi0), guess, eta, xi));
            const double errXi = maxDiff(xi, xi0);
            CAPTURE(L, s, guess, eta, sci(errXi));
            CHECK(eta >= 0.0);
            CHECK(eta < L);
            CHECK_THAT(eta, WithinAbs(s, 1e-9));
            CHECK(errXi < 1e-9);
        }

        // globalProject: the closest of the 720 samples to a pose at s = L - 2e-4 is s = 0 (the sample spacing is
        // L / 720 > 1e-3), so Newton steps back across the seam.
        const double s = L - 2e-4;
        double eta = -1.0;
        REQUIRE(T->globalProject(T->Psi(s, xi0), eta));
        CAPTURE(L, s, eta);
        CHECK(eta >= 0.0);
        CHECK(eta < L);
        CHECK_THAT(eta, WithinAbs(s, 1e-9));
    }
}

TEST_CASE("Tube: globalProject finds the closest point where Newton from L/2 does not")
{
    struct Case {
        const Tube* tube;
        Vec5 xi0;
        double wrongRoot;
    };
    // The wrong roots of phi_g that Newton reaches from L/2 (Online: projection, pitfall)
    for (const Case& c : {Case{&rigidTube(), rigidXi0(), 5.614}, Case{&pandaTube(), pandaXi0(), 0.485}}) {
        const Tube& T = *c.tube;
        const Pose g = T.Psi(0.1, c.xi0);

        double eta = -1.0;
        REQUIRE(T.globalProject(g, eta));
        CHECK_THAT(eta, WithinAbs(0.1, 1e-9));

        Vec5 xi = Vec5::Zero();
        REQUIRE(T.outputs(g, 0.5 * T.length(), eta, xi));
        CHECK_THAT(eta, WithinAbs(c.wrongRoot, 1e-3));
        CHECK(xi.norm() > c.xi0.norm());  // a root of phi_g, but not the closest point
    }

    // random poses in the Panda tube
    const Tube& T = pandaTube();
    for (int i = 0; i < 10; ++i) {
        const double s = uniform(0.0, T.length());
        double eta = -1.0;
        REQUIRE(T.globalProject(T.Psi(s, randomXi(uniform(0.0, 0.1))), eta));
        CAPTURE(s, eta);
        CHECK(eta >= 0.0);
        CHECK(eta < T.length());
        CHECK(std::abs(arcDiff(eta, s, T.length())) < 1e-9);
    }
}

TEST_CASE("Tube: Newton fails at a stationary point of phi_g (|phi'| < 1e-12)")
{
    // For the Panda start pose, phi_g has the right root 0.1 (phi' = -0.993) and the wrong root 0.485 (phi' =
    // +0.993) found from L/2, so phi' vanishes at some t* in between. Find t* by bisection on the sign of phi'.
    // From t*, the Newton step phi/phi' is about 1e14: outputs must report failure instead of landing anywhere.
    const Tube& T = pandaTube();
    const Pose g = T.Psi(0.1, pandaXi0());
    double phi = 0.0, dphi = 0.0;
    Vec6 e = Vec6::Zero();
    double a = 0.1, b = 0.484983;
    REQUIRE(T.projectionCondition(g, a, phi, dphi, e));
    REQUIRE(dphi < 0.0);
    REQUIRE(T.projectionCondition(g, b, phi, dphi, e));
    REQUIRE(dphi > 0.0);
    for (int k = 0; k < 60; ++k) {  // past the resolution of a double: a and b end up adjacent
        const double mid = 0.5 * (a + b);
        REQUIRE(T.projectionCondition(g, mid, phi, dphi, e));
        if (dphi < 0.0)
            a = mid;
        else
            b = mid;
    }
    const double tStar = a;
    REQUIRE(T.projectionCondition(g, tStar, phi, dphi, e));
    CAPTURE(tStar, phi, sci(dphi));
    CHECK(std::abs(dphi) < 1e-12);

    double eta = -1.0;
    Vec5 xi = Vec5::Zero();
    CHECK_FALSE(T.outputs(g, tStar, eta, xi));
    CHECK(eta == -1.0);  // outputs untouched on failure
    CHECK(maxDiff(xi, Vec5::Zero()) == 0.0);
}

TEST_CASE("Tube: outputs and globalProject report failure as false, not as an exception")
{
    const Tube& T = pandaTube();
    Pose g = T.path().pose(0.2);
    g(0, 3) = std::numeric_limits<double>::quiet_NaN();
    double eta = -1.0;
    Vec5 xi = Vec5::Zero();
    double phi, dphi;
    Vec6 e;
    CHECK_FALSE(T.projectionCondition(g, 0.2, phi, dphi, e));
    CHECK_FALSE(T.outputs(g, 0.2, eta, xi));
    CHECK_FALSE(T.globalProject(g, eta));
    CHECK(eta == -1.0);  // outputs untouched on failure
}

TEST_CASE("Tube: a non-finite s gives NaN or false, not a crash")
{
    // The spline under frame() would compute a node index from NaN / h (undefined behaviour, then an
    // out-of-range read). The online functions report it as a value instead.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const Tube* T : {&pandaTube(), &pandaTube1(), &rigidTube()}) {
        for (const double s : {nan, inf, -inf}) {
            CAPTURE(T->length(), s);
            Mat65 N, dN;
            T->frame(s, N, &dN);
            CHECK_FALSE(N.allFinite());
            CHECK_FALSE(dN.allFinite());
            CHECK_FALSE(T->Psi(s, Vec5::Zero()).allFinite());
            CHECK_FALSE(T->JPsi(s, Vec5::Zero()).allFinite());

            const Pose g = T->path().pose(0.1);
            double phi, dphi;
            Vec6 e;
            CHECK_FALSE(T->projectionCondition(g, s, phi, dphi, e));
            double eta = -1.0;
            Vec5 xi = Vec5::Zero();
            CHECK_FALSE(T->outputs(g, s, eta, xi));  // a non-finite guess
            CHECK(eta == -1.0);
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------
// TUBE-3: the body Jacobian J_Psi

TEST_CASE("Tube: J_Psi against central differences of Psi")
{
    for (const Tube* T : {&pandaTube(), &pandaTube1(), &rigidTube()}) {
        const double scale = T == &rigidTube() ? 0.3 : 0.015;
        for (int i = 0; i < 10; ++i) {
            const double eta = uniform(0.0, T->length()), h = 1e-6;
            const Vec5 xi = randomXi(scale);
            const Vec6 d = randn<Vec6>();  // direction (d eta, d xi)
            const Vec5 dxi = d.tail<5>();
            const Pose gP = T->Psi(eta + h * d[0], xi + h * dxi);
            const Pose gM = T->Psi(eta - h * d[0], xi - h * dxi);
            const Vec6 fd = vee(Matrix4d(inverse(T->Psi(eta, xi)) * (gP - gM) / (2.0 * h)));
            const double err = maxDiff(fd, Vec6(T->JPsi(eta, xi) * d));
            CAPTURE(T->length(), eta, sci(err));
            CHECK(err < 1e-6);
        }
    }
}

TEST_CASE("Tube: J_Psi(eta, 0) = [zeta_gamma(eta) | N(eta)]")
{
    for (const Tube* T : {&pandaTube(), &pandaTube1(), &rigidTube()}) {
        for (int i = 0; i < 10; ++i) {
            const double eta = uniform(0.0, T->length());
            Mat65 N;
            T->frame(eta, N);
            Mat6 expected;
            expected << T->path().velocity(eta), N;
            const double err = maxDiff(T->JPsi(eta, Vec5::Zero()), expected);
            CAPTURE(T->length(), eta, sci(err));
            CHECK(err < 1e-14);
        }
    }
}

TEST_CASE("Tube: time of one outputs() + JPsi() call (Panda)", "[.timing]")
{
    // Release build only: ./build-rel/tests "[.timing]"
#ifndef NDEBUG
    SKIP("timing is meaningful in a Release build only");
#endif
    const Tube& T = pandaTube();
    const Pose g = T.Psi(0.1, pandaXi0());
    const int calls = 20000;
    double sum = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < calls; ++i) {
        double eta = 0.0;
        Vec5 xi = Vec5::Zero();
        T.outputs(g, 0.1 + 1e-4 * (i % 2), eta, xi);  // a prediction off by up to 1e-4: one Newton update
        sum += eta + T.JPsi(eta, xi)(0, 0);
    }
    const double us = 1e6 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / calls;
    WARN("outputs() + JPsi(): " << us << " us per call (checksum " << sum << ")");
    CHECK(us < 20.0);  // a few microseconds in Release; a Debug build is about 200 times slower
}

// ---------------------------------------------------------------------------------------------------------------
// TUBE-4: the tube-radius estimate

TEST_CASE("estimateTubeRadius: coarse options give the reference value; bad options throw")
{
    // A cheap version of the test, so the default run exercises the code: 8 values of s, the ten directions +-e_i
    // (no random ones, so the result does not depend on the standard library's normal distribution), 100 path
    // samples, 8 bisections. The reference implementation gives 0.124658 and 2.122145 with these options.
    RadiusOptions opt;
    opt.footPoints = 8;
    opt.randomDirections = 0;
    opt.checkPoints = 100;
    opt.bisections = 8;
    CHECK_THAT(estimateTubeRadius(pandaTube(), opt), WithinAbs(0.124658, 1e-6));
    CHECK_THAT(estimateTubeRadius(rigidTube(), opt), WithinAbs(2.122145, 1e-6));

    opt.checkPoints = 0;
    CHECK_THROWS_AS(estimateTubeRadius(pandaTube(), opt), std::invalid_argument);
}

// The paper's numbers, so they run in the Release stage of tools/ci.sh (a few seconds); hidden in Debug (3 min there).
// To see the estimates: ./build-rel/tests "estimateTubeRadius: the paper*" -s
TEST_CASE("estimateTubeRadius: the paper's examples (under 1 s each in Release)", SLOW_IN_DEBUG)
{
    struct Case {
        const Tube* tube;
        double expected;   // Testing anchors (the paper reports 1.89, 0.125 and 0.279)
        double startXi;    // |xi0| of the example's start, which must be inside the tube
    };
    for (const Case& c : {Case{&rigidTube(), 1.93, rigidXi0().norm()}, Case{&pandaTube(), 0.125, pandaXi0().norm()},
                          Case{&pandaTube1(), 0.307, pandaXi0().norm()}}) {
        const auto t0 = std::chrono::steady_clock::now();
        const double eps = estimateTubeRadius(*c.tube);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        CAPTURE(eps, c.expected, seconds);
        CHECK(std::abs(eps - c.expected) < 0.1 * c.expected);
        CHECK(c.startXi < eps);
#ifdef NDEBUG
        CHECK(seconds < 5.0);  // the requirement; a Debug build is 30-100 times slower
#endif
    }
}
