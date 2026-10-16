#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "se3pf/lift.hpp"
#include "se3pf/panda_example.hpp"
#include "se3pf/se3.hpp"
#include "panda_fixture.hpp"
#include "test_helpers.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace Eigen;
using namespace se3pf;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using Catch::Matchers::WithinAbs;
using testing::maxDiff;
using testing::randn;
using testing::sci;
using testing::uniform;

namespace {
    using Mat7 = Mat<7, 7>;

    // The paper's Panda example (l = 0.1 m), built once for the whole test program (panda_fixture.hpp): its path,
    // tube, rest-posture lift Lift(path, q_r) on 2000 intervals, kernel frame and LiftRedundancy are the objects
    // under test here.
    const Metric& metric() { return testing::pandaExample().metric(); }
    const ConingCircle& path() { return testing::pandaExample().path(); }
    const Tube& tube() { return testing::pandaExample().tube(); }
    const Lift& lift() { return testing::pandaExample().lift(); }
    const KernelFrame& kernelFrame() { return testing::pandaExample().kernelFrame(); }  // X (not the tube's frame N)
    const LiftRedundancy& beta() { return testing::pandaExample().redundancy(); }

    // Small helpers
    Mat67 bodyJacobian(const Vec7& q)
    {
        Mat67 J = Mat67::Zero();
        PandaArm::output(q, &J);
        return J;
    }

    Vec7 kernelAt(const Vec7& q) { return kernelVector(bodyJacobian(q)); }

    Vec7 liftAt(double s)
    {
        Vec7 q = Vec7::Zero();
        lift().eval(s, q);
        return q;
    }

    // Pose error (log R, p) of g^-1 target: zero iff the poses agree. Written out here rather than taken from
    // lift.cpp, so the tests do not share code with what they check.
    double poseGap(const Pose& g, const Pose& target)
    {
        const Pose rel = inverse(g) * target;
        Vec6 e;
        e << logSO3(Matrix3d(rel.block<3,3>(0,0))), rel.block<3,1>(0,3);
        return e.cwiseAbs().maxCoeff();
    }

    // Jacobian of x(q) = kernelVector(J_b(q)) by central differences (step 1e-6), column j = dx/dq_j
    Mat7 kernelJacobian(const Vec7& q)
    {
        const double d = 1e-6;
        Mat7 D = Mat7::Zero();
        for (int j = 0; j < 7; ++j)
            D.col(j) = (kernelAt(Vec7(q + d * Vec7::Unit(j))) - kernelAt(Vec7(q - d * Vec7::Unit(j)))) / (2.0 * d);
        return D;
    }

    // q*'(s) from the implicit function theorem, independent of the lift's central differences: q* solves
    // FK(q) = gamma(s) and h(q) = x(q)^T (q - q_r) = 0 (the Lagrange condition), so [J_b; grad h^T] q*' =
    // (zeta_gamma, 0) with the full gradient grad h = Dx^T (q - q_r) + x.
    Vec7 exactLiftDerivative(double s, const Vec7& q)
    {
        const Vec7 x = kernelAt(q);
        const Vec7 gradH = kernelJacobian(q).transpose() * (q - lift().restPosture()) + x;
        Mat7 A;
        A << bodyJacobian(q), gradH.transpose();
        Vec7 rhs;
        rhs << path().velocity(s), 0.0;
        return A.partialPivLu().solve(rhs);
    }

    // eta(q), xi(q) and d eta(q) = row 0 of J_Psi^-1 J_b, as the controller computes them (eq:A)
    struct Outputs {
        double eta = 0.0;
        Vec5 xi = Vec5::Zero();
        Vec7 deta = Vec7::Zero();
    };

    Outputs outputsAt(const Vec7& q, double guess)
    {
        Mat67 Jb = Mat67::Zero();
        const Pose g = PandaArm::output(q, &Jb);
        Outputs o;
        REQUIRE(tube().outputs(g, guess, o.eta, o.xi));
        const Mat67 B = tube().JPsi(o.eta, o.xi).partialPivLu().solve(Jb);
        o.deta = B.row(0).transpose();
        return o;
    }

    // beta(q) with eta = eta(Pi(q)), the function whose differential dbeta is (eq:beta-flat)
    double betaAt(const Vec7& q, double guess)
    {
        const Outputs o = outputsAt(q, guess);
        double b = 0.0;
        Vec7 db = Vec7::Zero();
        beta().evaluate(q, o.eta, o.deta, b, db);
        return b;
    }

    // A configuration off the lift: pose Psi(s0, xi0) and beta = beta0
    Vec7 startAt(double s0, const Vec5& xi0, double beta0)
    {
        Vec7 q = Vec7::Zero();
        REQUIRE(startConfiguration(tube(), beta(), s0, xi0, beta0, liftAt(s0), q));
        return q;
    }

    Vec5 paperXi0() { return PandaExample::xi0(); }

    // A curve that does not close: the paper's coning circle cut off at 90 % of its length, so gamma(L) != gamma(0)
    // and the continuation ends away from q*(0), as it would if it had changed branches. Not a valid Path (gamma
    // must be periodic): it exists to reach the Lift's closure check. Keeps a copy of the circle, not a reference.
    class OpenArc final : public Path {
        public:
            explicit OpenArc(const ConingCircle& circle) : circle_(circle) {}
            double length() const override { return 0.9 * circle_.length(); }
            Pose pose(double s) const override { return circle_.pose(s); }
            Vec6 velocity(double s) const override { return circle_.velocity(s); }
            Vec6 acceleration(double s) const override { return circle_.acceleration(s); }

        private:
            ConingCircle circle_;
    };
}

// ---------------------------------------------------------------------------------------------------------------
// kernelVector: the generalized cross product

TEST_CASE("kernelVector: the orientation of the cross product on hand-made matrices")
{
    // In R^3, e_1 x e_2 = e_3 and e_2 x e_3 = e_1. The same pattern in R^7: rows e_1..e_6 give e_7, rows e_2..e_7
    // give e_1 (n_0 = (-1)^0 det(I) = +1).
    Mat67 J = Mat67::Zero();
    J.leftCols<6>().setIdentity();
    CHECK(maxDiff(kernelVector(J), Vec7(Vec7::Unit(6))) == 0.0);
    J.setZero();
    J.rightCols<6>().setIdentity();
    CHECK(maxDiff(kernelVector(J), Vec7(Vec7::Unit(0))) == 0.0);

    // J = [I | v]: n_6 = det(I) = 1, and n_j = (-1)^j (-1)^(5-j) v_j = -v_j (moving column v from position 5 to j
    // takes 5 - j swaps), so x = (-v, 1) / |(-v, 1)|.
    Vec6 v;
    v << 0.3, -1.2, 0.5, 2.0, -0.7, 0.1;
    J.leftCols<6>().setIdentity();
    J.col(6) = v;
    Vec7 n;
    n << -v, 1.0;
    const double err = maxDiff(kernelVector(J), Vec7(n.normalized()));
    CAPTURE(sci(err));
    CHECK(err < 1e-15);
}

TEST_CASE("kernelVector: a unit vector of ker J_b with the orientation of the generalized cross product")
{
    // w^T n = det([w^T; J_b]) for every w (expand the determinant along its first row). So n_j = det([e_j^T; J_b]),
    // computed here by 7 x 7 determinants instead of the 6 x 6 minors of kernelVector, and |n| = det([x^T; J_b]) > 0.
    for (int i = 0; i < 20; ++i) {
        const Mat67 J = (i < 10) ? randn<Mat67>() : bodyJacobian(liftAt(uniform(0.0, lift().length())));
        const Vec7 x = kernelVector(J);
        Mat7 A;
        A << x.transpose(), J;
        const double norm = A.determinant();
        Vec7 n = Vec7::Zero();
        for (int j = 0; j < 7; ++j) {
            A << Vec7::Unit(j).transpose(), J;
            n[j] = A.determinant();
        }
        const double errKernel = Vec6(J * x).cwiseAbs().maxCoeff();
        const double errUnit = std::abs(x.norm() - 1.0);
        const double errCross = maxDiff(Vec7(norm * x), n);
        CAPTURE(i, sci(errKernel), sci(errUnit), sci(errCross), norm);
        CHECK(errKernel < 1e-14);
        CHECK(errUnit < 1e-15);
        CHECK(norm > 0.0);
        CHECK(errCross < 1e-13 * std::max(1.0, n.norm()));
    }

    // rank 5: n = 0, no direction to return
    Mat67 J = randn<Mat67>();
    J.row(5) = J.row(4);
    CHECK_FALSE(kernelVector(J).allFinite());
}

// ---------------------------------------------------------------------------------------------------------------
// constrainedIK: Newton's method on [J_b; p^T]

TEST_CASE("constrainedIK: reaches the pose and the hyperplane or fails without touching q")
{
    // From the rest posture onto gamma(0.3) with the coupling q_3 = q_1 (about 0.3 rad away)
    const Pose target = path().pose(0.3);
    const Vec7 p = Vec7::Unit(2) - Vec7::Unit(0);
    Vec7 q = PandaArm::restPosture();
    int iterations = -1;
    REQUIRE(constrainedIK(target, p, 0.0, q, &iterations));
    const double errPose = poseGap(PandaArm::output(q, nullptr), target);
    CAPTURE(sci(errPose), iterations);
    CHECK(errPose < 1e-12);
    CHECK(std::abs(q[2] - q[0]) < 1e-12);
    CHECK(iterations >= 1);
    CHECK(iterations <= 10);  // quadratic convergence (the 0.2 rad cap adds a few updates far away)

    // Already a solution: no update
    int again = -1;
    REQUIRE(constrainedIK(target, p, 0.0, q, &again));
    CHECK(again == 0);

    // Out of reach (2 m away): false, and q is left as it was
    Pose far = target;
    far(0, 3) += 2.0;
    const Vec7 q0 = PandaArm::restPosture();
    Vec7 q1 = q0;
    CHECK_FALSE(constrainedIK(far, p, 0.0, q1));
    CHECK(maxDiff(q1, q0) == 0.0);
}

// ---------------------------------------------------------------------------------------------------------------
// The lift q* (Assumption 1)

TEST_CASE("Lift: LiftRedundancy refuses temporaries (it keeps references)")
{
    STATIC_REQUIRE(std::is_constructible_v<LiftRedundancy, const Lift&, const KernelFrame&>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<LiftRedundancy, Lift&&, const KernelFrame&>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<LiftRedundancy, const Lift&, KernelFrame&&>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<LiftRedundancy, Lift&&, KernelFrame&&>);
    // The lift and the kernel frame copy what they need, so they may be built from temporaries.
    STATIC_REQUIRE(std::is_constructible_v<Lift, ConingCircle&&>);
    STATIC_REQUIRE(std::is_constructible_v<KernelFrame, Lift&&>);
}

TEST_CASE("Lift: FK(q*(s)) = gamma(s) between and at the nodes")
{
    // 40 random s, then the nodes s_k, k = 0, 200, ..., 1800
    const double h = lift().length() / lift().gridIntervals();
    double worst = 0.0;
    for (int i = 0; i < 50; ++i) {
        const double s = (i < 40) ? uniform(0.0, lift().length()) : (i - 40) * 200 * h;
        const double err = poseGap(PandaArm::output(liftAt(s), nullptr), path().pose(s));
        worst = std::max(worst, err);
    }
    CAPTURE(sci(worst));
    CHECK(worst < 1e-8);
}

TEST_CASE("Lift: the branch closes on the node values")
{
    // Measured by the constructor before node n is replaced by node 0: eval wraps s = L to 0, so it cannot see the
    // gap. Both Newton tolerances are 1e-12, so a closed branch leaves a gap at that level (the reference's square
    // solver: 5.9e-13; the projected step here: about 1.5e-13).
    const double gap = lift().diagnostics().closureGap;
    CAPTURE(sci(gap));
    CHECK(gap < 1e-6);   // Assumption 1 (the constructor throws above it)
    CHECK(gap < 1e-11);  // what the stopping rule gives
    CHECK(gap > 0.0);    // exactly 0 would mean it was measured after node n was replaced
    CHECK(maxDiff(lift().nodes().back(), lift().nodes().front()) == 0.0);
    CHECK(static_cast<int>(lift().nodes().size()) == lift().gridIntervals() + 1);
    CHECK(lift().gridIntervals() == 2000);
}

TEST_CASE("Lift: the paper's lift table")
{
    const LiftDiagnostics& d = lift().diagnostics();

    Vec7 expected;
    expected << 0.01457, -0.13186, -0.02614, -2.21552, 0.22458, 2.24779, 0.64710;
    const double errStart = maxDiff(liftAt(0.0), expected);
    CAPTURE(sci(errStart));
    CHECK(errStart < 1e-5);

    // Smallest joint margin: joint 4 (zero-based 3) at its lower limit, near s = 0.367
    CHECK_THAT(d.minJointMargin, WithinAbs(0.3592, 1e-4));
    CHECK(d.minMarginJoint == 3);
    CHECK_THAT(d.minMarginS, WithinAbs(0.367, 1e-3));
    CHECK_THAT(d.minSigma, WithinAbs(0.2070, 1e-4));
    CHECK_THAT(d.maxDerivative, WithinAbs(5.28, 5e-3));

    // Projected Newton steps to 1e-12 (linear convergence): 11 from q_r at s = 0, 9 typical, 13 at most here and in
    // the design doc. Ranges, not equalities: a node right at the 1e-12 tolerance can take one update more or less
    // under other floating-point settings (FMA contraction, another compiler). They still catch a slower solver (a
    // half null-space step: 31 typical) and a continuation that restarts every node from q_r (12 typical, 17 max).
    CAPTURE(d.firstUpdates, d.typicalUpdates, d.maxUpdates);
    CHECK(d.firstUpdates >= 10);
    CHECK(d.firstUpdates <= 12);
    CHECK(d.typicalUpdates >= 8);
    CHECK(d.typicalUpdates <= 10);
    CHECK(d.maxUpdates >= d.typicalUpdates);
    CHECK(d.maxUpdates <= 14);
}

TEST_CASE("Lift: q* is a minimum of |q - q_r|^2 on each self-motion and not a maximum")
{
    // Second difference along the self-motion at every node, with q+- at x^T q = x^T q* +- 0.05: the reference
    // gives 4.1028e-3 at least, at s = 0.206383. Moving along the straight line q* +- 0.05 x instead would give
    // 2 (0.05)^2 = 5e-3 at every node, whatever q* is: the upper check catches that.
    double where = -1.0;
    const double d2 = minimumCheck(lift(), 0.05, &where);
    CAPTURE(sci(d2), where);
    CHECK(d2 >= 4.1e-3);
    CHECK(d2 < 4.2e-3);
    CHECK_THAT(where, WithinAbs(0.2064, 1e-3));

    // Fails closed: with delta = 100 rad, constrainedIK cannot reach the hyperplane (50 updates of at most 0.2 rad
    // move x^T q by at most 10), so q+ is not found and the result is NaN, reported at the first node. A check
    // d2 > 0 then fails.
    double whereFar = -1.0;
    CHECK(std::isnan(minimumCheck(lift(), 100.0, &whereFar)));
    CHECK(whereFar == 0.0);
}

TEST_CASE("Lift: q*' at the nodes agrees with the implicit-function derivative")
{
    // The node derivatives are fourth-order central differences: they agree with the exact q*' to about 1.4e-9
    // (largest over all nodes). A second-order stencil would be off by up to 1.9e-5.
    const int n = lift().gridIntervals();
    const double h = lift().length() / n;
    std::vector<int> nodes = {0, 1, n - 1, 168};  // 168: the largest |q*'|
    for (int i = 0; i < 40; ++i)
        nodes.push_back((i * 157) % n);
    double worst = 0.0;
    for (int k : nodes) {
        Vec7 q = Vec7::Zero(), dq = Vec7::Zero();
        lift().eval(k * h, q, &dq);
        worst = std::max(worst, maxDiff(dq, exactLiftDerivative(k * h, q)));
    }
    CAPTURE(sci(worst));
    CHECK(worst < 5e-9);
}

TEST_CASE("Lift: q*' between the nodes is the derivative of the interpolant; eval wraps and rejects NaN")
{
    for (int i = 0; i < 20; ++i) {
        const double s = uniform(0.0, lift().length()), e = 1e-6;
        Vec7 q = Vec7::Zero(), dq = Vec7::Zero();
        lift().eval(s, q, &dq);
        const double err = maxDiff(dq, Vec7((liftAt(s + e) - liftAt(s - e)) / (2.0 * e)));
        CAPTURE(s, sci(err));
        CHECK(err < 1e-6);
        CHECK(maxDiff(liftAt(s + lift().length()), q) < 1e-12);
        CHECK(maxDiff(liftAt(s - 3.0 * lift().length()), q) < 1e-12);
    }
    Vec7 q = Vec7::Zero(), dq = Vec7::Zero();
    lift().eval(std::nan(""), q, &dq);
    CHECK_FALSE(q.allFinite());
    CHECK_FALSE(dq.allFinite());
}

TEST_CASE("Lift: construction fails loudly")
{
    // Each failure by its own message, so a test cannot pass on the wrong throw. (The KernelFrame's rank-loss throw
    // is not tested: no Lift can reach it, because a Lift that was built has a Cholesky factor of J_b J_b^T at every
    // node, so J_b(q*) has full rank there and kernelVector is finite. It stays as a safeguard.)

    // Too coarse a grid for the central differences, caught before any Newton step (the spline would throw the same
    // type later, so the message tells them apart)
    CHECK_THROWS_MATCHES(Lift(path(), PandaArm::restPosture(), 3), std::invalid_argument,
                         MessageMatches(ContainsSubstring("Lift: need at least 4 grid intervals")));

    // A path out of reach: Newton's method fails at the first node
    const ConingCircle far(Matrix3d::Identity(), Vector3d(2.0, 0.0, 0.45), 0.12, 0.25, metric());
    CHECK_THROWS_MATCHES(Lift(far), std::runtime_error, MessageMatches(ContainsSubstring("Newton")));

    // Reachable, but only outside the joint limits: the paper's path turned 2.3 rad about the flange's z-axis
    // (the joint-7 axis) needs q_7 = 3.23 > 2.8973 at s = 0 (0.65 on the paper's path).
    const Matrix3d R0 = PandaArm::output(PandaArm::restPosture(), nullptr).block<3,3>(0,0);
    const ConingCircle twisted(Matrix3d(R0 * expSO3(Vector3d(0.0, 0.0, 2.3))), Vector3d(0.40, 0.0, 0.45), 0.12,
                               0.25, metric());
    CHECK_THROWS_MATCHES(Lift(twisted), std::runtime_error,
                         MessageMatches(ContainsSubstring("joint 6 (zero-based) left its range")));

    // A curve that does not close: q*(L) != q*(0). (A coarse grid keeps it fast; the continuation still converges.)
    const OpenArc open(path());
    CHECK_THROWS_MATCHES(Lift(open, PandaArm::restPosture(), 200), std::runtime_error,
                         MessageMatches(ContainsSubstring("does not close")));
    // The same coarse grid on the closed path builds: the throw above is the closure check, not the grid.
    CHECK_NOTHROW(Lift(path(), PandaArm::restPosture(), 200));
}

// ---------------------------------------------------------------------------------------------------------------
// The kernel frame X (Lemma 3, r = 1)

TEST_CASE("KernelFrame: X spans ker J_b(q*) between the nodes")
{
    // At s_k + ds/4, ds/2, 3ds/4: the spline X against x(q*(s)) = kernelVector(J_b(q*(s))) with q*(s) from the
    // lift's spline. Both interpolate smooth functions with h = 3.9e-4, so they agree to about 3e-13.
    const int n = lift().gridIntervals();
    const double h = lift().length() / n;
    double errX = 0.0, errKernel = 0.0, errUnit = 0.0, errNodes = 0.0;
    for (int k = 0; k < n; ++k) {
        Vec7 X = Vec7::Zero();
        kernelFrame().eval(k * h, X);
        errNodes = std::max(errNodes, maxDiff(X, kernelAt(lift().nodes()[k])));
        for (double f : {0.25, 0.5, 0.75}) {
            const double s = (k + f) * h;
            kernelFrame().eval(s, X);
            const Mat67 J = bodyJacobian(liftAt(s));
            errX = std::max(errX, maxDiff(X, kernelVector(J)));
            errKernel = std::max(errKernel, Vec6(J * X).cwiseAbs().maxCoeff());
            errUnit = std::max(errUnit, std::abs(X.norm() - 1.0));
        }
    }
    CAPTURE(sci(errNodes), sci(errX), sci(errKernel), sci(errUnit));
    CHECK(errNodes < 1e-15);
    CHECK(errX < 1e-12);
    CHECK(errKernel < 1e-12);
    CHECK(errUnit < 1e-12);
}

TEST_CASE("KernelFrame: X' against central differences; no sign flips; X(0) of the reference")
{
    const int n = lift().gridIntervals();
    const double h = lift().length() / n;

    // At the nodes, against the exact derivative d/ds x(q*(s)) = Dx(q*) q*', Dx by central differences: this tests
    // the node derivatives themselves (a second-order stencil would be off by about 1e-5).
    double errNodes = 0.0;
    for (int i = 0; i < 40; ++i) {
        const int k = (i * 211) % n;
        Vec7 X = Vec7::Zero(), dX = Vec7::Zero();
        kernelFrame().eval(k * h, X, &dX);
        const Vec7 q = lift().nodes()[k];
        errNodes = std::max(errNodes, maxDiff(dX, Vec7(kernelJacobian(q) * exactLiftDerivative(k * h, q))));
    }
    // Between the nodes, against a central difference of the X spline (step 1e-6)
    double errSpline = 0.0;
    for (int i = 0; i < 40; ++i) {
        const double s = uniform(0.0, lift().length()), e = 1e-6;
        Vec7 X = Vec7::Zero(), dX = Vec7::Zero(), Xp = Vec7::Zero(), Xm = Vec7::Zero();
        kernelFrame().eval(s, X, &dX);
        kernelFrame().eval(s + e, Xp);
        kernelFrame().eval(s - e, Xm);
        errSpline = std::max(errSpline, maxDiff(dX, Vec7((Xp - Xm) / (2.0 * e))));
    }
    CAPTURE(sci(errNodes), sci(errSpline));
    CHECK(errNodes < 1e-8);
    CHECK(errSpline < 1e-9);

    // No sign flips between neighbouring nodes: min X(s_k)^T X(s_k+1) = 0.99999966
    CHECK(kernelFrame().minNeighbourDot() > 0.9999);
    CHECK_THAT(kernelFrame().minNeighbourDot(), WithinAbs(0.99999966, 1e-8));

    // The orientation is fixed by the formula, so X(0) is the reference implementation's value, sign included.
    Vec7 X0 = Vec7::Zero(), expected;
    kernelFrame().eval(0.0, X0);
    expected << -0.71023893, 0.00352626, 0.68845422, 0.00505790, 0.11836401, -0.01915896, -0.08471320;
    CHECK(maxDiff(X0, expected) < 1e-8);

    Vec7 Xn = Vec7::Zero();
    kernelFrame().eval(std::nan(""), Xn);
    CHECK_FALSE(Xn.allFinite());
}

// ---------------------------------------------------------------------------------------------------------------
// Redundancy outputs (eq:beta-flat, eq:dbeta)

TEST_CASE("LiftRedundancy: beta = 0 and dbeta X = +1 on the lift")
{
    // On the lift, d eta X = 0 (the self-motion does not move the flange), so dbeta X = X^T X = 1 exactly. Not
    // |dbeta X| = 1: that would hide a sign error.
    for (int i = 0; i < 30; ++i) {
        const double s = uniform(0.0, lift().length());
        const Vec7 q = liftAt(s);
        const Outputs o = outputsAt(q, s);
        Vec7 X = Vec7::Zero();
        kernelFrame().eval(o.eta, X);
        double b = 1.0;
        Vec7 db = Vec7::Zero();
        beta().evaluate(q, o.eta, o.deta, b, db);
        const double errProjection = std::abs(o.eta - s) + o.xi.cwiseAbs().maxCoeff();
        const double errOne = std::abs(db.dot(X) - 1.0);
        CAPTURE(s, sci(errProjection), sci(b), sci(errOne));
        CHECK(errProjection < 1e-11);
        CHECK(std::abs(b) < 1e-12);
        CHECK(errOne < 1e-12);
    }
}

TEST_CASE("LiftRedundancy: dbeta against a central difference of beta off the lift")
{
    // Off the lift, q - q* != 0, so the X' term of eq:dbeta is active (on the lift it drops out and a sign error in
    // it would go unseen). eta(q) comes from the tube's projection, as in the controller.
    struct Start {
        double s0;
        Vec5 xi0;
        double beta0;
    };
    const std::vector<Start> starts = {{0.1, paperXi0(), 0.0}, {0.1, paperXi0(), 0.2}, {0.37, paperXi0(), -0.1},
                                       {0.6, Vec5(-paperXi0()), 0.05}};
    const double h = 1e-6;
    for (const Start& st : starts) {
        const Vec7 q = startAt(st.s0, st.xi0, st.beta0);
        const Outputs o = outputsAt(q, st.s0);
        double b = 0.0;
        Vec7 db = Vec7::Zero();
        beta().evaluate(q, o.eta, o.deta, b, db);

        // The size of the X' term, which a sign error would flip
        Vec7 qs = Vec7::Zero(), X = Vec7::Zero(), dX = Vec7::Zero();
        lift().eval(o.eta, qs);
        kernelFrame().eval(o.eta, X, &dX);
        const double xPrimeTerm = std::abs(dX.dot(q - qs)) * o.deta.norm();
        CAPTURE(st.s0, st.beta0, xPrimeTerm);
        CHECK(xPrimeTerm > 0.01);  // big enough for the checks below to see it

        for (int i = 0; i < 5; ++i) {
            const Vec7 w = randn<Vec7>().normalized();
            const double fd = (betaAt(Vec7(q + h * w), st.s0) - betaAt(Vec7(q - h * w), st.s0)) / (2.0 * h);
            const double err = std::abs(fd - db.dot(w));
            CAPTURE(sci(err));
            CHECK(err < 1e-8);
        }
    }
}

TEST_CASE("RedundancyOutput: levelSet is the hyperplane {beta = value}")
{
    // For fixed eta: beta(q) - value = p^T q - b for every q, for both outputs.
    const CouplingRedundancy coupling;
    for (const RedundancyOutput* out : {static_cast<const RedundancyOutput*>(&beta()),
                                        static_cast<const RedundancyOutput*>(&coupling)}) {
        for (int i = 0; i < 10; ++i) {
            const double eta = uniform(0.0, lift().length()), value = uniform(-0.2, 0.2);
            const Vec7 q = liftAt(eta) + randn<Vec7>(0.1);
            double b0 = 0.0, b = 0.0;
            Vec7 db = Vec7::Zero(), p = Vec7::Zero();
            out->evaluate(q, eta, randn<Vec7>(), b0, db);
            out->levelSet(eta, value, p, b);
            const double err = std::abs((b0 - value) - (p.dot(q) - b));
            CAPTURE(eta, value, sci(err));
            CHECK(err < 1e-12);
        }
    }

    // LiftRedundancy: p = X(eta), and the hyperplane contains q*(eta) + value X(eta)
    const double eta = 0.3;
    Vec7 p = Vec7::Zero(), X = Vec7::Zero();
    double b = 0.0;
    beta().levelSet(eta, 0.07, p, b);
    kernelFrame().eval(eta, X);
    CHECK(maxDiff(p, X) == 0.0);
    CHECK(std::abs(p.dot(Vec7(liftAt(eta) + 0.07 * X)) - b) < 1e-14);
}

TEST_CASE("CouplingRedundancy: beta = q_3 - q_1 and dbeta = e_3 - e_1 for every eta")
{
    const CouplingRedundancy coupling;
    const Vec7 e31 = Vec7::Unit(2) - Vec7::Unit(0);
    for (int i = 0; i < 5; ++i) {
        const Vec7 q = randn<Vec7>();
        double b = 0.0;
        Vec7 db = Vec7::Zero();
        coupling.evaluate(q, uniform(0.0, 1.0), randn<Vec7>(), b, db);
        CHECK(b == q[2] - q[0]);
        CHECK(maxDiff(db, e31) == 0.0);
    }
    Vec7 p = Vec7::Zero();
    double b = 0.0;
    coupling.levelSet(0.4, 0.03, p, b);
    CHECK(maxDiff(p, e31) == 0.0);
    CHECK(b == 0.03);
}

// ---------------------------------------------------------------------------------------------------------------
// Start configurations

TEST_CASE("startConfiguration: the requested pose and beta")
{
    struct Start {
        double s0;
        Vec5 xi0;
        double beta0;
    };
    const std::vector<Start> starts = {{0.1, paperXi0(), 0.0}, {0.1, Vec5::Zero(), 0.2},
                                       {0.5, Vec5(-paperXi0()), -0.05}, {0.7, Vec5(0.5 * paperXi0()), 0.1}};
    for (const Start& st : starts) {
        const Vec7 q = startAt(st.s0, st.xi0, st.beta0);
        const Pose target = tube().Psi(st.s0, st.xi0);
        const Outputs o = outputsAt(q, st.s0);
        double b = 0.0;
        Vec7 db = Vec7::Zero();
        beta().evaluate(q, o.eta, o.deta, b, db);
        const double errPose = poseGap(PandaArm::output(q, nullptr), target);
        const double errOutputs = std::abs(o.eta - st.s0) + maxDiff(o.xi, st.xi0);
        CAPTURE(st.s0, st.beta0, sci(errPose), sci(errOutputs), sci(b - st.beta0));
        CHECK(errPose < 1e-10);
        CHECK(errOutputs < 1e-9);
        CHECK(std::abs(b - st.beta0) < 1e-10);
        CHECK(PandaArm::jointMargin(q) > 0.1);
    }

    // With the coupling: the same pose, and q_3 - q_1 = beta0
    const CouplingRedundancy coupling;
    Vec7 q = Vec7::Zero();
    REQUIRE(startConfiguration(tube(), coupling, 0.1, paperXi0(), 0.02, liftAt(0.1), q));
    CHECK(poseGap(PandaArm::output(q, nullptr), tube().Psi(0.1, paperXi0())) < 1e-10);
    CHECK(std::abs(q[2] - q[0] - 0.02) < 1e-10);

    // Out of reach, about 3 m from the path: false, and q is left as it was. (xi = N^T W d is the normal part of
    // the twist d, a translation along the flange's z-axis.)
    Mat65 N = Mat65::Zero();
    tube().frame(0.1, N);
    Vec6 d = Vec6::Zero();
    d[5] = 3.0;
    const Vec5 far = N.transpose() * metric().flat(d);
    const Vec7 before = q;
    CHECK_FALSE(startConfiguration(tube(), beta(), 0.1, far, 0.0, liftAt(0.1), q));
    CHECK(maxDiff(q, before) == 0.0);
}
