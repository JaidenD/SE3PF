#include "se3pf/lift.hpp"
#include "se3pf/se3.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Eigen;

namespace se3pf {

    namespace {
        using Mat7 = Mat<7, 7>;

        constexpr double kTol = 1e-12;           // stopping tolerance of both Newton methods (design doc)
        constexpr int kMaxLiftUpdates = 100;     // projected step: linear convergence, 9 updates per node, 13 at most
        constexpr int kMaxSquareUpdates = 50;    // square system: quadratic convergence, a few updates
        constexpr double kMaxStep = 0.2;         // cap on |dq| for the square system (rad)

        // Pose error e = Log_I(g^-1 target)^vee of the product metric: (log R, p) of the relative pose, in the flange
        // frame. If FK(q + dq) = FK(q) exp((J_b dq)^), then e falls by J_b dq to first order near e = 0, so Newton's
        // methods solve J_b dq = e. Any chart of SE(3) with that property (the SE(3) log, too) has the same root
        // FK(q) = target, and so gives the same lift; only the iterates differ. Written out here instead of calling
        // Metric::Log, because the lift needs no metric: the root does not depend on the chart, and this closed
        // form is the same for every product-type W.
        Vec6 poseError(const Pose& g, const Pose& target)
        {
            const Pose rel = inverse(g) * target;
            Vec6 e;
            e << logSO3(Matrix3d(rel.block<3,3>(0,0))), rel.block<3,1>(0,3);
            return e;
        }

        // Smallest singular value of J_b (Assumption 1: J_b(q*(s)) has full rank)
        double sigmaMin(const Mat67& J)
        {
            const JacobiSVD<Mat67> svd(J);
            return svd.singularValues().minCoeff();
        }

        // One node of the lift: projected Newton steps from q until the stopping rule holds (Offline: the lift q*).
        // Returns the number of updates, or -1 on failure (q is then unspecified).
        int projectedNewton(const Pose& target, const Vec7& qr, Vec7& q)
        {
            for (int k = 0; k <= kMaxLiftUpdates; ++k) {
                Mat67 J = Mat67::Zero();
                const Pose g = PandaArm::output(q, &J);
                const Vec6 e = poseError(g, target);

                // J^+ y = J^T (J J^T)^-1 y; J J^T is symmetric positive definite while J has full rank.
                const LLT<Mat6> JJt(J * J.transpose());
                if (JJt.info() != Success)
                    return -1;
                const Vec7 d = qr - q;
                const Vec7 nullStep = d - J.transpose() * JJt.solve(Vec6(J * d));  // (I - J^+ J)(q_r - q)

                // nullStep = 0 is the Lagrange condition: q - q_r is in the range of J^T.
                if (e.norm() < kTol && nullStep.norm() < kTol)
                    return k;
                if (k == kMaxLiftUpdates)
                    break;
                const Vec7 dq = J.transpose() * JJt.solve(e) + nullStep;
                if (!dq.allFinite())
                    return -1;
                q += dq;
            }
            return -1;
        }

        // f(q+) - 2 f(q*) + f(q-) for f(q) = |q - q_r|^2, with q+- the configurations with pose FK(q*) on the
        // self-motion through q*, at x^T q = x^T q* +- delta (minimumCheck). NaN if q+ or q- is not found.
        double secondDifference(const Vec7& qs, const Vec7& qr, double delta)
        {
            Mat67 J = Mat67::Zero();
            const Pose target = PandaArm::output(qs, &J);  // = gamma(s_k) to 1e-12: the same self-motion
            const Vec7 x = kernelVector(J);
            Vec7 qp = qs, qm = qs;
            if (!constrainedIK(target, x, x.dot(qs) + delta, qp) || !constrainedIK(target, x, x.dot(qs) - delta, qm))
                return std::numeric_limits<double>::quiet_NaN();
            return (qp - qr).squaredNorm() - 2.0 * (qs - qr).squaredNorm() + (qm - qr).squaredNorm();
        }
    }  // namespace

    // ---------------------------------------------------------------------------------------------------------
    // Kinematic building blocks

    Vec7 kernelVector(const Mat67& Jb)
    {
        Vec7 n = Vec7::Zero();
        for (int j = 0; j < 7; ++j) {
            Mat6 M = Mat6::Zero();  // J_b without column j
            M.leftCols(j) = Jb.leftCols(j);
            M.rightCols(6 - j) = Jb.rightCols(6 - j);
            n[j] = (j % 2 == 0 ? 1.0 : -1.0) * M.determinant();
        }
        return n / n.norm();
    }

    bool constrainedIK(const Pose& target, const Vec7& p, double b, Vec7& q, int* iterations)
    {
        Vec7 qk = q;  // the iterate; q is written only on success
        for (int k = 0; k <= kMaxSquareUpdates; ++k) {
            Mat67 J = Mat67::Zero();
            const Pose g = PandaArm::output(qk, &J);
            const Vec6 e = poseError(g, target);
            const double c = b - p.dot(qk);
            if (!e.allFinite() || !std::isfinite(c))
                return false;
            if (e.norm() < kTol && std::abs(c) < kTol) {
                q = qk;
                if (iterations != nullptr)
                    *iterations = k;
                return true;
            }
            if (k == kMaxSquareUpdates)
                break;

            // [J_b; p^T] dq = (e, b - p^T q): six rows move the flange onto the target, the seventh moves q onto
            // the hyperplane, which picks one point of the self-motion manifold.
            Mat7 A;
            A << J, p.transpose();
            Vec7 rhs;
            rhs << e, c;
            Vec7 dq = A.partialPivLu().solve(rhs);  // a singular A gives inf or NaN, caught below
            const double size = dq.norm();
            if (!std::isfinite(size))
                return false;
            if (size > kMaxStep)
                dq *= kMaxStep / size;
            qk += dq;
        }
        return false;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Offline: the lift q* (Assumption 1)

    Lift::Lift(const Path& path, const Vec7& restPosture, int gridIntervals)
        : L_(path.length()), n_(gridIntervals), qr_(restPosture)
    {
        if (n_ < 4)
            throw std::invalid_argument("Lift: need at least 4 grid intervals (fourth-order differences)");
        const double h = L_ / n_;

        // Continuation: node 0 from q_r (a large first step, about 0.97 rad on the paper's path), every later node
        // from the previous solution (steps of about 0.002 rad), so the lift stays on one branch.
        std::vector<Vec7> Q(n_ + 1, Vec7::Zero());
        std::vector<int> updates(n_ + 1, 0);
        Vec7 q = qr_;
        for (int k = 0; k <= n_; ++k) {
            const double s = k * h;
            updates[k] = projectedNewton(path.pose(s), qr_, q);
            if (updates[k] < 0)
                throw std::runtime_error("Lift: Newton's method failed at s = " + std::to_string(s) +
                                         " (is the path reachable?)");
            // Inside the joint limits at every node; a NaN margin fails too.
            int joint = -1;
            if (!(PandaArm::jointMargin(q, &joint) >= 0.0))
                throw std::runtime_error("Lift: joint " + std::to_string(joint) +
                                         " (zero-based) left its range at s = " + std::to_string(s));
            Q[k] = q;
        }

        // Closure on the node values, before the spline replaces node n by node 0.
        diag_.closureGap = (Q[n_] - Q[0]).cwiseAbs().maxCoeff();
        if (!(diag_.closureGap <= 1e-6))
            throw std::runtime_error("Lift: q*(L) and q*(0) differ by " + std::to_string(diag_.closureGap) +
                                     " rad, the lift does not close (Assumption 1 fails on this path)");
        Q[n_] = Q[0];

        // Node derivatives q*' by fourth-order periodic central differences (they agree with the exact derivative
        // from the implicit function theorem to about 1e-9).
        std::vector<Vec7> dQ = PeriodicHermite<Vec7>::centralDifferences(Q, L_);

        // Diagnostics over the nodes k = 0, ..., n - 1 (node n is node 0)
        diag_.minJointMargin = std::numeric_limits<double>::infinity();
        diag_.minSigma = std::numeric_limits<double>::infinity();
        for (int k = 0; k < n_; ++k) {
            const double s = k * h;
            int joint = -1;
            const double margin = PandaArm::jointMargin(Q[k], &joint);
            if (margin < diag_.minJointMargin) {
                diag_.minJointMargin = margin;
                diag_.minMarginJoint = joint;
                diag_.minMarginS = s;
            }
            Mat67 J = Mat67::Zero();
            PandaArm::output(Q[k], &J);
            const double sigma = sigmaMin(J);
            if (sigma < diag_.minSigma) {
                diag_.minSigma = sigma;
                diag_.minSigmaS = s;
            }
            diag_.maxDerivative = std::max(diag_.maxDerivative, dQ[k].norm());
        }

        diag_.firstUpdates = updates[0];
        std::vector<int> later(updates.begin() + 1, updates.end());
        std::sort(later.begin(), later.end());
        diag_.typicalUpdates = later[later.size() / 2];
        diag_.maxUpdates = later.back();

        nodes_ = Q;
        q_.build(L_, std::move(Q), std::move(dQ));
    }

    void Lift::eval(double s, Vec7& q, Vec7* dq) const
    {
        q_.eval(s, q, dq);  // NaN for a non-finite s
    }

    double minimumCheck(const Lift& lift, double delta, double* where)
    {
        const std::vector<Vec7>& Q = lift.nodes();
        const int n = lift.gridIntervals();
        const double h = lift.length() / n;
        double smallest = std::numeric_limits<double>::infinity();
        double at = 0.0;
        for (int k = 0; k < n; ++k) {
            const double d2 = secondDifference(Q[k], lift.restPosture(), delta);
            if (std::isnan(d2)) {  // no later node can make the check pass
                smallest = d2;
                at = k * h;
                break;
            }
            if (d2 < smallest) {
                smallest = d2;
                at = k * h;
            }
        }
        if (where != nullptr)
            *where = at;
        return smallest;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Offline: the kernel frame X (Lemma 3, r = 1)

    KernelFrame::KernelFrame(const Lift& lift)
    {
        const std::vector<Vec7>& Q = lift.nodes();
        const int n = lift.gridIntervals();
        std::vector<Vec7> X(n + 1, Vec7::Zero());
        for (int k = 0; k <= n; ++k) {
            Mat67 J = Mat67::Zero();
            PandaArm::output(Q[k], &J);
            X[k] = kernelVector(J);
            if (!X[k].allFinite())
                throw std::runtime_error("KernelFrame: J_b(q*) is singular at a node");
        }

        // No sign flips between neighbouring nodes: the formula of kernelVector fixes the orientation.
        minDot_ = 1.0;
        for (int k = 0; k < n; ++k)
            minDot_ = std::min(minDot_, X[k].dot(X[k + 1]));

        std::vector<Vec7> dX = PeriodicHermite<Vec7>::centralDifferences(X, lift.length());
        X_.build(lift.length(), std::move(X), std::move(dX));
    }

    void KernelFrame::eval(double s, Vec7& X, Vec7* dX) const
    {
        X_.eval(s, X, dX);  // NaN for a non-finite s
    }

    // ---------------------------------------------------------------------------------------------------------
    // Online: redundancy outputs

    void LiftRedundancy::evaluate(const Vec7& q, double eta, const Vec7& deta, double& beta, Vec7& dbeta) const
    {
        Vec7 qs = Vec7::Zero(), dqs = Vec7::Zero(), X = Vec7::Zero(), dX = Vec7::Zero();
        lift_.eval(eta, qs, &dqs);
        kernel_.eval(eta, X, &dX);
        const Vec7 e = q - qs;
        beta = X.dot(e);                               // eq:beta-flat
        dbeta = X + (dX.dot(e) - X.dot(dqs)) * deta;   // eq:dbeta: q moves directly, and through eta(q)
    }

    void LiftRedundancy::levelSet(double eta, double value, Vec7& p, double& b) const
    {
        Vec7 qs = Vec7::Zero();
        lift_.eval(eta, qs);
        kernel_.eval(eta, p);
        b = p.dot(qs) + value;  // X^T (q - q*) = value  <=>  X^T q = X^T q* + value
    }

    void CouplingRedundancy::evaluate(const Vec7& q, double, const Vec7&, double& beta, Vec7& dbeta) const
    {
        beta = q[2] - q[0];
        dbeta = Vec7::Unit(2) - Vec7::Unit(0);
    }

    void CouplingRedundancy::levelSet(double, double value, Vec7& p, double& b) const
    {
        p = Vec7::Unit(2) - Vec7::Unit(0);
        b = value;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Start configurations

    bool startConfiguration(const Tube& tube, const RedundancyOutput& beta, double s0, const Vec5& xi0, double beta0,
                            const Vec7& guess, Vec7& q)
    {
        Vec7 p = Vec7::Zero();
        double b = 0.0;
        beta.levelSet(s0, beta0, p, b);
        Vec7 qk = guess;  // constrainedIK starts from its argument; q is written only on success
        if (!constrainedIK(tube.Psi(s0, xi0), p, b, qk))
            return false;
        q = qk;
        return true;
    }

}  // namespace se3pf
