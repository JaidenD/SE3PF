#pragma once
#include "se3pf/types.hpp"
#include "se3pf/panda.hpp"
#include "se3pf/path.hpp"
#include "se3pf/redundancy.hpp"
#include "se3pf/spline.hpp"
#include "se3pf/tube.hpp"

#include <vector>

namespace se3pf {

    // ---------------------------------------------------------------------------------------------------------
    // Kinematic building blocks for the Panda (n = 7, r = 1)

    // Unit vector spanning ker J_b of a 6 x 7 Jacobian of rank 6 (Kernel frame, r = 1): the normalized generalized
    // cross product of the six rows, n_j = (-1)^j det(J_b without column j), j = 0, ..., 6, and x = n / |n|.
    // J_b n = 0, because entry i of J_b n is the determinant of J_b with row i repeated on top. More generally
    // w^T n = det([w^T; J_b]) for every w, so the formula fixes the orientation: det([x^T; J_b]) = |n| > 0. x is a
    // smooth function of J_b, so along a closed lift it is periodic with no sign to choose. NaN if the rank is below
    // 6 (then n = 0). Fixed-size, no allocation.
    Vec7 kernelVector(const Mat67& Jb);

    // Newton's method on the square system [J_b(q); p^T] dq = (e, b - p^T q) for FK(q) = target and p^T q = b,
    // started from q (Start configurations). e = Log_I(FK(q)^-1 target) is the pose error in the flange frame (the
    // product-metric Log_I). Each update is capped at 0.2 rad. Stops when |e| < 1e-12 and |p^T q - b| < 1e-12, and
    // then writes the solution to q; iterations (if not null) receives the number of updates. False, with q
    // unchanged, if the system is singular or there is no convergence in 50 updates. Fixed-size, no allocation.
    bool constrainedIK(const Pose& target, const Vec7& p, double b, Vec7& q, int* iterations = nullptr);

    // ---------------------------------------------------------------------------------------------------------
    // Offline: the lift q* (Assumption 1)

    // What the construction of the lift found, for tests and diagnostics. Everything is measured on the node
    // values s_k = k L / n.
    struct LiftDiagnostics {
        double closureGap = 0.0;           // max |q*(s_n) - q*(s_0)| (s_n = L), before node n is replaced by node 0
        double minJointMargin = 0.0;       // smallest distance of a joint to its nearer limit, rad
        int minMarginJoint = -1;           // that joint (zero-based)
        double minMarginS = 0.0;           // and where
        double minSigma = 0.0;             // smallest singular value of J_b(q*(s_k)) (Assumption 1 needs it > 0)
        double minSigmaS = 0.0;
        double maxDerivative = 0.0;        // largest |q*'(s_k)| (2-norm), rad per unit arc length
        int firstUpdates = 0;              // Newton updates at s = 0, started from q_r
        int typicalUpdates = 0;            // median over the later nodes, each started from the previous solution
        int maxUpdates = 0;                // largest over the later nodes
    };

    // The lift q*(s) of Assumption 1 for the Panda: on the self-motion manifold over gamma(s), the configuration
    // nearest the rest posture q_r in the norm |q - q_r|_W^2 = (q - q_r)^T W (q - q_r), W = diag(weights) > 0
    // (the paper's lift: W = I). Computed on the uniform grid s_k = k L / n by continuation: node 0 from q_r,
    // every later node from the previous solution, by the projected Newton step
    //     dq = J^+ e + (I - J^+ J)(q_r - q),   J = J_b(q),  J^+ = W^-1 J^T (J W^-1 J^T)^-1,
    // until |e| < 1e-12 and |(I - J^+ J)(q - q_r)| < 1e-12, the Lagrange condition of min |q - q_r|_W subject to
    // FK(q) = gamma(s): W (q - q_r) lies in the range of J^T. Node derivatives q*' by fourth-order periodic central
    // differences; q* and q*' are stored in a cubic Hermite spline. A weight above 1 on a joint keeps it nearer
    // its rest value and moves the self-motion to the others: on the robot (CTRL-1, 2026-10-19) a weight on joint 7
    // hands the flange's yaw from joint 7 to joint 1, which the paper's lift leaves almost still (0.28 rad of travel
    // per loop against joint 7's 1.3), so still that its static friction made it stick-slip.
    // The lift copies what it needs, so the path may be destroyed after construction. Construction is offline:
    // it throws std::invalid_argument if gridIntervals < 4 (the fourth-order differences need two neighbours on
    // each side) or a weight is not finite and positive, and std::runtime_error if Newton's method fails at a
    // node, a joint leaves its range, or the branch does not close (|q*(L) - q*(0)| > 1e-6: the continuation
    // ended on another self-motion, and Assumption 1 fails). Whether q* is a minimum and not a maximum is a
    // separate check, minimumCheck below. eval is online: no allocation, no exceptions.
    class Lift {
        public:
            explicit Lift(const Path& path, const Vec7& restPosture = PandaArm::restPosture(),
                          int gridIntervals = 2000, const Vec7& weights = Vec7::Ones());

            // q*(s) and, if dq is not null, q*'(s), by cubic Hermite interpolation of the node values; s is
            // wrapped into [0, L). A non-finite s gives NaN entries (the controller's allFinite checks catch them).
            void eval(double s, Vec7& q, Vec7* dq = nullptr) const;

            double length() const { return L_; }
            int gridIntervals() const { return n_; }
            const Vec7& restPosture() const { return qr_; }
            const Vec7& weights() const { return w_; }  // the diagonal of W
            // q*(s_k) for k = 0, ..., n; node n is node 0 (the spline's convention)
            const std::vector<Vec7>& nodes() const { return nodes_; }
            const LiftDiagnostics& diagnostics() const { return diag_; }

        private:
            double L_;
            int n_;
            Vec7 qr_;
            Vec7 w_;
            std::vector<Vec7> nodes_;
            PeriodicHermite<Vec7> q_;
            LiftDiagnostics diag_;
    };

    // Is q* a minimum of f(q) = |q - q_r|_W^2 on each self-motion, and not a maximum? At every node s_k, q+- are the
    // configurations with the pose FK(q*(s_k)) at x^T q = x^T q*(s_k) +- delta, x = kernelVector(J_b(q*(s_k))),
    // found by constrainedIK from q*(s_k). Returns the smallest second difference f(q+) - 2 f(q*) + f(q-) over the
    // nodes: positive for a minimum. (q* +- delta x itself would not do: f is quadratic, so its second difference
    // along a straight line is 2 delta^2 whatever q* is.) where, if not null, receives the s_k of the result. NaN
    // if q+ or q- cannot be found at some node, so a check of the value fails. Offline, and kept out of the Lift
    // constructor because the controller does not need it: 2n constrained solves (about 10 ms in Release, 2 s in
    // Debug).
    double minimumCheck(const Lift& lift, double delta = 0.05, double* where = nullptr);

    // ---------------------------------------------------------------------------------------------------------
    // Offline: the kernel frame X (Lemma 3, r = 1)

    // X(s) spans ker J_b(q*(s)): for r = 1 no transport is needed (Implementation, Lift), X(s_k) is
    // kernelVector(J_b(q*(s_k))) at the lift's nodes, with X' by fourth-order periodic central differences, stored
    // in a cubic Hermite spline. It depends only on q*(s), so it is periodic (no holonomy to cancel). Copies what
    // it needs from the lift. Construction throws std::runtime_error if J_b(q*) loses rank at a node; eval is
    // online.
    class KernelFrame {
        public:
            explicit KernelFrame(const Lift& lift);

            // X(s) and, if dX is not null, X'(s); a non-finite s gives NaN entries.
            void eval(double s, Vec7& X, Vec7* dX = nullptr) const;

            // min over k of X(s_k)^T X(s_k+1): close to 1 when no neighbouring nodes have opposite signs
            double minNeighbourDot() const { return minDot_; }

        private:
            PeriodicHermite<Vec7> X_;
            double minDot_ = 0.0;
    };

    // ---------------------------------------------------------------------------------------------------------
    // Online: redundancy outputs (eq:beta-flat, eq:dbeta)

    // beta(q) = X(eta)^T (q - q*(eta)) (eq:beta-flat), the signed distance from the lift along the self-motion
    // direction, with
    //     dbeta(q) = X(eta)^T + (X'(eta)^T (q - q*(eta)) - X(eta)^T q*'(eta)) deta(q)   (eq:dbeta).
    // Keeps references to the lift and the kernel frame, so both must outlive it.
    class LiftRedundancy final : public RedundancyOutput {
        public:
            LiftRedundancy(const Lift& lift, const KernelFrame& kernel) : lift_(lift), kernel_(kernel) {}
            LiftRedundancy(const Lift&& lift, const KernelFrame& kernel) = delete;   // no temporaries
            LiftRedundancy(const Lift& lift, const KernelFrame&& kernel) = delete;
            LiftRedundancy(const Lift&& lift, const KernelFrame&& kernel) = delete;

            void evaluate(const Vec7& q, double eta, const Vec7& deta, double& beta, Vec7& dbeta) const override;

            // p = X(eta), b = X(eta)^T q*(eta) + value
            void levelSet(double eta, double value, Vec7& p, double& b) const override;

        private:
            const Lift& lift_;
            const KernelFrame& kernel_;  // X of Lemma 3 (not the normal frame N of the tube)
    };

    // The joint coupling beta = q_3 - q_1 (zero-based q[2] - q[0]), a redundancy output for a lift with q_3 = q_1
    // (the coupled six-DOF arm). Linear in q: dbeta = (e_3 - e_1)^T, and the level set {beta = value} is
    // p^T q = b with p = e_3 - e_1, b = value. Does not depend on eta.
    class CouplingRedundancy final : public RedundancyOutput {
        public:
            void evaluate(const Vec7& q, double eta, const Vec7& deta, double& beta, Vec7& dbeta) const override;
            void levelSet(double eta, double value, Vec7& p, double& b) const override;
    };

    // ---------------------------------------------------------------------------------------------------------
    // Start configurations

    // A configuration q with pose Psi(s0, xi0) and beta(q) = beta0: with (p, b) from beta.levelSet(s0, beta0),
    // constrainedIK solves FK(q) = Psi(s0, xi0) and p^T q = b from guess (usually q*(s0)). Correct because
    // eta = s0 for that pose (Proposition 1). False if Newton's method fails; q is written only on success.
    bool startConfiguration(const Tube& tube, const RedundancyOutput& beta, double s0, const Vec5& xi0, double beta0,
                            const Vec7& guess, Vec7& q);

}  // namespace se3pf
