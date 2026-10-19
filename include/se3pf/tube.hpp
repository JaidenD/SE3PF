#pragma once
#include "se3pf/types.hpp"
#include "se3pf/metric.hpp"
#include "se3pf/path.hpp"
#include "se3pf/spline.hpp"

#include <array>

namespace se3pf {

    // A rotation H in SO(5) written as rotations in at most two orthogonal planes: H u_i = cos(theta_i) u_i +
    // sin(theta_i) v_i and H v_i = -sin(theta_i) u_i + cos(theta_i) v_i, and H is the identity on the rest. Then
    // log H = sum_i theta_i (v_i u_i^T - u_i v_i^T), and exp(t log H) rotates each plane by t theta_i. This is how
    // the tube takes the logarithm Phi of the holonomy (Lemma 1) and evaluates exp(-(s/L) Phi) (eq:holonomy).
    // Why not a plain Phi = log(H): Eigen's matrix logarithm lives in unsupported/ and has no real logarithm for a
    // rotation by pi; the planes handle that case (test "rotationPlanes: a rotation by pi has a logarithm too") and
    // keep a tiny angle such as the Panda's 8.0e-5 (test "Tube: holonomy H is a rotation ..."), and they give
    // exp(-(s/L) Phi) at each of the 2001 grid nodes in closed form, with no matrix exponential.
    struct RotationPlanes {
        int count = 0;  // number of planes, 0 to 2
        std::array<Vec5, 2> u = {Vec5::Zero(), Vec5::Zero()};
        std::array<Vec5, 2> v = {Vec5::Zero(), Vec5::Zero()};
        std::array<double, 2> theta = {0.0, 0.0};  // angles in (0, pi], largest first

        Mat5 log() const;          // Phi = log H, skew-symmetric
        Mat5 exp(double t) const;  // exp(t Phi)
    };

    // Invariant planes of a rotation H (orthogonal, det H = 1), from the eigenvectors of its symmetric part, which
    // has the eigenvalue cos(theta_i) twice for each plane and 1 on the fixed directions. A direction counts as
    // fixed only if 1 - cos(theta) < 1e-12. For a rotation by pi the planes, and so log H, are not unique: any
    // pairing inside the (-1)-eigenspace is returned, which is still a logarithm of H. Offline: throws
    // std::invalid_argument if det H <= 0, and std::runtime_error if H is not orthogonal enough to split.
    RotationPlanes rotationPlanes(const Mat5& H);

    // What the construction of the normal frame found, for tests and diagnostics.
    struct TubeDiagnostics {
        int planes = 0;                              // rotation planes of H (two for both example paths)
        std::array<double, 2> angles = {0.0, 0.0};   // their rotation angles, largest first (0 if absent)
        double detH = 0.0;                           // det H
        double orthogonalityError = 0.0;             // max |H^T H - I|, before projecting H onto SO(5)
        double logError = 0.0;                       // max |exp(Phi) - H_r|, H_r the nearest rotation to H
        double skewError = 0.0;                      // max |Phi + Phi^T|
        double closureGap = 0.0;                     // max |N(s_n) - N(s_0)| on the node values (s_n = L)
    };

    // The tube around a closed path (Section IV): the normal frame N(s) of Lemma 1, the tubular map Psi (eq:Psi),
    // the outputs (eta, xi) = Psi^-1(g) of Proposition 1 and the body Jacobian J_Psi.
    // The tube keeps a reference to the path, so the path must outlive it (build both as named objects).
    // Construction is offline (it allocates and may throw); frame, Psi, projectionCondition, outputs and JPsi are
    // online: they do not allocate and report failure by returning false.
    class Tube {
        public:
            // Builds N(s) on gridIntervals uniform intervals; throws std::runtime_error if the construction fails.
            Tube(const Path& path, const Metric& metric, int gridIntervals = 2000);
            Tube(const Path&& path, const Metric& metric, int gridIntervals = 2000) = delete;  // no temporaries

            double length() const { return L_; }
            int gridIntervals() const { return n_; }
            const Path& path() const { return path_; }
            const Metric& metric() const { return metric_; }
            const TubeDiagnostics& diagnostics() const { return diag_; }

            double wrap(double s) const;  // s mod L, in [0, L)

            // N(s) (Lemma 1) and, if dN is not null, N'(s), by cubic Hermite interpolation of the node values.
            // A non-finite s gives NaN entries, so Psi and JPsi return NaN too (the caller checks allFinite).
            void frame(double s, Mat65& N, Mat65* dN = nullptr) const;

            // Psi(s, xi) = gamma(s) Exp_I(N(s) xi) (eq:Psi)
            Pose Psi(double s, const Vec5& xi) const;

            // e = e_g(s) = Log_I(gamma(s)^-1 g) (eq:projection-error), phi = phi_g(s) = <e, zeta_gamma(s)>_W
            // (eq:projection-condition) and dphi = phi_g'(s). False if Log_I fails or a value is not finite.
            bool projectionCondition(const Pose& g, double s, double& phi, double& dphi, Vec6& e) const;

            // (eta, xi) = Psi^-1(g) (Proposition 1): eta by Newton's method on phi_g from s = guess, then
            // xi = N(eta)^T W e_g(eta). iterations (if not null) receives the number of Newton updates, 1 to 20.
            // False if the pose has left the tube: a value not finite, |phi'| < 1e-12, or no convergence.
            bool outputs(const Pose& g, double guess, double& eta, Vec5& xi, int* iterations = nullptr) const;

            // Closest point for initial conditions (Proposition 1(c)): the nearest of `samples` uniform points of
            // the path, refined by Newton's method. Not for the control loop: it evaluates the path samples times.
            bool globalProject(const Pose& g, double& eta, int samples = 720) const;

            // Body Jacobian of Psi, Psi(y)^-1 dPsi = (J_Psi(y) dy)^ with y = (eta, xi): columns (eta, xi_1..xi_5).
            Mat6 JPsi(double eta, const Vec5& xi) const;

        private:
            // Newton's method on phi_g (Online stage); e receives e_g(eta).
            bool newton(const Pose& g, double guess, double& eta, Vec6& e, int* iterations) const;

            const Path& path_;
            Metric metric_;
            double L_;
            int n_;
            PeriodicHermite<Mat65> N_;
            TubeDiagnostics diag_;
    };

    // Options of the tube-radius estimate (Implementation, Tube radius). The defaults are the paper's test.
    struct RadiusOptions {
        int footPoints = 32;          // values of s, at the midpoints of 32 equal pieces of [0, L)
        int randomDirections = 6;     // unit directions d in R^5, in addition to +-e_1, ..., +-e_5
        unsigned seed = 1;            // of the random directions, so the estimate is reproducible
        int checkPoints = 1000;       // uniform samples of the path searched for a point closer than epsilon
        double relTol = 1e-7;         // "closer" means dist < epsilon (1 - relTol)
        double startFraction = 0.01;  // the doubling starts at epsilon = startFraction * L
        int bisections = 14;
    };

    // Estimate of the tube radius epsilon_0 of Theorem 1: epsilon passes if no sampled point of the path is closer
    // than epsilon to any pose Psi(s, epsilon d) (the test of Proposition 1(c)). Doubles epsilon until it fails,
    // then bisects; returns the largest passing epsilon. Offline, under a second in Release.
    double estimateTubeRadius(const Tube& tube, const RadiusOptions& options = {});

}  // namespace se3pf
