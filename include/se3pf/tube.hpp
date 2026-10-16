#pragma once
#include "se3pf/types.hpp"
#include "se3pf/metric.hpp"
#include "se3pf/path.hpp"
#include "se3pf/spline.hpp"

#include <array>

namespace se3pf {
    // A rotation H in SO(5) can be written as the rotation of at most two orthogonal planes:
    // H u_i = cos(theta_i) u_i + sin(theta_i) v_i
    // H v_i = -sin(theta_i) u_i + cos(theta_i) v_i
    struct RotationPlanes {
        int count = 0;  // number of planes, 0 to 2
        std::array<Vec5, 2> u = {Vec5::Zero(), Vec5::Zero()};
        std::array<Vec5, 2> v = {Vec5::Zero(), Vec5::Zero()};
        std::array<double, 2> theta = {0.0, 0.0};  // angles in (0, pi], largest first

        Mat5 log() const;          // Phi = log H, skew-symmetric
        Mat5 exp(double t) const;  // exp(t Phi)
    };

    RotationPlanes rotationPlanes(const Mat5& H);

    // What the construction of the normal frame found, for tests and diagnostics.
    struct TubeDiagnostics {
        int planes = 0;                              // rotation planes of H
        std::array<double, 2> angles = {0.0, 0.0};   // their rotation angles, largest first
        double detH = 0.0;                           // det H
        double orthogonalityError = 0.0;             // max |H^T H - I|
        double logError = 0.0;                       // max |exp(Phi) - H_r|
        double skewError = 0.0;                      // max |Phi + Phi^T|
        double closureGap = 0.0;                     // max |N(s_n) - N(s_0)|
    };

    // The tube around a closed path (Section IV): the normal frame N(s) of Lemma 1, the tubular map Psi (eq:Psi),
    // the outputs (eta, xi) = Psi^-1(g) of Proposition 1 and the body Jacobian J_Psi.
    class Tube {
        public:
            // Build N(s) on gridIntervals many uniform intervals
            Tube(const Path& path, const Metric& metric, int gridIntervals = 2000);
            Tube(const Path&& path, const Metric& metric, int gridIntervals = 2000) = delete;  // no temporaries

            double length() const { return L_; }
            int gridIntervals() const { return n_; }
            
            const Path& path() const { return path_; }
            const Metric& metric() const { return metric_; }
            const TubeDiagnostics& diagnostics() const { return diag_; }

            double wrap(double s) const;  // s mod L, in [0, L)

            // N(s) and N'(s), by cubic Hermite interpolation of the node values.
            void frame(double s, Mat65& N, Mat65* dN = nullptr) const;

            // Psi(s, xi) = gamma(s) Exp_I(N(s) xi)
            Pose Psi(double s, const Vec5& xi) const;

            // e = e_g(s) = Log_I(gamma(s)^-1 g) (eq:projection-error), phi = phi_g(s) = <e, zeta_gamma(s)>_W
            // (eq:projection-condition) and dphi = phi_g'(s). Returns false if Log_I fails.
            bool projectionCondition(const Pose& g, double s, double& phi, double& dphi, Vec6& e) const;

            // (eta, xi) = Psi^-1(g). eta is computed by newtons method on phi_g from guessed s value,
            // and then xi is computed as xi = N(eta)^T W e_g(eta).
            // If iterations isnt null then it receives the number of newton updates.
            // Returns false if the pose has left the tube.
            bool outputs(const Pose& g, double guess, double& eta, Vec5& xi, int* iterations = nullptr) const;

            // Find point on the path closest to pose g when you have no starting guess.
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
    // Estimste the tube radius by checking if any sampled points of the path are closer than epsilon to any pose
    // Psi(s, epsilon d). Double epsilon until it fails, then bisect.
    double estimateTubeRadius(const Tube& tube, const RadiusOptions& options = {});

}  // namespace se3pf
