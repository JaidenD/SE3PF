#pragma once
#include "se3pf/types.hpp"
#include "se3pf/controller.hpp"
#include "se3pf/metric.hpp"
#include "se3pf/path.hpp"
#include "se3pf/rigid_body.hpp"
#include "se3pf/tube.hpp"

#include <vector>

namespace se3pf {

    // The paper's rigid-body example (design doc, Rigid-body simulation and Appendix: Rigid body): the free body
    // J = diag(1.0, 1.2, 0.8) kg m^2, m = 2 kg follows the tumbling circle (r = 1 m) in the metric W = I. It starts
    // at rest at Psi(0.1, xi0), xi0 = (0.4, -0.3, 0.2, -0.1, 0.3), with transverse poles at -lambda = -4, under
    // a pose reference (k_p = 9, k_d = 6) or a velocity reference (k_t = 5).
    //
    // The object owns everything the closed loop needs, in dependency order: the metric, the path, the tube (which
    // keeps a reference to the path) and the body. It cannot be copied or moved, since a copy's tube would still
    // refer to the original's path. Construction builds the tube (offline, a few milliseconds in Release).
    class RigidBodyExample {
        public:
            // The tangential law: a position eta_d(t) (eq:tan-pos) or a velocity nu_d(t) (eq:tan-vel), as in
            // TangentialReference::Mode
            enum class Reference { Position, Velocity };

            static constexpr double kS0 = 0.1;       // start point on the path
            static constexpr double kLambda = 4.0;   // transverse poles at -lambda, rad/s
            static constexpr long kMaxSteps = 10'000'000;  // the longest run simulate accepts, in control periods

            RigidBodyExample();
            RigidBodyExample(const RigidBodyExample&) = delete;
            RigidBodyExample& operator=(const RigidBodyExample&) = delete;

            // Not on a temporary example (deleted &&-overloads): a controller keeps a reference to the tube.
            const Metric& metric() const& { return metric_; }
            const Path& path() const& { return path_; }
            const Tube& tube() const& { return tube_; }
            const RigidBody& body() const { return body_; }
            const Metric& metric() const&& = delete;
            const Path& path() const&& = delete;
            const Tube& tube() const&& = delete;

            static Vec5 xi0();                // start offset (0.4, -0.3, 0.2, -0.1, 0.3), |xi0| = 0.6245
            Pose startPose() const;           // g0 = Psi(0.1, xi0); since a pose Psi(s0, xi) has eta = s0, reset(0.1)

            // lambda = 4 for the transverse gains; k_t = 5, k_p = 9, k_d = 6 for the tangential laws
            static Gains gains();

            // eta_d(t) = L/4 + (L/5) (1 - cos(pi t / 10)) for Reference::Position,
            // nu_d(t) = 0.5 + 0.2 sin(0.5 t) for Reference::Velocity
            TangentialReference reference(Reference which) const;

            // The designed response of the transverse output from rest: |xi(t)| = |xi0| (1 + lambda t) e^{-lambda t}
            // (eq:transverse-cl with K_p = lambda^2 I, K_d = 2 lambda I and xi_dot(0) = 0).
            static double designedXiNorm(double t);

        private:
            Metric metric_;
            TumblingCircle path_;
            Tube tube_;
            RigidBody body_;
    };

    // One row of the closed-loop log, at a control instant t_k = k dt (before the plant moves on).
    struct RigidBodyLogRow {
        double t = 0.0;
        double eta = 0.0, etaDot = 0.0;
        double ref = 0.0;            // eta_d (pose reference) or nu_d (velocity reference)
        double eTan = 0.0;           // e_eta or e_v
        Vec5 xi = Vec5::Zero();
        double xiNorm = 0.0;
        double xiDesign = 0.0;       // |xi0| (1 + 4t) e^{-4t}
        Vec6 u = Vec6::Zero();       // body wrench (tau, f), held over [t, t + dt)
    };

    // The log and the summary numbers of one closed-loop run (the Closed loops table of the design doc).
    struct RigidBodyRun {
        RigidBodyExample::Reference reference = RigidBodyExample::Reference::Position;
        double dt = 0.0;
        int substeps = 0;
        std::vector<RigidBodyLogRow> log;

        double startAttitudeDeg = 0.0;   // angle between the attitudes of g0 and gamma(0.1), degrees
        double startPosition = 0.0;      // distance between their positions, m
        double maxDeviation = 0.0;       // max over the log of ||xi| - |xi0| (1 + 4t) e^{-4t}|
        double maxDeviationTime = 0.0;
        double finalXiNorm = 0.0;        // |xi| at the last row
        double finalTanError = 0.0;      // e_eta or e_v at the last row
        double maxInput = 0.0;           // max over the log of |u|_inf
        double maxInputTime = 0.0;
        int maxIterations = 0;           // most Newton updates in one projection
        bool ok = false;                 // false if a step failed (the pose left the tube)
        double failureTime = -1.0;       // time of that step, or -1
    };

    // The paper's closed loop with sample-and-hold: start at rest at g0 = Psi(0.1, xi0) and reset(0.1); at every
    // t_k = k dt, k = 0, ..., stepCount(T, dt) (the instants k dt <= T): step gives qdd*, u = M qdd* + C zeta (eq:fb
    // with G = 0) is held for the period while RigidBody::integrate takes `substeps` equal steps. Stops at the first
    // failed step. Offline: the log is reserved before the loop; throws std::invalid_argument unless T >= 0 and
    // dt > 0 are finite, T / dt <= RigidBodyExample::kMaxSteps and substeps >= 1.
    RigidBodyRun simulate(const RigidBodyExample& example, RigidBodyExample::Reference which, double T = 20.0,
                          double dt = 0.002, int substeps = 4);

}  // namespace se3pf
