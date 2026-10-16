#pragma once
#include "se3pf/types.hpp"
#include "se3pf/redundancy.hpp"
#include "se3pf/tube.hpp"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>

namespace se3pf {

    // Gains of the controller (Section V, Feedback linearization). The defaults are the paper's: transverse and
    // redundancy poles at -lambda with lambda = 4 rad/s (K_p = lambda^2 I = 16 I, K_d = 2 lambda I = 8 I), k_t = 5
    // for a velocity reference, and k_p = 9, k_d = 6 for a pose reference (tangential poles at -3).
    struct Gains {
        double Kp = 16.0, Kd = 8.0;          // transverse, eq:transverse: K_p = Kp I, K_d = Kd I
        double kt = 5.0;                     // tangential velocity, eq:tan-vel
        double kpEta = 9.0, kdEta = 6.0;     // tangential position, eq:tan-pos
        double KpBeta = 16.0, KdBeta = 8.0;  // redundancy (n = 7 only), eq:beta-law

        // Transverse and redundancy poles at -lambda (critically damped): Kp = KpBeta = lambda^2, Kd = KdBeta =
        // 2 lambda. The tangential gains keep their defaults.
        static Gains fromLambda(double lambda)
        {
            Gains g;
            g.Kp = g.KpBeta = lambda * lambda;
            g.Kd = g.KdBeta = 2.0 * lambda;
            return g;
        }
    };

    // What the tangential output eta should do (Section V): follow a velocity nu_d(t, eta) along the path
    // (eq:tan-vel), or a position eta_d(t) on it (eq:tan-pos). eta_d is real-valued (it may run past L); the
    // controller wraps the error eta - eta_d into (-L/2, L/2]. Only the function of the selected mode is used.
    // Calling a std::function does not allocate, so the 1 kHz step may call it; building one may allocate, so do
    // it before the run.
    struct TangentialReference {
        enum class Mode { Velocity, Position };

        // nu = nu_d(t, eta) and its partial derivatives d nu_d / dt and d nu_d / d eta
        using VelocityFn = std::function<void(double t, double eta, double& nu, double& dnuDt, double& dnuDeta)>;
        // etaD = eta_d(t) and its first two time derivatives
        using PositionFn = std::function<void(double t, double& etaD, double& detaD, double& ddetaD)>;

        Mode mode = Mode::Velocity;
        VelocityFn velocity;
        PositionFn position;

        static TangentialReference fromVelocity(VelocityFn f)
        {
            TangentialReference r;
            r.mode = Mode::Velocity;
            r.velocity = std::move(f);
            return r;
        }

        static TangentialReference fromPosition(PositionFn f)
        {
            TangentialReference r;
            r.mode = Mode::Position;
            r.position = std::move(f);
            return r;
        }

        // nu_d = nu, a constant speed along the path
        static TangentialReference constantVelocity(double nu)
        {
            return fromVelocity([nu](double, double, double& v, double& dvDt, double& dvDeta) {
                v = nu;
                dvDt = 0.0;
                dvDeta = 0.0;
            });
        }

        // True if the function of the selected mode is set.
        bool valid() const
        {
            return mode == Mode::Velocity ? static_cast<bool>(velocity) : static_cast<bool>(position);
        }
    };

    // Quintic smoothstep r(x) = 10 x^3 - 15 x^4 + 6 x^5 for x in [0, 1], clamped to 0 below and 1 above. r and its
    // first two derivatives vanish at x = 0, and r' and r'' vanish at x = 1, so a command multiplied by r(t / T)
    // starts at exactly zero and reaches its full value smoothly (Rules for the 1 kHz loop: hardware only).
    inline double quinticRamp(double x)
    {
        if (!(x > 0.0))  // also NaN: no command rather than a NaN command
            return 0.0;
        if (x >= 1.0)
            return 1.0;
        return x * x * x * (10.0 + x * (-15.0 + 6.0 * x));
    }

    // x mod L as the representative in (-L/2, L/2]: the error e_eta of eq:tan-pos between two points of R/LZ.
    inline double wrapSymmetric(double x, double L)
    {
        double d = std::fmod(x, L);  // in (-L, L)
        if (d > 0.5 * L)
            d -= L;
        if (d <= -0.5 * L)
            d += L;
        return d;
    }

    // Number of control periods of a simulated run of length T with period dt: the instants t_k = k dt <= T are
    // k = 0, ..., stepCount(T, dt, maxSteps), so a run has stepCount + 1 rows (20001 for T = 20 s, dt = 1 ms). The
    // 1e-9 allows for rounding, so that t = T is included when T is a multiple of dt. Shared by the simulators
    // (rigid body, Panda, MuJoCo). Offline: throws std::invalid_argument unless T and dt are finite, T >= 0, dt > 0
    // and T / dt <= maxSteps (which also keeps the result in range of long, and turns a mistyped T or dt into an
    // error instead of a huge allocation).
    inline long stepCount(double T, double dt, long maxSteps)
    {
        if (!std::isfinite(T) || !std::isfinite(dt) || !(T >= 0.0) || !(dt > 0.0))
            throw std::invalid_argument("stepCount: need finite T >= 0 and dt > 0");
        if (!(T / dt <= static_cast<double>(maxSteps)))  // also T / dt = inf for a tiny dt
            throw std::invalid_argument("stepCount: T / dt exceeds the step limit of the run");
        return static_cast<long>(std::floor(T / dt + 1e-9));
    }

    // The path following controller of Section V for a fully actuated system with pose output (Algorithm 1).
    //
    // System provides (see RigidBody and PandaArm, and the design doc, The system interface):
    //     static constexpr int kDof;                                       n = 6 (rigid body) or 7 (Panda)
    //     using Config;                                                    q: a Pose, or a Vec7 of joint angles
    //     static Pose output(const Config& q, Mat<6, kDof>* Jb);           Pi(q) and the body Jacobian J_b(q)
    //     static Config retract(const Config& q, const Vec<kDof>& dq, double h);  a curve through q with velocity dq
    //
    // The outputs are ybar = (eta, xi) for n = 6 and ybar = (eta, xi, beta) for n = 7, with ybar_dot = A(q) qdot
    // (eq:A). step returns the generalized acceleration qdd* = A^-1 (mu - Adot qdot) that makes ybar_dd = mu. The
    // caller turns it into the input of eq:fb, u = M qdd* + C qdot + G, so the controller needs no dynamics.
    //
    // The controller keeps a reference to the tube, so the tube (and its path) must outlive it, and a pointer to the
    // redundancy output, which must outlive it too. Construction is offline and may throw; outputs and step are
    // online: fixed-size types only, no allocation, no exceptions (failure is a return value of false).
    template <class System>
    class Controller {
        public:
            static constexpr int n = System::kDof;
            static_assert(n == 6 || n == 7, "Controller: n = 6 (r = 0) or n = 7 (r = 1) generalized velocities");
            using Config = typename System::Config;
            using VecN = Vec<n>;
            using MatN = Mat<n, n>;

            // Outputs at one configuration (eq:A)
            struct OutputState {
                double eta = 0.0;          // tangential output, in [0, L)
                Vec5 xi = Vec5::Zero();    // transverse output
                double beta = 0.0;         // redundancy output (n = 7 only)
                MatN A = MatN::Zero();     // ybar_dot = A qdot; rows (eta, xi_1..xi_5[, beta])
                int iterations = 0;        // Newton updates of the projection
            };

            // What one step saw and did, for logs and tests
            struct Status {
                double t = 0.0;
                double eta = 0.0, etaDot = 0.0;
                double ref = 0.0;          // nu_d(t, eta) (velocity reference) or eta_d(t) (position reference)
                double eTan = 0.0;         // e_v = eta_dot - nu_d, or e_eta = eta - eta_d wrapped into (-L/2, L/2]
                Vec5 xi = Vec5::Zero(), xiDot = Vec5::Zero();
                double beta = 0.0, betaDot = 0.0;  // n = 7 only
                double ramp = 1.0;         // the factor r(t / T) on mu (1 when the ramp is off)
                double sigmaMin = 0.0;     // smallest singular value of A (Proposition 2: positive near the path)
                int iterations = 0;        // Newton updates of the projection at q
                bool ok = false;
            };

            // Throws std::invalid_argument if the reference has no function for its mode, or if the redundancy
            // output is null for n = 7 (an arm with a redundant joint needs beta) or set for n = 6 (no redundancy).
            Controller(const Tube& tube, const Gains& gains, TangentialReference reference,
                       const RedundancyOutput* redundancy = nullptr)
                : tube_(tube), gains_(gains), reference_(std::move(reference)), redundancy_(redundancy)
            {
                if (!reference_.valid())
                    throw std::invalid_argument("Controller: the tangential reference has no function for its mode");
                if (n == 7 && redundancy_ == nullptr)
                    throw std::invalid_argument("Controller: a redundant system (n = 7) needs a redundancy output");
                if (n == 6 && redundancy_ != nullptr)
                    throw std::invalid_argument("Controller: a system with n = 6 has no redundancy output");
            }
            // no temporary tube: the controller keeps a reference to it
            Controller(const Tube&& tube, const Gains& gains, TangentialReference reference,
                       const RedundancyOutput* redundancy = nullptr) = delete;

            const Tube& tube() const { return tube_; }
            const Gains& gains() const { return gains_; }
            const TangentialReference& reference() const { return reference_; }
            const RedundancyOutput* redundancy() const { return redundancy_; }
            double prediction() const { return etaHat_; }  // eta-hat, the start of the next projection
            double rampTime() const { return rampTime_; }

            // Offline: throws std::invalid_argument if the reference has no function for its mode.
            void setReference(TangentialReference reference)
            {
                if (!reference.valid())
                    throw std::invalid_argument("Controller: the tangential reference has no function for its mode");
                reference_ = std::move(reference);
            }

            // mu ramp for hardware-like runs: mu is multiplied by quinticRamp(t / T), so the first command at rest
            // is exactly zero. T = 0 turns it off (the default, and the paper's closed loop). Offline: throws
            // std::invalid_argument unless T is finite and T >= 0.
            void setRampTime(double T)
            {
                if (!(T >= 0.0) || !std::isfinite(T))
                    throw std::invalid_argument("Controller: the ramp time must be finite and >= 0");
                rampTime_ = T;
            }

            // Sets the prediction eta-hat = eta0 wrapped into [0, L). Call it before the first step of every run,
            // with eta0 near the closest point (Tube::globalProject, or s0 for a start pose Psi(s0, xi)).
            void reset(double eta0) { etaHat_ = tube_.wrap(eta0); }

            // ybar and A at q (eq:A), with the projection started at guess. False if the pose has left the tube or
            // a value is not finite.
            bool outputs(const Config& q, double guess, OutputState& out) const
            {
                Mat<6, n> Jb = Mat<6, n>::Zero();
                const Pose g = System::output(q, &Jb);
                if (!tube_.outputs(g, guess, out.eta, out.xi, &out.iterations))
                    return false;

                // J_Psi^-1 J_b by one LU factorization of J_Psi; the inverse is never formed. Its rows are
                // d eta (row 0) and d xi (rows 1..5).
                const Mat6 JPsi = tube_.JPsi(out.eta, out.xi);
                const Mat<6, n> top = JPsi.partialPivLu().solve(Jb);
                if (!top.allFinite())
                    return false;
                out.A.template topRows<6>() = top;

                // eq:A, last row: d beta, which needs d eta = row 0 of J_Psi^-1 J_b (eq:dbeta)
                if constexpr (n == 7) {
                    const Vec7 deta = top.row(0).transpose();
                    Vec7 dbeta = Vec7::Zero();
                    redundancy_->evaluate(q, out.eta, deta, out.beta, dbeta);
                    out.A.row(6) = dbeta.transpose();
                    if (!std::isfinite(out.beta) || !dbeta.allFinite())
                        return false;
                }
                return true;
            }

            // One control cycle (Algorithm 1) at time t with period dt: qdd receives qdd* = A^-1 (mu - Adot qdot).
            // False (and qdd unchanged) if the projection fails, A is nearly singular (sigma_min(A) < 1e-9) or a
            // result is not finite; a simulation then stops, the robot driver damps the arm to rest.
            bool step(const Config& q, const VecN& dq, double t, double dt, VecN& qdd, Status& st)
            {
                st = Status();
                st.t = t;

                // Steps 1-6 and the first part of step 7, all in outputs: g = Pi(q), eta by Newton's method from the
                // prediction, xi, J_Psi, beta (n = 7) and A. Then ybar_dot = A qdot = (eta_dot, xi_dot[, beta_dot]).
                OutputState out;
                if (!outputs(q, etaHat_, out))
                    return false;
                const VecN yDot = out.A * dq;
                const Vec5 xiDot = yDot.template segment<5>(1);

                // Step 7 (continued): Adot qdot = d/dh A(c(h)) qdot at h = 0 along the curve c(h) = retract(q, qdot,
                // h), by a central difference. Both projections start from the eta just found. At rest
                // Adot qdot = 0, and the two projections are skipped.
                VecN Adq = VecN::Zero();
                if (!(dq.array() == 0.0).all()) {
                    constexpr double h = 1e-5;
                    OutputState plus, minus;
                    if (!outputs(System::retract(q, dq, h), out.eta, plus)
                        || !outputs(System::retract(q, dq, -h), out.eta, minus))
                        return false;
                    Adq = ((plus.A - minus.A) * dq) / (2.0 * h);
                }

                // Step 8: mu = (mu_par, mu_transverse[, mu_beta])
                double muPar = 0.0;
                if (reference_.mode == TangentialReference::Mode::Velocity) {
                    if (!reference_.velocity)
                        return false;
                    double nu = 0.0, dnuDt = 0.0, dnuDeta = 0.0;
                    reference_.velocity(t, out.eta, nu, dnuDt, dnuDeta);
                    st.ref = nu;
                    st.eTan = yDot[0] - nu;  // e_v
                    muPar = dnuDt + dnuDeta * yDot[0] - gains_.kt * st.eTan;  // eq:tan-vel
                } else {
                    if (!reference_.position)
                        return false;
                    double etaD = 0.0, detaD = 0.0, ddetaD = 0.0;
                    reference_.position(t, etaD, detaD, ddetaD);
                    st.ref = etaD;
                    st.eTan = wrapSymmetric(out.eta - etaD, tube_.length());  // e_eta in (-L/2, L/2]
                    muPar = ddetaD - gains_.kdEta * (yDot[0] - detaD) - gains_.kpEta * st.eTan;  // eq:tan-pos
                }

                VecN mu = VecN::Zero();
                mu[0] = muPar;
                mu.template segment<5>(1) = -gains_.Kp * out.xi - gains_.Kd * xiDot;  // eq:transverse
                if constexpr (n == 7)
                    mu[6] = -gains_.KpBeta * out.beta - gains_.KdBeta * yDot[6];  // eq:beta-law

                // Optional ramp on mu only: at rest at t = 0, Adot qdot = 0 too, so the first command is zero.
                const double r = rampTime_ > 0.0 ? quinticRamp(t / rampTime_) : 1.0;
                mu *= r;

                // Step 9 (the controller's part of eq:fb): qdd* = A^-1 (mu - Adot qdot). The LU factorization does
                // not see rank loss (a nearly singular A gives huge but finite values), so check sigma_min first.
                const double sigmaMin = Eigen::JacobiSVD<MatN>(out.A).singularValues().minCoeff();
                if (!(sigmaMin >= 1e-9))  // also NaN
                    return false;
                const VecN result = out.A.partialPivLu().solve(mu - Adq);
                if (!result.allFinite())
                    return false;
                qdd = result;

                // Step 9 (end): the prediction for the next cycle, eta-hat = (eta + eta_dot dt) mod L
                etaHat_ = tube_.wrap(out.eta + yDot[0] * dt);

                st.eta = out.eta;
                st.etaDot = yDot[0];
                st.xi = out.xi;
                st.xiDot = xiDot;
                if constexpr (n == 7) {
                    st.beta = out.beta;
                    st.betaDot = yDot[6];
                }
                st.ramp = r;
                st.sigmaMin = sigmaMin;
                st.iterations = out.iterations;
                st.ok = true;
                return true;
            }

        private:
            const Tube& tube_;
            Gains gains_;
            TangentialReference reference_;
            const RedundancyOutput* redundancy_;
            double etaHat_ = 0.0;     // prediction eta-hat (Algorithm 1)
            double rampTime_ = 0.0;   // T of the mu ramp; 0 = off
    };

}  // namespace se3pf
