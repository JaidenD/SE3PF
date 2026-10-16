// Variant L of the joint-space controller (reviewer's "strongest J"): all coordinates still live in joint space
// around the closed lift q*(s), but they are the linearization of the paper's outputs at the lift, precomputed
// offline along q*:
//     a(s)^T = d eta (row 0 of A(q*(s))),  B(s) = (d xi; d beta) (rows 1..6 of A(q*(s))),
//     eta_L : root of phi(s) = a(s)^T (q - q*(s)) = 0  (the level set of eta, linearized at q*(s)),
//     z     = B(eta_L) (q - q*(eta_L))  in R^6        (xi and beta, linearized at q*(eta_L)).
// So z = (xi, beta) + O(|q - q*|^2) and eta_L = eta + O(|q - q*|^2): the same units as the paper's outputs, the same
// gains act in task space, and the frame carries the task / self-motion split to first order. No forward kinematics,
// SE(3) map or Jacobian is evaluated in the loop. Offline construction allocates; project/outputs/step do not.
#pragma once
#include "joint_tube.hpp"

namespace cc {
    template <class T>
    void nodeDerivativesT(const std::vector<T>& y, double L, std::vector<T>& d1, std::vector<T>& d2)
    {
        const int n = static_cast<int>(y.size()) - 1;
        const double h = L / n;
        auto at = [&](int k) -> const T& { return y[((k % n) + n) % n]; };
        d1.assign(n + 1, T::Zero());
        d2.assign(n + 1, T::Zero());
        for (int k = 0; k < n; ++k) {
            d1[k] = (at(k - 2) - 8.0 * at(k - 1) + 8.0 * at(k + 1) - at(k + 2)) / (12.0 * h);
            d2[k] = (-at(k - 2) + 16.0 * at(k - 1) - 30.0 * at(k) + 16.0 * at(k + 1) - at(k + 2)) / (12.0 * h * h);
        }
        d1[n] = d1[0];
        d2[n] = d2[0];
    }

    struct LinTubeDiagnostics {
        double maxATQ1 = 0;   // max |a^T q*' - 1| at the nodes (0 for an exact lift: eta(q*(s)) = s)
        double maxBQ1 = 0;    // max |B q*'| at the nodes (0 for an exact lift: xi = beta = 0 along q*)
        double minSigmaA = 0; // min sigma_min([a^T; B]) at the nodes
    };

    class LinTube {
        public:
            struct Out {
                double eta = 0;
                Vec6 xi = Vec6::Zero();   // z = (xi_lin, beta_lin)
                Mat7 A = Mat7::Zero();
                double denom = 0;         // a^T q*' - a'^T (q - q*)  (> 0 inside the chart)
                int iterations = 0;
            };

            // Qnodes: lift nodes q*(s_k), k = 0..n (node n = node 0); Anodes: A(q*(s_k)) of the paper's controller.
            LinTube(const std::vector<Vec7>& Qnodes, const std::vector<Mat7>& Anodes, double L)
                : L_(L), n_(static_cast<int>(Qnodes.size()) - 1)
            {
                std::vector<Vec7> Y = Qnodes, dY, ddY, a(n_ + 1), da, dda;
                std::vector<Mat67> B(n_ + 1), dB, ddB;
                Y[n_] = Y[0];
                for (int k = 0; k <= n_; ++k) {
                    const int kk = k == n_ ? 0 : k;
                    a[k] = Anodes[kk].row(0).transpose();
                    B[k] = Anodes[kk].bottomRows<6>();
                }
                nodeDerivativesT(Y, L_, dY, ddY);
                nodeDerivativesT(a, L_, da, dda);
                nodeDerivativesT(B, L_, dB, ddB);
                c_.build(L_, Y, dY, ddY);
                a_.build(L_, a, da, dda);
                B_.build(L_, B, dB, ddB);
                diag_.minSigmaA = std::numeric_limits<double>::infinity();
                for (int k = 0; k < n_; ++k) {
                    diag_.maxATQ1 = std::max(diag_.maxATQ1, std::abs(a[k].dot(dY[k]) - 1.0));
                    diag_.maxBQ1 = std::max(diag_.maxBQ1, (B[k] * dY[k]).cwiseAbs().maxCoeff());
                    diag_.minSigmaA = std::min(diag_.minSigmaA,
                                               Eigen::JacobiSVD<Mat7>(Anodes[k]).singularValues().minCoeff());
                }
            }

            double length() const { return L_; }
            const LinTubeDiagnostics& diagnostics() const { return diag_; }
            double wrap(double s) const
            {
                double x = std::fmod(s, L_);
                if (x < 0.0) x += L_;
                if (x >= L_) x -= L_;
                return x;
            }
            void curve(double s, Vec7& Q, Vec7* Q1 = nullptr) const { c_.eval(s, Q, Q1); }

            bool phiAt(const Vec7& q, double s, double& phi, double& dphi) const
            {
                if (!std::isfinite(s)) return false;
                Vec7 Q, Q1, a, a1;
                c_.eval(s, Q, &Q1);
                a_.eval(s, a, &a1);
                const Vec7 e = q - Q;
                phi = a.dot(e);
                dphi = a1.dot(e) - a.dot(Q1);
                return std::isfinite(phi) && std::isfinite(dphi);
            }

            // Newton on phi from guess, at least one update, |phi| < 1e-10, at most 20 updates (as JointTube)
            bool project(const Vec7& q, double guess, double& eta, int* iterations = nullptr) const
            {
                double s = wrap(guess), phi = 0, dphi = 0;
                if (!phiAt(q, s, phi, dphi)) return false;
                for (int k = 1; k <= 20; ++k) {
                    if (std::abs(dphi) < 1e-12) return false;
                    s = wrap(s - phi / dphi);
                    if (!phiAt(q, s, phi, dphi)) return false;
                    if (std::abs(phi) < 1e-10) {
                        eta = s;
                        if (iterations) *iterations = k;
                        return true;
                    }
                }
                return false;
            }

            bool globalProject(const Vec7& q, double& eta) const
            {
                double best = std::numeric_limits<double>::infinity(), sBest = 0;
                for (int k = 0; k < 4 * n_; ++k) {
                    const double s = L_ * k / (4.0 * n_);
                    Vec7 Q;
                    c_.eval(s, Q);
                    const double d = (q - Q).norm();
                    if (d < best) { best = d; sBest = s; }
                }
                return project(q, sBest, eta);
            }

            bool outputs(const Vec7& q, double guess, Out& out) const
            {
                if (!project(q, guess, out.eta, &out.iterations)) return false;
                Vec7 Q, Q1, a, a1;
                Mat67 B, B1;
                c_.eval(out.eta, Q, &Q1);
                a_.eval(out.eta, a, &a1);
                B_.eval(out.eta, B, &B1);
                const Vec7 e = q - Q;
                out.denom = a.dot(Q1) - a1.dot(e);
                const Eigen::Matrix<double, 1, 7> deta = a.transpose() / out.denom;
                out.xi = B * e;
                out.A.row(0) = deta;
                out.A.bottomRows<6>() = B + (B1 * e - B * Q1) * deta;
                return out.xi.allFinite() && out.A.allFinite() && out.denom > 0.0;
            }

        private:
            double L_;
            int n_;
            PeriodicQuintic<Vec7> c_, a_;
            PeriodicQuintic<Mat67> B_;
            LinTubeDiagnostics diag_;
    };

    // The same loop as JointController (laws, gains, prediction, Adot by a central difference, sigma_min check),
    // over any tube with outputs(q, guess, Out) and wrap(s).
    template <class TubeT>
    class GenericJointController {
        public:
            using Status = JointController::Status;
            GenericJointController(const TubeT& tube, const Gains& g, double nu) : tube_(tube), g_(g), nu_(nu) {}
            void reset(double eta0) { etaHat_ = tube_.wrap(eta0); }

            bool step(const Vec7& q, const Vec7& dq, double /*t*/, double dt, Vec7& qdd, Status& st)
            {
                st = Status();
                typename TubeT::Out out;
                if (!tube_.outputs(q, etaHat_, out)) return false;
                const Vec7 yDot = out.A * dq;
                const Vec6 xiDot = yDot.tail<6>();
                Vec7 Adq = Vec7::Zero();
                if (!(dq.array() == 0.0).all()) {
                    constexpr double h = 1e-5;
                    typename TubeT::Out plus, minus;
                    if (!tube_.outputs(q + h * dq, out.eta, plus) || !tube_.outputs(q - h * dq, out.eta, minus))
                        return false;
                    Adq = ((plus.A - minus.A) * dq) / (2.0 * h);
                }
                Vec7 mu;
                st.eTan = yDot[0] - nu_;
                mu[0] = -g_.kt * st.eTan;
                mu.tail<6>() = -g_.Kp * out.xi - g_.Kd * xiDot;
                const double sigmaMin = Eigen::JacobiSVD<Mat7>(out.A).singularValues().minCoeff();
                if (!(sigmaMin >= 1e-9)) return false;
                const Vec7 r = out.A.partialPivLu().solve(mu - Adq);
                if (!r.allFinite()) return false;
                qdd = r;
                etaHat_ = tube_.wrap(out.eta + yDot[0] * dt);
                st.eta = out.eta;
                st.etaDot = yDot[0];
                st.xi = out.xi;
                st.xiDot = xiDot;
                st.sigmaMin = sigmaMin;
                st.denom = out.denom;
                st.iterations = out.iterations;
                st.ok = true;
                return true;
            }

        private:
            const TubeT& tube_;
            Gains g_;
            double nu_;
            double etaHat_ = 0;
    };

}  // namespace cc
