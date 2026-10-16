// Alternative controller J: joint-space tubular coordinates around a closed curve q*(s) in R^7 (Euclidean metric).
// eta_J = closest-point parameter, xi_J = N_J(eta_J)^T (q - q*(eta_J)) in R^6, N_J a closed rotation-minimizing
// normal frame (Lemma 1 transcribed to R^7). Offline construction allocates; project/outputs/step do not.
#pragma once
#include "se3pf/types.hpp"
#include "se3pf/spline.hpp"
#include "se3pf/controller.hpp"

#include <unsupported/Eigen/MatrixFunctions>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cc {
    using namespace se3pf;
    using Mat7 = Mat<7, 7>;
    using Mat76 = Mat<7, 6>;

    // Periodic quintic Hermite interpolation on a uniform grid (C2: value, first and second derivative continuous).
    template <class T>
    class PeriodicQuintic {
        public:
            void build(double L, std::vector<T> y, std::vector<T> d1, std::vector<T> d2)
            {
                n_ = static_cast<int>(y.size()) - 1;
                L_ = L;
                h_ = L / n_;
                y_ = y; d1_ = d1; d2_ = d2;
                y_[n_] = y_[0]; d1_[n_] = d1_[0]; d2_[n_] = d2_[0];
            }

            void eval(double s, T& y, T* dy = nullptr, T* ddy = nullptr) const
            {
                if (!std::isfinite(s)) {
                    y.setConstant(std::numeric_limits<double>::quiet_NaN());
                    if (dy) dy->setConstant(std::numeric_limits<double>::quiet_NaN());
                    if (ddy) ddy->setConstant(std::numeric_limits<double>::quiet_NaN());
                    return;
                }
                double x = std::fmod(s, L_);
                if (x < 0.0) x += L_;
                int k = static_cast<int>(x / h_);
                if (k >= n_) k = n_ - 1;
                const double t = x / h_ - k, t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
                const double h = h_, h2 = h * h;
                const T &y0 = y_[k], &y1 = y_[k + 1], &a0 = d1_[k], &a1 = d1_[k + 1], &c0 = d2_[k], &c1 = d2_[k + 1];
                const double H0 = 1 - 10 * t3 + 15 * t4 - 6 * t5, H1 = t - 6 * t3 + 8 * t4 - 3 * t5,
                             H2 = 0.5 * t2 - 1.5 * t3 + 1.5 * t4 - 0.5 * t5, H3 = 10 * t3 - 15 * t4 + 6 * t5,
                             H4 = -4 * t3 + 7 * t4 - 3 * t5, H5 = 0.5 * t3 - t4 + 0.5 * t5;
                y = H0 * y0 + (h * H1) * a0 + (h2 * H2) * c0 + H3 * y1 + (h * H4) * a1 + (h2 * H5) * c1;
                if (dy) {
                    const double D0 = -30 * t2 + 60 * t3 - 30 * t4, D1 = 1 - 18 * t2 + 32 * t3 - 15 * t4,
                                 D2 = t - 4.5 * t2 + 6 * t3 - 2.5 * t4, D3 = -D0,
                                 D4 = -12 * t2 + 28 * t3 - 15 * t4, D5 = 1.5 * t2 - 4 * t3 + 2.5 * t4;
                    *dy = (D0 / h) * y0 + D1 * a0 + (h * D2) * c0 + (D3 / h) * y1 + D4 * a1 + (h * D5) * c1;
                }
                if (ddy) {
                    const double E0 = -60 * t + 180 * t2 - 120 * t3, E1 = -36 * t + 96 * t2 - 60 * t3,
                                 E2 = 1 - 9 * t + 18 * t2 - 10 * t3, E3 = -E0,
                                 E4 = -24 * t + 84 * t2 - 60 * t3, E5 = 3 * t - 12 * t2 + 10 * t3;
                    *ddy = (E0 / h2) * y0 + (E1 / h) * a0 + E2 * c0 + (E3 / h2) * y1 + (E4 / h) * a1 + E5 * c1;
                }
            }

        private:
            double L_ = 1, h_ = 1;
            int n_ = 0;
            std::vector<T> y_, d1_, d2_;
    };

    // Fourth-order periodic central differences (first and second derivative) of node values y_0..y_n, y_n = y_0.
    inline void nodeDerivatives(const std::vector<Vec7>& y, double L, std::vector<Vec7>& d1, std::vector<Vec7>& d2)
    {
        const int n = static_cast<int>(y.size()) - 1;
        const double h = L / n;
        auto at = [&](int k) -> const Vec7& { return y[((k % n) + n) % n]; };
        d1.assign(n + 1, Vec7::Zero());
        d2.assign(n + 1, Vec7::Zero());
        for (int k = 0; k < n; ++k) {
            d1[k] = (at(k - 2) - 8.0 * at(k - 1) + 8.0 * at(k + 1) - at(k + 2)) / (12.0 * h);
            d2[k] = (-at(k - 2) + 16.0 * at(k - 1) - 30.0 * at(k) + 16.0 * at(k + 1) - at(k + 2)) / (12.0 * h * h);
        }
        d1[n] = d1[0];
        d2[n] = d2[0];
    }

    struct JointTubeDiagnostics {
        double closureGap = 0;        // max |N_J(L) - N_J(0)| on the node values, before the spline closes them
        double holonomyAngle = 0;     // largest rotation angle of H in SO(6)
        double detH = 0;
        double orthH = 0;             // max |H^T H - I|
        double logError = 0;          // max |exp(Phi) - H_r|
        double maxOrthoNodes = 0;     // max |N^T N - I| over nodes
        double maxNormalNodes = 0;    // max |N^T t| over nodes (t = q*'/|q*'|)
        double maxOrthoInterp = 0;    // the same at off-node points of the interpolant
        double maxNormalInterp = 0;
        double minSpeed = 0, maxSpeed = 0;  // min/max |q*'| (rad per unit of s)
        double maxCurvature = 0;            // max kappa = |t'| / |q*'| (1/rad): focal distance 1/kappa
    };

    // The tube around the closed joint curve.
    class JointTube {
        public:
            struct Out {
                double eta = 0;
                Vec6 xi = Vec6::Zero();
                Mat7 A = Mat7::Zero();   // rows (eta, xi_1..xi_6)
                double denom = 0;        // |q*'|^2 - (q - q*)^T q*''  (> 0 inside the tube)
                int iterations = 0;
            };

            JointTube(const std::vector<Vec7>& nodes, double L) : L_(L), n_(static_cast<int>(nodes.size()) - 1)
            {
                std::vector<Vec7> Y = nodes, d1, d2;
                Y[n_] = Y[0];
                nodeDerivatives(Y, L_, d1, d2);
                c_.build(L_, Y, d1, d2);
                const double h = L_ / n_;

                // initial frame: Gram-Schmidt on t(0), e_1..e_7
                std::vector<Mat76> N0(n_ + 1, Mat76::Zero());
                {
                    Vec7 Q, Q1;
                    c_.eval(0.0, Q, &Q1);
                    std::vector<Vec7> basis;
                    auto add = [&](Vec7 x) {
                        for (const Vec7& b : basis) x -= b.dot(x) * b;
                        if (x.norm() > 1e-8) basis.push_back(x.normalized());
                    };
                    add(Q1);
                    for (int i = 0; i < 7; ++i) add(Vec7::Unit(i));
                    if (basis.size() != 7) throw std::runtime_error("JointTube: Gram-Schmidt failed");
                    for (int j = 0; j < 6; ++j) N0[0].col(j) = basis[j + 1];
                }
                // rotation-minimizing transport N0' = -t t'^T N0 by RK4, one step per grid interval
                for (int k = 0; k < n_; ++k) {
                    const double s = k * h;
                    const Mat76& N = N0[k];
                    const Mat76 k1 = rate(s, N);
                    const Mat76 k2 = rate(s + 0.5 * h, N + (0.5 * h) * k1);
                    const Mat76 k3 = rate(s + 0.5 * h, N + (0.5 * h) * k2);
                    const Mat76 k4 = rate(s + h, N + h * k3);
                    N0[k + 1] = N + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                }
                // holonomy H in SO(6), Phi = log H_r
                const Mat6 H = N0[0].transpose() * N0[n_];
                diag_.detH = H.determinant();
                diag_.orthH = (H.transpose() * H - Mat6::Identity()).cwiseAbs().maxCoeff();
                const Eigen::JacobiSVD<Mat6> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
                const Mat6 Hr = svd.matrixU() * svd.matrixV().transpose();
                Mat6 Phi = Hr.log();
                Phi = 0.5 * (Phi - Phi.transpose());
                diag_.logError = (Mat6(Phi.exp()) - Hr).cwiseAbs().maxCoeff();
                {
                    const Eigen::EigenSolver<Mat6> es(Phi);
                    diag_.holonomyAngle = es.eigenvalues().imag().cwiseAbs().maxCoeff();
                }
                Phi_ = Phi;
                std::vector<Mat76> Yn(n_ + 1), Dn(n_ + 1);
                for (int k = 0; k <= n_; ++k) {
                    const double s = k * h;
                    const Mat6 E = Mat6(-(s / L_) * Phi).exp();
                    Yn[k] = N0[k] * E;
                    Dn[k] = rate(s, Yn[k]) - (1.0 / L_) * (Yn[k] * Phi);
                }
                diag_.closureGap = (Yn[n_] - Yn[0]).cwiseAbs().maxCoeff();
                for (int k = 0; k < n_; ++k) {
                    Vec7 Q, Q1;
                    c_.eval(k * h, Q, &Q1);
                    diag_.maxOrthoNodes = std::max(diag_.maxOrthoNodes,
                        (Yn[k].transpose() * Yn[k] - Mat6::Identity()).cwiseAbs().maxCoeff());
                    diag_.maxNormalNodes = std::max(diag_.maxNormalNodes, (Yn[k].transpose() * Q1.normalized()).cwiseAbs().maxCoeff());
                }
                N_.build(L_, Yn, Dn);

                // off-node checks and curve statistics
                diag_.minSpeed = std::numeric_limits<double>::infinity();
                for (int k = 0; k < 4 * n_; ++k) {
                    const double s = (k + 0.37) * h / 4.0;
                    Vec7 Q, Q1, Q2;
                    c_.eval(s, Q, &Q1, &Q2);
                    Mat76 N;
                    N_.eval(s, N);
                    diag_.maxOrthoInterp = std::max(diag_.maxOrthoInterp, (N.transpose() * N - Mat6::Identity()).cwiseAbs().maxCoeff());
                    diag_.maxNormalInterp = std::max(diag_.maxNormalInterp, (N.transpose() * Q1.normalized()).cwiseAbs().maxCoeff());
                    const double v = Q1.norm();
                    diag_.minSpeed = std::min(diag_.minSpeed, v);
                    diag_.maxSpeed = std::max(diag_.maxSpeed, v);
                    const Vec7 t = Q1 / v;
                    const Vec7 tp = (Q2 - t * t.dot(Q2)) / v;
                    diag_.maxCurvature = std::max(diag_.maxCurvature, tp.norm() / v);
                }
            }

            double length() const { return L_; }
            int gridIntervals() const { return n_; }
            const JointTubeDiagnostics& diagnostics() const { return diag_; }

            double wrap(double s) const
            {
                double x = std::fmod(s, L_);
                if (x < 0.0) x += L_;
                if (x >= L_) x -= L_;
                return x;
            }

            void curve(double s, Vec7& Q, Vec7* Q1 = nullptr, Vec7* Q2 = nullptr) const { c_.eval(s, Q, Q1, Q2); }
            void frame(double s, Mat76& N, Mat76* dN = nullptr) const { N_.eval(s, N, dN); }

            // phi(s) = (q - q*(s))^T q*'(s), dphi = -|q*'|^2 + (q - q*)^T q*''
            bool phiAt(const Vec7& q, double s, double& phi, double& dphi) const
            {
                if (!std::isfinite(s)) return false;
                Vec7 Q, Q1, Q2;
                c_.eval(s, Q, &Q1, &Q2);
                const Vec7 e = q - Q;
                phi = e.dot(Q1);
                dphi = -Q1.squaredNorm() + e.dot(Q2);
                return std::isfinite(phi) && std::isfinite(dphi);
            }

            // Newton on phi from guess; at least one update (as Tube::outputs); |phi| < 1e-10; at most 20 updates.
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

            // (eta_J, xi_J) and A_J = d(eta_J, xi_J)/dq (exact derivative of the interpolated objects)
            bool outputs(const Vec7& q, double guess, Out& out) const
            {
                if (!project(q, guess, out.eta, &out.iterations)) return false;
                Vec7 Q, Q1, Q2;
                c_.eval(out.eta, Q, &Q1, &Q2);
                Mat76 N, dN;
                N_.eval(out.eta, N, &dN);
                const Vec7 e = q - Q;
                out.denom = Q1.squaredNorm() - e.dot(Q2);
                const Eigen::Matrix<double, 1, 7> deta = Q1.transpose() / out.denom;
                out.xi = N.transpose() * e;
                out.A.row(0) = deta;
                out.A.bottomRows<6>() = N.transpose() + (dN.transpose() * e - N.transpose() * Q1) * deta;
                return out.xi.allFinite() && out.A.allFinite() && out.denom > 0.0;
            }

        private:
            Mat76 rate(double s, const Mat76& N) const
            {
                Vec7 Q, Q1, Q2;
                c_.eval(s, Q, &Q1, &Q2);
                const double v = Q1.norm();
                const Vec7 t = Q1 / v;
                const Vec7 tp = (Q2 - t * t.dot(Q2)) / v;
                return -t * (tp.transpose() * N);
            }

            double L_;
            int n_;
            PeriodicQuintic<Vec7> c_;
            PeriodicHermite<Mat76> N_;
            Mat6 Phi_ = Mat6::Zero();
            JointTubeDiagnostics diag_;
    };

    // Controller J: same structure, laws and gains as Controller<PandaArm>, on (eta_J, xi_J) in place of
    // (eta, xi, beta). qdd* = A_J^-1 (mu - Adot_J qdot).
    class JointController {
        public:
            struct Status {
                double eta = 0, etaDot = 0, eTan = 0, sigmaMin = 0, denom = 0;
                Vec6 xi = Vec6::Zero(), xiDot = Vec6::Zero();
                int iterations = 0;
                bool ok = false;
            };

            JointController(const JointTube& tube, const Gains& g, double nu) : tube_(tube), g_(g), nu_(nu) {}
            void reset(double eta0) { etaHat_ = tube_.wrap(eta0); }

            bool step(const Vec7& q, const Vec7& dq, double /*t*/, double dt, Vec7& qdd, Status& st)
            {
                st = Status();
                JointTube::Out out;
                if (!tube_.outputs(q, etaHat_, out)) return false;
                const Vec7 yDot = out.A * dq;
                const Vec6 xiDot = yDot.tail<6>();

                Vec7 Adq = Vec7::Zero();
                if (!(dq.array() == 0.0).all()) {
                    constexpr double h = 1e-5;
                    JointTube::Out plus, minus;
                    if (!tube_.outputs(q + h * dq, out.eta, plus) || !tube_.outputs(q - h * dq, out.eta, minus))
                        return false;
                    Adq = ((plus.A - minus.A) * dq) / (2.0 * h);
                }

                Vec7 mu;
                st.eTan = yDot[0] - nu_;
                mu[0] = -g_.kt * st.eTan;                           // dnu/dt = 0, dnu/deta = 0
                mu.tail<6>() = -g_.Kp * out.xi - g_.Kd * xiDot;     // all six joint-space transverse coordinates

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
            const JointTube& tube_;
            Gains g_;
            double nu_;
            double etaHat_ = 0;
    };

}  // namespace cc
