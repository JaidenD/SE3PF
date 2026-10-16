// Fair comparison: paper's controller P (SE(3) normal tubular coordinates + redundancy output beta) against
// controller J (joint-space tubular coordinates around the lifted curve q*(s), Euclidean, rotation-minimizing frame;
// joint_tube.hpp) and controller L (joint-space coordinates around q*(s) linearized from the paper's outputs at the
// lift, so the same gains act in task units; joint_lin.hpp). L was added by the review as the strongest version of J.
// Build (rev/include is a snapshot of SE3_pathFollowing/include with the spline.hpp the library was built from;
// rev/lib/libse3pf_rev.a is the library rebuilt from the current sources with -O3 -DNDEBUG -DEIGEN_RUNTIME_NO_MALLOC):
//   g++ -O3 -DNDEBUG -std=c++17 -DEIGEN_RUNTIME_NO_MALLOC -I rev/include -I /usr/include/eigen3 -I . compare.cpp
//       rev/lib/libse3pf_rev.a -o compare_rel && ./compare_rel out_rel > out_rel/log_rel.txt
// (-O2 without -DNDEBUG gives compare_o2, which also forbids heap allocation around every controller step.)
#include "joint_tube.hpp"
#include "joint_lin.hpp"
#include "se3pf/panda_example.hpp"
#include "se3pf/se3.hpp"
#include "se3pf/lift.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

using namespace cc;

namespace {
    constexpr double dt = 1e-3;
    constexpr double nuD = PandaExample::nuD;
    std::string outDir = ".";

#ifndef NDEBUG
    struct NoMalloc {
        NoMalloc() { Eigen::internal::set_is_malloc_allowed(false); }
        ~NoMalloc() { Eigen::internal::set_is_malloc_allowed(true); }
    };
#else
    struct NoMalloc {};
#endif

    Vec6 poseError(const Pose& g, const Pose& target)
    {
        const Pose rel = inverse(g) * target;
        Vec6 e;
        e << logSO3(Matrix3d(rel.block<3, 3>(0, 0))), rel.block<3, 1>(0, 3);
        return e;
    }

    // The example's lift nodes, polished by up to 4 more projected Newton updates (the same update as Lift's
    // continuation) so that the fourth-order second differences used for q*'' are not dominated by the 1e-12
    // stopping tolerance. maxChange receives the largest change of a node (rad).
    std::vector<Vec7> polishedLiftNodes(const PandaExample& ex, double& maxChange)
    {
        const Lift& lift = ex.lift();
        std::vector<Vec7> Q = lift.nodes();
        const int n = lift.gridIntervals();
        const double h = lift.length() / n;
        const Vec7 qr = lift.restPosture();
        maxChange = 0;
        for (int k = 0; k < n; ++k) {
            const Pose target = ex.path().pose(k * h);
            Vec7 q = Q[k];
            for (int it = 0; it < 4; ++it) {
                Mat67 J;
                const Pose g = PandaArm::output(q, &J);
                const Vec6 e = poseError(g, target);
                const Eigen::LLT<Mat6> JJt(J * J.transpose());
                const Vec7 d = qr - q;
                const Vec7 nullStep = d - J.transpose() * JJt.solve(Vec6(J * d));
                q += J.transpose() * JJt.solve(e) + nullStep;
            }
            maxChange = std::max(maxChange, (q - Q[k]).cwiseAbs().maxCoeff());
            Q[k] = q;
        }
        Q[n] = Q[0];
        return Q;
    }

    // Redundancy output for P built from a (wrong) lift given by its nodes: q*_used by cubic Hermite with
    // fourth-order differences (as Lift), X_used = kernelVector(J_b(q*_used)) at the nodes (as KernelFrame).
    class NodeLiftRedundancy final : public RedundancyOutput {
        public:
            NodeLiftRedundancy(std::vector<Vec7> Q, double L)
            {
                Q.back() = Q.front();
                std::vector<Vec7> dQ = PeriodicHermite<Vec7>::centralDifferences(Q, L);
                std::vector<Vec7> X(Q.size());
                for (std::size_t k = 0; k < Q.size(); ++k) {
                    Mat67 J;
                    PandaArm::output(Q[k], &J);
                    X[k] = kernelVector(J);
                }
                std::vector<Vec7> dX = PeriodicHermite<Vec7>::centralDifferences(X, L);
                q_.build(L, Q, dQ);
                X_.build(L, X, dX);
            }
            void evaluate(const Vec7& q, double eta, const Vec7& deta, double& beta, Vec7& dbeta) const override
            {
                Vec7 qs, dqs, X, dX;
                q_.eval(eta, qs, &dqs);
                X_.eval(eta, X, &dX);
                const Vec7 e = q - qs;
                beta = X.dot(e);
                dbeta = X + (dX.dot(e) - X.dot(dqs)) * deta;
            }
            void levelSet(double eta, double value, Vec7& p, double& b) const override
            {
                Vec7 qs;
                q_.eval(eta, qs);
                X_.eval(eta, p);
                b = p.dot(qs) + value;
            }

        private:
            PeriodicHermite<Vec7> q_, X_;
    };

    struct Row {
        double t = 0;
        double xiSE3 = 0, etaSE3 = 0, etaDotSE3 = 0, betaP = 0;  // task outputs (true tube, true lift's beta)
        double posErr = 0, attErrDeg = 0;                       // Pi(q) against gamma(eta)
        double jdist = 0, etaJ = 0;                              // joint-space distance to the true lift
        double ownNorm = 0, ownEtaDot = 0, ownBeta = 0, sigma = 0, denom = 0;
        Vec6 ownVec = Vec6::Zero(), ownVecDot = Vec6::Zero();
        double velFrac = 0, accFrac = 0, margin = 0;
        int iterations = 0;
    };

    struct RunOut {
        std::string label;
        std::vector<Row> rows;
        bool ok = true;
        double failT = -1;
        double maxStepMicro = 0;
    };

    using StepFn = std::function<bool(const Vec7&, const Vec7&, double, Vec7&, Row&)>;
    using KickFn = std::function<void(long, const Vec7&, Vec7&)>;

    RunOut runLoop(const std::string& label, const PandaExample& ex, const JointTube& trueJT, Vec7 q, Vec7 qd,
                   double T, const StepFn& step, const KickFn& kick = nullptr)
    {
        RunOut run;
        run.label = label;
        Controller<PandaArm> meas(ex.tube(), PandaExample::gains(), PandaExample::reference(), &ex.redundancy());
        double gS = 0, gJ = 0;
        if (!ex.tube().globalProject(PandaArm::output(q, nullptr), gS) || !trueJT.globalProject(q, gJ))
            throw std::runtime_error("runLoop: initial projection failed");
        const long steps = stepCount(T, dt, 100000000);
        run.rows.reserve(steps + 1);
        const Vec7 dqMax = PandaArm::dqMax(), ddqMax = PandaArm::ddqMax();
        for (long k = 0; k <= steps; ++k) {
            const double t = k * dt;
            if (kick) kick(k, q, qd);
            Row row;
            row.t = t;
            Controller<PandaArm>::OutputState mo;
            if (!meas.outputs(q, gS, mo)) { run.ok = false; run.failT = t; break; }
            const Vec7 yS = mo.A * qd;
            row.xiSE3 = mo.xi.norm();
            row.etaSE3 = mo.eta;
            row.etaDotSE3 = yS[0];
            row.betaP = mo.beta;
            gS = ex.tube().wrap(mo.eta + yS[0] * dt);
            {
                const Pose g = PandaArm::output(q, nullptr);
                const Pose rel = inverse(ex.path().pose(mo.eta)) * g;
                row.posErr = rel.block<3, 1>(0, 3).norm();
                row.attErrDeg = logSO3(Matrix3d(rel.block<3, 3>(0, 0))).norm() * 180.0 / EIGEN_PI;
            }
            JointTube::Out jo;
            if (!trueJT.outputs(q, gJ, jo)) { run.ok = false; run.failT = t; break; }
            row.jdist = jo.xi.norm();
            row.etaJ = jo.eta;
            gJ = trueJT.wrap(jo.eta + jo.A.row(0).dot(qd) * dt);

            Vec7 qdd = Vec7::Zero();
            bool ok;
            {
                NoMalloc guard;
                (void)guard;
                ok = step(q, qd, t, qdd, row);
            }
            if (!ok) { run.ok = false; run.failT = t; break; }
            row.velFrac = (qd.cwiseAbs().array() / dqMax.array()).maxCoeff();
            row.accFrac = (qdd.cwiseAbs().array() / ddqMax.array()).maxCoeff();
            row.margin = PandaArm::jointMargin(q);
            run.rows.push_back(row);
            q += dt * qd + (0.5 * dt * dt) * qdd;
            qd += dt * qdd;
        }
        return run;
    }

    StepFn stepP(Controller<PandaArm>& c)
    {
        return [&c](const Vec7& q, const Vec7& qd, double t, Vec7& qdd, Row& row) {
            Controller<PandaArm>::Status st;
            if (!c.step(q, qd, t, dt, qdd, st)) return false;
            row.ownNorm = st.xi.norm();
            row.ownEtaDot = st.etaDot;
            row.ownBeta = st.beta;
            row.sigma = st.sigmaMin;
            row.ownVec.head<5>() = st.xi;
            row.ownVec[5] = st.beta;
            row.ownVecDot.head<5>() = st.xiDot;
            row.ownVecDot[5] = st.betaDot;
            row.iterations = st.iterations;
            return true;
        };
    }

    StepFn stepJ(JointController& c)
    {
        return [&c](const Vec7& q, const Vec7& qd, double t, Vec7& qdd, Row& row) {
            JointController::Status st;
            if (!c.step(q, qd, t, dt, qdd, st)) return false;
            row.ownNorm = st.xi.norm();
            row.ownEtaDot = st.etaDot;
            row.sigma = st.sigmaMin;
            row.denom = st.denom;
            row.ownVec = st.xi;
            row.ownVecDot = st.xiDot;
            row.iterations = st.iterations;
            return true;
        };
    }

    StepFn stepL(GenericJointController<LinTube>& c)
    {
        return [&c](const Vec7& q, const Vec7& qd, double t, Vec7& qdd, Row& row) {
            JointController::Status st;
            if (!c.step(q, qd, t, dt, qdd, st)) return false;
            row.ownNorm = st.xi.norm();
            row.ownEtaDot = st.etaDot;
            row.ownBeta = st.xi[5];
            row.sigma = st.sigmaMin;
            row.denom = st.denom;
            row.ownVec = st.xi;
            row.ownVecDot = st.xiDot;
            row.iterations = st.iterations;
            return true;
        };
    }

    // max over rows k0..k1 of |y(t) - y_design(t)| (vector norm), y_design from y(t_k0), y'(t_k0) and the
    // critically damped law y'' = -16 y - 8 y'; dims = number of components compared (6 for J, 5 for P's xi)
    double designDeviation(const RunOut& r, std::size_t k0, std::size_t k1, int dims, double* normDev = nullptr)
    {
        const Vec6 x0 = r.rows[k0].ownVec, v0 = r.rows[k0].ownVecDot;
        double m = 0, mn = 0;
        for (std::size_t k = k0; k <= k1 && k < r.rows.size(); ++k) {
            const double tau = r.rows[k].t - r.rows[k0].t;
            const Vec6 yd = (x0 + (v0 + 4.0 * x0) * tau) * std::exp(-4.0 * tau);
            m = std::max(m, (r.rows[k].ownVec - yd).head(dims).norm());
            mn = std::max(mn, std::abs(r.rows[k].ownVec.head(dims).norm() - yd.head(dims).norm()));
        }
        if (normDev) *normDev = mn;
        return m;
    }

    void writeCsv(const RunOut& r, const std::string& file)
    {
        std::ofstream f(outDir + "/" + file);
        f << "t,xiSE3,posErr,attErrDeg,etaDotSE3,betaP,jdist,ownNorm,ownEtaDot,ownBeta,velFrac,accFrac,margin\n";
        f.precision(10);
        for (std::size_t k = 0; k < r.rows.size(); ++k) {
            const Row& w = r.rows[k];
            const bool nearKick = std::abs(w.t - 5.0) < 0.3;
            if (k % 10 != 0 && !nearKick && k + 1 != r.rows.size()) continue;
            f << w.t << ',' << w.xiSE3 << ',' << w.posErr << ',' << w.attErrDeg << ',' << w.etaDotSE3 << ','
              << w.betaP << ',' << w.jdist << ',' << w.ownNorm << ',' << w.ownEtaDot << ',' << w.ownBeta << ','
              << w.velFrac << ',' << w.accFrac << ',' << w.margin << '\n';
        }
    }

    struct Window { double a, b; };

    template <class F>
    double maxOver(const RunOut& r, Window w, F f, double* at = nullptr)
    {
        double m = -1, tm = 0;
        for (const Row& x : r.rows)
            if (x.t >= w.a - 1e-9 && x.t <= w.b + 1e-9 && f(x) > m) { m = f(x); tm = x.t; }
        if (at) *at = tm;
        return m;
    }

    void summarize(std::ostream& o, const std::string& sc, const RunOut& r, double T, double kickT = -1)
    {
        const Row& f = r.rows.back();
        const Row& s = r.rows.front();
        double tMax = 0;
        o << "[" << sc << " " << r.label << "] ok=" << r.ok << " failT=" << r.failT << " rows=" << r.rows.size() << "\n";
        o << "  tool |xi_SE3|: start=" << s.xiSE3 << " max=" << maxOver(r, {0, T}, [](const Row& x) { return x.xiSE3; }, &tMax)
          << " @t=" << tMax << " final=" << f.xiSE3
          << " steady(max last 5 s)=" << maxOver(r, {T - 5, T}, [](const Row& x) { return x.xiSE3; }) << "\n";
        o << "  tool pos err: max=" << maxOver(r, {0, T}, [](const Row& x) { return x.posErr; })
          << " m, steady=" << maxOver(r, {T - 5, T}, [](const Row& x) { return x.posErr; })
          << " m; att err: max=" << maxOver(r, {0, T}, [](const Row& x) { return x.attErrDeg; })
          << " deg, steady=" << maxOver(r, {T - 5, T}, [](const Row& x) { return x.attErrDeg; }) << " deg\n";
        if (kickT > 0) {
            o << "  pre-kick max |xi_SE3| (1..5 s)=" << maxOver(r, {1, kickT - 1e-3}, [](const Row& x) { return x.xiSE3; })
              << "  post-kick max=" << maxOver(r, {kickT, T}, [](const Row& x) { return x.xiSE3; }, &tMax) << " @t=" << tMax
              << "  post-kick max posErr=" << maxOver(r, {kickT, T}, [](const Row& x) { return x.posErr; })
              << "  max attErrDeg=" << maxOver(r, {kickT, T}, [](const Row& x) { return x.attErrDeg; }) << "\n";
        }
        o << "  tool eta_dot: final=" << f.etaDotSE3 << " max|eta_dot - nu_d| last 5 s="
          << maxOver(r, {T - 5, T}, [](const Row& x) { return std::abs(x.etaDotSE3 - nuD); })
          << " max|eta_dot - nu_d| after 2 s=" << maxOver(r, {2, T}, [](const Row& x) { return std::abs(x.etaDotSE3 - nuD); }) << "\n";
        o << "  joint dist to lift: start=" << s.jdist << " max=" << maxOver(r, {0, T}, [](const Row& x) { return x.jdist; })
          << " final=" << f.jdist << "  paper beta(true lift): start=" << s.betaP << " final=" << f.betaP << "\n";
        o << "  peak |qd|/qdmax=" << maxOver(r, {0, T}, [](const Row& x) { return x.velFrac; })
          << " peak |qdd*|/qddmax=" << maxOver(r, {0, T}, [](const Row& x) { return x.accFrac; }, &tMax) << " @t=" << tMax
          << " min joint margin=" << -maxOver(r, {0, T}, [](const Row& x) { return -x.margin; }) << " rad\n";
        o << "  own outputs: final |own xi|=" << f.ownNorm << " own eta_dot=" << f.ownEtaDot << " own beta=" << f.ownBeta
          << " min sigma(A)=" << -maxOver(r, {0, T}, [](const Row& x) { return -x.sigma; })
          << " max iterations=" << maxOver(r, {0, T}, [](const Row& x) { return double(x.iterations); }) << "\n";
    }

    template <class Ctrl, class St>
    double medianStepMicro(Ctrl& c, Vec7 q, Vec7 qd, int N)
    {
        std::vector<double> times;
        times.reserve(N);
        for (int k = 0; k < N + 200; ++k) {
            Vec7 qdd;
            St st;
            const auto t0 = std::chrono::steady_clock::now();
            const bool ok = c.step(q, qd, k * dt, dt, qdd, st);
            const auto t1 = std::chrono::steady_clock::now();
            if (!ok) throw std::runtime_error("timing: step failed");
            if (k >= 200) times.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
            q += dt * qd + (0.5 * dt * dt) * qdd;
            qd += dt * qdd;
        }
        std::nth_element(times.begin(), times.begin() + times.size() / 2, times.end());
        return times[times.size() / 2];
    }
}  // namespace

int main(int argc, char** argv)
{
    if (argc > 1) outDir = argv[1];
    const bool timingOnly = argc > 2 && std::string(argv[2]) == "timing";
    std::ostream& o = std::cout;
    o.precision(6);

    const auto tb0 = std::chrono::steady_clock::now();
    PandaExample ex;
    const double L = ex.tube().length();
    double polish = 0;
    const std::vector<Vec7> Qn = polishedLiftNodes(ex, polish);
    JointTube jt(Qn, L);
    const int n = ex.lift().gridIntervals();
    const Gains G = PandaExample::gains();

    // Wrong lift for S4
    std::vector<Vec7> Qw = Qn;
    for (int k = 0; k <= n; ++k) {
        const double s = k * L / n;
        Qw[k] += 0.02 * (std::sin(2 * EIGEN_PI * s / L) * Vec7::Unit(1) + std::cos(2 * EIGEN_PI * s / L) * Vec7::Unit(3));
    }
    Qw[n] = Qw[0];
    JointTube jtWrong(Qw, L);
    NodeLiftRedundancy betaWrong(Qw, L);

    // Variant L: the paper's A at the lift nodes (offline), rows (d eta; d xi; d beta)
    auto liftA = [&](const std::vector<Vec7>& Q, const RedundancyOutput& red) {
        Controller<PandaArm> c(ex.tube(), G, PandaExample::reference(), &red);
        std::vector<Mat7> A(n + 1);
        for (int k = 0; k < n; ++k) {
            Controller<PandaArm>::OutputState out;
            if (!c.outputs(Q[k], k * L / n, out)) throw std::runtime_error("liftA: outputs failed");
            A[k] = out.A;
        }
        A[n] = A[0];
        return A;
    };
    LinTube lt(Qn, liftA(Qn, ex.redundancy()), L);
    LinTube ltWrong(Qw, liftA(Qw, betaWrong), L);
    const auto tb1 = std::chrono::steady_clock::now();
    o << "Build (example + polished lift + joint tube J + linearized tubes L): "
      << std::chrono::duration<double>(tb1 - tb0).count() << " s; L = " << L << ", grid n = " << jt.gridIntervals() << "\n";

    Vec7 qs01, dqs01;
    ex.lift().eval(PandaExample::s0, qs01, &dqs01);

    // ---------------- S5 timing ----------------
    {
        Controller<PandaArm> P(ex.tube(), G, PandaExample::reference(), &ex.redundancy());
        JointController J(jt, G, nuD);
        GenericJointController<LinTube> Lc(lt, G, nuD);
        double e0 = 0;
        P.reset(PandaExample::s0);
        const double tP = medianStepMicro<Controller<PandaArm>, Controller<PandaArm>::Status>(P, ex.paperStart(), Vec7::Zero(), 5000);
        jt.globalProject(ex.paperStart(), e0);
        J.reset(e0);
        const double tJ = medianStepMicro<JointController, JointController::Status>(J, ex.paperStart(), Vec7::Zero(), 5000);
        lt.globalProject(ex.paperStart(), e0);
        Lc.reset(e0);
        const double tL = medianStepMicro<GenericJointController<LinTube>, JointController::Status>(Lc, ex.paperStart(), Vec7::Zero(), 5000);
        P.reset(PandaExample::s0);
        const double tP3 = medianStepMicro<Controller<PandaArm>, Controller<PandaArm>::Status>(P, qs01, nuD * dqs01, 5000);
        jt.globalProject(qs01, e0);
        J.reset(e0);
        const double tJ3 = medianStepMicro<JointController, JointController::Status>(J, qs01, nuD * dqs01, 5000);
        lt.globalProject(qs01, e0);
        Lc.reset(e0);
        const double tL3 = medianStepMicro<GenericJointController<LinTube>, JointController::Status>(Lc, qs01, nuD * dqs01, 5000);
        o << "S5 timing (median over 5000 steps, us): P=" << tP << " J=" << tJ << " L=" << tL << " (S1 start); P=" << tP3
          << " J=" << tJ3 << " L=" << tL3 << " (S3 start)\n";
        if (timingOnly) return 0;
    }

    // ---------------- Validation of J ----------------
    {
        const JointTubeDiagnostics& d = jt.diagnostics();
        o << "JointTube: polish max node change=" << polish << " rad; closureGap=" << d.closureGap
          << " holonomyAngle=" << d.holonomyAngle << " rad detH=" << d.detH << " orthH=" << d.orthH
          << " logError=" << d.logError << "\n  max|N^T N - I| nodes=" << d.maxOrthoNodes << " interp=" << d.maxOrthoInterp
          << "  max|N^T t| nodes=" << d.maxNormalNodes << " interp=" << d.maxNormalInterp
          << "\n  |q*'| in [" << d.minSpeed << ", " << d.maxSpeed << "] rad per unit s; max curvature="
          << d.maxCurvature << " 1/rad (focal distance " << 1.0 / d.maxCurvature << " rad)\n";
        const JointTubeDiagnostics& dw = jtWrong.diagnostics();
        o << "JointTube (wrong lift): closureGap=" << dw.closureGap << " holonomyAngle=" << dw.holonomyAngle
          << " max|N^T N - I| interp=" << dw.maxOrthoInterp << " max|N^T t| interp=" << dw.maxNormalInterp << "\n";

        double dQ = 0, dQ1 = 0;
        for (int k = 0; k < 4 * n; ++k) {
            const double s = (k + 0.5) * L / (4.0 * n);
            Vec7 a, a1, b, b1;
            jt.curve(s, a, &a1);
            ex.lift().eval(s, b, &b1);
            dQ = std::max(dQ, (a - b).cwiseAbs().maxCoeff());
            dQ1 = std::max(dQ1, (a1 - b1).cwiseAbs().maxCoeff());
        }
        o << "quintic q* vs Lift::eval: max|dq*|=" << dQ << " rad, max|dq*'|=" << dQ1 << "\n";

        // A_J and A_L against central differences
        auto fdCheck = [&](auto& tube, const char* name) {
            std::mt19937 rng(7);
            std::normal_distribution<double> nd(0, 1);
            double maxErr = 0, maxRel = 0;
            int fails = 0;
            for (int trial = 0; trial < 20; ++trial) {
                const double s = L * (trial + 0.3) / 20.0;
                Vec7 Q;
                tube.curve(s, Q);
                Vec7 r;
                for (int i = 0; i < 7; ++i) r[i] = nd(rng);
                const Vec7 q = Q + (0.05 + 0.15 * trial / 19.0) * r.normalized();
                typename std::decay_t<decltype(tube)>::Out out;
                if (!tube.outputs(q, s, out)) { ++fails; continue; }
                Mat7 Afd;
                const double h = 1e-6;
                for (int j = 0; j < 7; ++j) {
                    typename std::decay_t<decltype(tube)>::Out p, m;
                    tube.outputs(q + h * Vec7::Unit(j), out.eta, p);
                    tube.outputs(q - h * Vec7::Unit(j), out.eta, m);
                    Afd(0, j) = wrapSymmetric(p.eta - m.eta, L) / (2 * h);
                    Afd.block<6, 1>(1, j) = (p.xi - m.xi) / (2 * h);
                }
                maxErr = std::max(maxErr, (Afd - out.A).cwiseAbs().maxCoeff());
                maxRel = std::max(maxRel, (Afd - out.A).cwiseAbs().maxCoeff() / out.A.cwiseAbs().maxCoeff());
            }
            o << name << " vs central differences (h=1e-6, 20 points, |q - q*| 0.05..0.2 rad): max abs err=" << maxErr
              << " rel=" << maxRel << " projection failures=" << fails << "\n";
        };
        fdCheck(jt, "A_J");
        fdCheck(lt, "A_L");

        const LinTubeDiagnostics& dl = lt.diagnostics();
        o << "LinTube: max|a^T q*' - 1|=" << dl.maxATQ1 << " max|B q*'|=" << dl.maxBQ1 << " min sigma(A(q*))=" << dl.minSigmaA
          << "; wrong lift: max|a^T q*' - 1|=" << ltWrong.diagnostics().maxATQ1 << " max|B q*'|="
          << ltWrong.diagnostics().maxBQ1 << "\n";

        // L's coordinates against the paper's at distance d from the lift: the difference must be O(d^2)
        Controller<PandaArm> Pm(ex.tube(), G, PandaExample::reference(), &ex.redundancy());
        for (double dd : {0.01, 0.03, 0.1, 0.3}) {
            std::mt19937 rng(3);
            std::normal_distribution<double> nd(0, 1);
            double mz = 0, me = 0, mzJ = 0;
            for (int trial = 0; trial < 40; ++trial) {
                const double s = L * (trial + 0.6) / 40.0;
                Vec7 Q;
                lt.curve(s, Q);
                Vec7 r;
                for (int i = 0; i < 7; ++i) r[i] = nd(rng);
                const Vec7 q = Q + dd * r.normalized();
                Controller<PandaArm>::OutputState po;
                LinTube::Out lo;
                const bool okP = Pm.outputs(q, s, po), okL = lt.outputs(q, s, lo);
                if (!okP || !okL) {
                    o << "  O(d^2) check at d = " << dd << ", s = " << s << ": projection failed (P " << okP << ", L " << okL
                      << "); skipped\n";
                    continue;
                }
                Vec6 zP;
                zP << po.xi, po.beta;
                mz = std::max(mz, (lo.xi - zP).norm());
                me = std::max(me, std::abs(wrapSymmetric(lo.eta - po.eta, L)));
            }
            o << "  L vs paper outputs at |q - q*| = " << dd << ": max|z - (xi, beta)|=" << mz << " (/d^2 = " << mz / (dd * dd)
              << "), max|eta_L - eta|=" << me << " (/d^2 = " << me / (dd * dd) << ")\n";
            (void)mzJ;
        }
    }

    {
        const auto t0 = std::chrono::steady_clock::now();
        const double eps = estimateTubeRadius(ex.tube());
        o << "SE(3) tube radius estimate (Panda path, W = diag(0.01 I, I)): " << eps << " ("
          << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s)\n";
    }

    auto report = [&](const std::string& sc, RunOut& r, double T, double kickT = -1) {
        summarize(o, sc, r, T, kickT);
        writeCsv(r, sc + "_" + r.label + ".csv");
    };

    struct Trio { RunOut P, J, L; };
    // The three controllers from the same q0, qd0 with the same loop. P with redundancy output redP; J and L on the
    // given tubes; each starts its projection at its own closest point.
    auto runTrio = [&](const Vec7& q0, const Vec7& qd0, double T, const KickFn& kick, const RedundancyOutput& redP,
                       const JointTube& jtUse, const LinTube& ltUse) {
        Trio tr;
        Controller<PandaArm> P(ex.tube(), G, PandaExample::reference(), &redP);
        P.reset(PandaExample::s0);
        tr.P = runLoop("P", ex, jt, q0, qd0, T, stepP(P), kick);
        JointController J(jtUse, G, nuD);
        double e0 = 0;
        if (!jtUse.globalProject(q0, e0)) throw std::runtime_error("runTrio: J projection");
        J.reset(e0);
        tr.J = runLoop("J", ex, jt, q0, qd0, T, stepJ(J), kick);
        GenericJointController<LinTube> Lc(ltUse, G, nuD);
        if (!ltUse.globalProject(q0, e0)) throw std::runtime_error("runTrio: L projection");
        Lc.reset(e0);
        tr.L = runLoop("L", ex, jt, q0, qd0, T, stepL(Lc), kick);
        return tr;
    };
    auto peak = [](const RunOut& r, double a, double b, double* at = nullptr) {
        return maxOver(r, {a, b}, [](const Row& x) { return x.xiSE3; }, at);
    };

    // ---------------- S1: paper start ----------------
    {
        const Vec7 q0 = ex.paperStart();
        Trio tr = runTrio(q0, Vec7::Zero(), 20.0, nullptr, ex.redundancy(), jt, lt);
        report("S1", tr.P, 20.0);
        report("S1", tr.J, 20.0);
        report("S1", tr.L, 20.0);
        double nd = 0;
        const double devP = designDeviation(tr.P, 0, tr.P.rows.size() - 1, 5, &nd);
        o << "  VALID S1 P: max|xi - design| (vector)=" << devP << " max||xi| - |xi0|(1+4t)e^-4t|=" << nd << "\n";
        const double devJ = designDeviation(tr.J, 0, tr.J.rows.size() - 1, 6, &nd);
        o << "  VALID S1 J: |xi_J0|=" << tr.J.rows[0].ownNorm << " max|xi_J - design| (vector)=" << devJ
          << " max||xi_J| - |xi_J0|(1+4t)e^-4t|=" << nd << " final eta_J_dot=" << tr.J.rows.back().ownEtaDot << "\n";
        const double devL = designDeviation(tr.L, 0, tr.L.rows.size() - 1, 6, &nd);
        o << "  VALID S1 L: |z0|=" << tr.L.rows[0].ownNorm << " (paper |(xi0, beta0)| = " << tr.P.rows[0].ownVec.norm()
          << ") max|z - design| (vector)=" << devL << " max||z| - |z0|(1+4t)e^-4t|=" << nd
          << " final eta_L_dot=" << tr.L.rows.back().ownEtaDot << "\n";
        auto ratio = [&](const RunOut& r, double& lo, double& hi) {
            lo = 1e9; hi = 0;
            for (const Row& w : r.rows)
                if (w.t <= 1.5) {
                    const double des = PandaExample::designedResponse(r.rows[0].xiSE3, w.t);
                    lo = std::min(lo, w.xiSE3 / des);
                    hi = std::max(hi, w.xiSE3 / des);
                }
        };
        auto devTool = [&](const RunOut& r) {
            double m = 0;
            for (const Row& w : r.rows) m = std::max(m, std::abs(w.xiSE3 - PandaExample::designedResponse(r.rows[0].xiSE3, w.t)));
            return m;
        };
        for (const RunOut* r : {&tr.P, &tr.J, &tr.L}) {
            double lo = 0, hi = 0;
            ratio(*r, lo, hi);
            o << "  S1 " << r->label << " tool error / designed response over 0..1.5 s in [" << lo << ", " << hi
              << "]; max | |xi_SE3| - design | = " << devTool(*r) << "; |xi_SE3| at t = 0.5, 1, 2 s: "
              << r->rows[500].xiSE3 << ", " << r->rows[1000].xiSE3 << ", " << r->rows[2000].xiSE3 << "\n";
        }
        o << "  S1 design at t = 0.5, 1, 2 s: " << PandaExample::designedResponse(tr.P.rows[0].xiSE3, 0.5) << ", "
          << PandaExample::designedResponse(tr.P.rows[0].xiSE3, 1.0) << ", "
          << PandaExample::designedResponse(tr.P.rows[0].xiSE3, 2.0) << "\n";
        {
            double e0 = 0, eL = 0;
            jt.globalProject(q0, e0);
            lt.globalProject(q0, eL);
            Vec7 Q, Q1, Q2;
            jt.curve(e0, Q, &Q1, &Q2);
            const Vec7 t = Q1.normalized();
            const double kappa = ((Q2 - t * t.dot(Q2)) / Q1.norm()).norm() / Q1.norm();
            o << "  S1 J: eta_J0=" << e0 << ", L: eta_L0=" << eL << " (eta0=" << PandaExample::s0
              << "), J local focal distance 1/kappa=" << 1.0 / kappa << " rad, start distance " << tr.J.rows[0].jdist
              << " rad, denom0 J=" << tr.J.rows[0].denom << " L=" << tr.L.rows[0].denom << "\n";
        }
    }

    // ---------------- S2: elbow off the lift, tool on the path ----------------
    auto qd0P = [&](const Vec7& q0) {
        Controller<PandaArm> P(ex.tube(), G, PandaExample::reference(), &ex.redundancy());
        Controller<PandaArm>::OutputState out;
        if (!P.outputs(q0, PandaExample::s0, out)) throw std::runtime_error("S2: outputs");
        return Vec7(out.A.partialPivLu().solve(Vec7(nuD * Vec7::Unit(0))));
    };
    auto runS2 = [&](const std::string& sc, const Vec7& q0, bool jNative = false) {
        Vec7 qd0 = qd0P(q0);
        if (jNative) {  // J's own "moving along, transverse at rest": A_J(q0) qd0 = (nu_d, 0, ..., 0)
            double e0 = 0;
            jt.globalProject(q0, e0);
            JointTube::Out jo;
            jt.outputs(q0, e0, jo);
            qd0 = jo.A.partialPivLu().solve(Vec7(nuD * Vec7::Unit(0)));
        }
        Trio tr = runTrio(q0, qd0, 20.0, nullptr, ex.redundancy(), jt, lt);
        report(sc, tr.P, 20.0);
        report(sc, tr.J, 20.0);
        report(sc, tr.L, 20.0);
        double nd = 0;
        const double beta0 = tr.P.rows[0].ownBeta;
        double bdev = 0;
        for (const Row& w : tr.P.rows) bdev = std::max(bdev, std::abs(w.ownBeta - PandaExample::designedResponse(beta0, w.t)));
        o << "  VALID " << sc << " P: beta0=" << beta0 << " max|beta - beta0(1+4t)e^-4t|=" << bdev
          << " max|xi|=" << maxOver(tr.P, {0, 20}, [](const Row& x) { return x.ownNorm; }) << "\n";
        const double devJ = designDeviation(tr.J, 0, tr.J.rows.size() - 1, 6, &nd);
        o << "  VALID " << sc << " J: |xi_J0|=" << tr.J.rows[0].ownNorm << " |xi_J_dot0|=" << tr.J.rows[0].ownVecDot.norm()
          << " max|xi_J - design| (vector)=" << devJ << " final eta_J_dot=" << tr.J.rows.back().ownEtaDot << "\n";
        const double devL = designDeviation(tr.L, 0, tr.L.rows.size() - 1, 6, &nd);
        o << "  VALID " << sc << " L: z0=(|xi part| " << tr.L.rows[0].ownVec.head<5>().norm() << ", beta part "
          << tr.L.rows[0].ownVec[5] << ") |z_dot0|=" << tr.L.rows[0].ownVecDot.norm()
          << " max|z - design| (vector)=" << devL << " final eta_L_dot=" << tr.L.rows.back().ownEtaDot << "\n";
    };
    runS2("S2", ex.selfMotionStart());
    runS2("S2small", ex.start(ex.redundancy(), Vec5::Zero(), 0.04));
    runS2("S2jnative", ex.selfMotionStart(), true);
    for (double b0 : {0.4, 0.6, 0.8}) {
        const Vec7 q0 = ex.start(ex.redundancy(), Vec5::Zero(), b0);
        o << "S2large beta0=" << b0 << ": start joint margin " << PandaArm::jointMargin(q0) << "\n";
        runS2("S2large" + std::to_string(int(b0 * 10)), q0);
    }

    // S2 sweep: peak tool error against the elbow offset beta0 (4 s runs; the peak is before 1 s)
    {
        std::ofstream f(outDir + "/sweep.csv");
        f << "beta0,P,J,L,P_mm,J_mm,L_mm,P_deg,J_deg,L_deg,P_acc,J_acc,L_acc\n";
        f.precision(10);
        o << "S2 sweep (peak tool error W-norm; mm; deg; peak qdd fraction):\n";
        for (double b0 : {0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8}) {
            const Vec7 q0 = ex.start(ex.redundancy(), Vec5::Zero(), b0);
            Trio tr = runTrio(q0, qd0P(q0), 4.0, nullptr, ex.redundancy(), jt, lt);
            f << b0;
            const RunOut* rs[3] = {&tr.P, &tr.J, &tr.L};
            for (const RunOut* r : rs) f << ',' << peak(*r, 0, 4);
            for (const RunOut* r : rs) f << ',' << 1e3 * maxOver(*r, {0, 4}, [](const Row& x) { return x.posErr; });
            for (const RunOut* r : rs) f << ',' << maxOver(*r, {0, 4}, [](const Row& x) { return x.attErrDeg; });
            for (const RunOut* r : rs) f << ',' << maxOver(*r, {0, 4}, [](const Row& x) { return x.accFrac; });
            f << '\n';
            o << "  beta0=" << b0 << ": P " << peak(tr.P, 0, 4) << "  J " << peak(tr.J, 0, 4) << " ("
              << 1e3 * maxOver(tr.J, {0, 4}, [](const Row& x) { return x.posErr; }) << " mm, "
              << maxOver(tr.J, {0, 4}, [](const Row& x) { return x.attErrDeg; }) << " deg)  L " << peak(tr.L, 0, 4) << " ("
              << 1e3 * maxOver(tr.L, {0, 4}, [](const Row& x) { return x.posErr; }) << " mm, "
              << maxOver(tr.L, {0, 4}, [](const Row& x) { return x.attErrDeg; }) << " deg)  ok " << tr.P.ok << tr.J.ok
              << tr.L.ok << "  acc P/J/L " << maxOver(tr.P, {0, 4}, [](const Row& x) { return x.accFrac; }) << "/"
              << maxOver(tr.J, {0, 4}, [](const Row& x) { return x.accFrac; }) << "/"
              << maxOver(tr.L, {0, 4}, [](const Row& x) { return x.accFrac; }) << "\n";
        }
    }

    // ---------------- S3: elbow push ----------------
    {
        const Vec7 qd0 = nuD * dqs01;
        auto kick = [](long k, const Vec7& q, Vec7& qd) {
            if (k == 5000) {
                Mat67 J;
                PandaArm::output(q, &J);
                qd += 0.3 * kernelVector(J);
            }
        };
        Trio tr = runTrio(qs01, qd0, 15.0, kick, ex.redundancy(), jt, lt);
        report("S3", tr.P, 15.0, 5.0);
        report("S3", tr.J, 15.0, 5.0);
        report("S3", tr.L, 15.0, 5.0);
        double nd = 0;
        const double devJ = designDeviation(tr.J, 5000, tr.J.rows.size() - 1, 6, &nd);
        o << "  VALID S3 J: before kick max|xi_J|=" << maxOver(tr.J, {0, 4.999}, [](const Row& x) { return x.ownNorm; })
          << "; after kick xi_J_dot=" << tr.J.rows[5000].ownVecDot.norm() << " max|xi_J - design| (vector)=" << devJ << "\n";
        const double devL = designDeviation(tr.L, 5000, tr.L.rows.size() - 1, 6, &nd);
        o << "  VALID S3 L: before kick max|z|=" << maxOver(tr.L, {0, 4.999}, [](const Row& x) { return x.ownNorm; })
          << "; after kick z_dot=(xi part " << tr.L.rows[5000].ownVecDot.head<5>().norm() << ", beta part "
          << tr.L.rows[5000].ownVecDot[5] << ") max|z - design| (vector)=" << devL << "\n";
        const double beta5 = tr.P.rows[5000].ownBeta, bd5 = tr.P.rows[5000].ownVecDot[5];
        double bdev = 0;
        for (std::size_t k = 5000; k < tr.P.rows.size(); ++k) {
            const double tau = tr.P.rows[k].t - 5.0;
            bdev = std::max(bdev, std::abs(tr.P.rows[k].ownBeta - (beta5 + (bd5 + 4 * beta5) * tau) * std::exp(-4 * tau)));
        }
        o << "  VALID S3 P: beta_dot after kick=" << bd5 << " max|beta - design|=" << bdev
          << " xi_dot after kick=" << tr.P.rows[5000].ownVecDot.head<5>().norm() << "\n";
    }

    // ---------------- S3r: generic pushes (random joint-velocity directions) ----------------
    {
        std::ofstream f(outDir + "/S3r.csv");
        f << "i,kernelComp,P,J,L,relJ,relL,devP\n";
        f.precision(10);
        std::mt19937 rng(11);
        std::normal_distribution<double> nd(0, 1);
        const Vec7 qd0 = nuD * dqs01;
        std::vector<double> rJ, rL, dJ, dL, pP;
        double maxDevP = 0;
        for (int i = 0; i < 20; ++i) {
            Vec7 dir;
            for (int j = 0; j < 7; ++j) dir[j] = nd(rng);
            dir.normalize();
            double kernelComp = 0;
            auto kick = [&](long k, const Vec7& q, Vec7& qd) {
                if (k == 1000) {
                    Mat67 Jb;
                    PandaArm::output(q, &Jb);
                    kernelComp = std::abs(kernelVector(Jb).dot(dir));
                    qd += 0.3 * dir;
                }
            };
            Trio tr = runTrio(qs01, qd0, 4.0, kick, ex.redundancy(), jt, lt);
            const double p = peak(tr.P, 1, 4), j = peak(tr.J, 1, 4), l = peak(tr.L, 1, 4);
            // trajectory difference of the tool error from P's (which follows the designed response)
            double mj = 0, ml = 0;
            for (std::size_t k = 1000; k < tr.P.rows.size(); ++k) {
                mj = std::max(mj, std::abs(tr.J.rows[k].xiSE3 - tr.P.rows[k].xiSE3));
                ml = std::max(ml, std::abs(tr.L.rows[k].xiSE3 - tr.P.rows[k].xiSE3));
            }
            double ndv = 0;
            const double devP = designDeviation(tr.P, 1000, tr.P.rows.size() - 1, 5, &ndv);
            maxDevP = std::max(maxDevP, devP / p);
            f << i << ',' << kernelComp << ',' << p << ',' << j << ',' << l << ',' << mj / p << ',' << ml / p << ','
              << devP << '\n';
            rJ.push_back(j / p);
            rL.push_back(l / p);
            dJ.push_back(mj / p);
            dL.push_back(ml / p);
            pP.push_back(p);
            if (!tr.P.ok || !tr.J.ok || !tr.L.ok) o << "  S3r " << i << ": a run failed\n";
        }
        auto med = [](std::vector<double> v) {
            std::sort(v.begin(), v.end());
            return v[v.size() / 2];
        };
        auto mx = [](const std::vector<double>& v) { return *std::max_element(v.begin(), v.end()); };
        auto mn = [](const std::vector<double>& v) { return *std::min_element(v.begin(), v.end()); };
        o << "S3r (20 random joint-velocity pushes of 0.3 rad/s at t = 1 s): P peak tool error median " << med(pP)
          << " [" << mn(pP) << ", " << mx(pP) << "]; P max|xi - design| / peak = " << maxDevP << "\n"
          << "  peak ratio J/P median " << med(rJ) << " [" << mn(rJ) << ", " << mx(rJ) << "]; L/P median " << med(rL)
          << " [" << mn(rL) << ", " << mx(rL) << "]\n"
          << "  max_t | |xi_SE3| - |xi_SE3,P| | / P peak: J median " << med(dJ) << " max " << mx(dJ) << "; L median "
          << med(dL) << " max " << mx(dL) << "\n";
    }

    // ---------------- S4: lift error ----------------
    {
        const Vec7 qd0 = nuD * dqs01;
        Trio tr = runTrio(qs01, qd0, 20.0, nullptr, betaWrong, jtWrong, ltWrong);
        report("S4", tr.P, 20.0);
        report("S4", tr.J, 20.0);
        report("S4", tr.L, 20.0);
        double nd = 0;
        const double devJ = designDeviation(tr.J, 0, tr.J.rows.size() - 1, 6, &nd);
        o << "  VALID S4 J (own wrong tube): |xi_J0|=" << tr.J.rows[0].ownNorm << " max|xi_J - design|=" << devJ
          << " final |xi_J|=" << tr.J.rows.back().ownNorm << "\n";
        const double devL = designDeviation(tr.L, 0, tr.L.rows.size() - 1, 6, &nd);
        o << "  VALID S4 L (own wrong chart): |z0|=" << tr.L.rows[0].ownNorm << " max|z - design|=" << devL
          << " final |z|=" << tr.L.rows.back().ownNorm << "\n";
        double fkErr = 0;
        for (int k = 0; k < n; ++k) {
            const Pose g = PandaArm::output(Qw[k], nullptr);
            double eta = 0;
            Vec5 xi;
            ex.tube().outputs(g, k * L / n, eta, xi);
            fkErr = std::max(fkErr, xi.norm());
        }
        o << "  S4: max over s of the tool's distance to the path at q*_used(s): " << fkErr << "\n";
    }
    return 0;
}
