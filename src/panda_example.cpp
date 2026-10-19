#include "se3pf/panda_example.hpp"
#include "se3pf/se3.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

using namespace Eigen;

namespace se3pf {

    namespace {
        // W = diag(l^2 I, I), checked first: every later member depends on it.
        Metric exampleMetric(double l)
        {
            if (!(l > 0.0))
                throw std::invalid_argument("PandaExample: need l > 0");
            return Metric::product(l * l, 1.0);
        }

        // R0 = the flange attitude at the rest posture (approach axis pointing down)
        Matrix3d restAttitude()
        {
            const Pose g = PandaArm::output(PandaArm::restPosture(), nullptr);
            return g.block<3,3>(0,0);
        }
    }  // namespace

    // Members are initialized in the order of their declaration, which is the order of the list below.
    PandaExample::PandaExample(double l, const Vec7& liftWeights)
        : metric_(exampleMetric(l)),
          R0_(restAttitude()),
          path_(R0_, center(), radius, cone, metric_),
          tube_(path_, metric_),
          lift_(path_, PandaArm::restPosture(), 2000, liftWeights),
          kernel_(lift_),
          beta_(lift_, kernel_),
          coupling_(),
          paperStart_(start(beta_, xi0(), 0.0)),
          selfMotionStart_(start(beta_, Vec5::Zero(), selfMotionBeta0))
    {
    }

    Vec7 PandaExample::hardwareLiftWeights()
    {
        // Weight 10 on joint 7 (CTRL-1, 2026-10-19): per loop at nu_d = 0.05 joint 1 then travels 1.26 rad at up to
        // 0.22 rad/s and is below 10 mrad/s for 2 % of the loop (the paper's lift: 0.28 rad, 0.047 rad/s, 42 %),
        // joint 7 still moves 0.90 rad; the joint margin stays 0.358 rad, sigma_min(J_b) 0.18 (0.21), the lift is a
        // minimum (check 4.1e-3), peak joint speeds 13 % of their limits (joint 3). 20 would put joint 1 at 15 %
        // and |q*'| at 10.7 for little more gain.
        return (Vec7() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0).finished();
    }

    Vec7 PandaExample::start(const RedundancyOutput& beta, const Vec5& xi, double beta0) const
    {
        Vec7 guess = Vec7::Zero();
        lift_.eval(s0, guess);
        Vec7 q = Vec7::Zero();
        if (!startConfiguration(tube_, beta, s0, xi, beta0, guess, q))
            throw std::runtime_error("PandaExample: no start configuration with this pose and beta");
        return q;
    }

    Gains PandaExample::gains()
    {
        Gains g = Gains::fromLambda(lambda);  // K_p = 16, K_d = 8 (and K_p,beta, K_d,beta, set again below)
        g.kt = kt;
        g.KpBeta = KpBeta;
        g.KdBeta = KdBeta;
        return g;
    }

    TangentialReference PandaExample::reference()
    {
        return TangentialReference::constantVelocity(nuD);
    }

    double PandaExample::designedResponse(double y0, double t)
    {
        return y0 * (1.0 + lambda * t) * std::exp(-lambda * t);
    }

    double PandaExample::duration(RunKind which) const
    {
        return which == RunKind::Paper ? T : tube_.length() / nuD;
    }

    PandaRun simulate(const PandaExample& example, PandaExample::RunKind which, std::optional<double> T, double dt)
    {
        const double duration = T ? *T : example.duration(which);
        // Checks T and dt (finite, T >= 0, dt > 0, T / dt <= maxSteps) before anything is built or allocated.
        const long steps = stepCount(duration, dt, PandaExample::maxSteps);  // the instants k dt <= T

        PandaRun run;
        run.kind = which;
        run.T = duration;
        run.dt = dt;

        Controller<PandaArm> controller(example.tube(), PandaExample::gains(), PandaExample::reference(),
                                        &example.redundancy());

        // The start and the designed responses: |xi| from |xi0| and beta from beta0, both with zero initial rate.
        const bool paper = which == PandaExample::RunKind::Paper;
        Vec7 q = paper ? example.paperStart() : example.selfMotionStart();
        Vec7 qd = Vec7::Zero();
        const double xi0 = paper ? PandaExample::xi0().norm() : 0.0;
        const double beta0 = paper ? 0.0 : PandaExample::selfMotionBeta0;
        if (!paper) {
            // Moving along the path from the first instant: A(q0) qd0 = (nu_d, 0, ..., 0) gives eta_dot = nu_d and
            // xi_dot = beta_dot = 0 (eq:A), so the transverse and redundancy outputs start at rest.
            Controller<PandaArm>::OutputState out;
            if (!controller.outputs(q, PandaExample::s0, out))
                throw std::runtime_error("simulate: no outputs at the self-motion start");
            qd = out.A.partialPivLu().solve(Vec7(PandaExample::nuD * Vec7::Unit(0)));
            if (!qd.allFinite())
                throw std::runtime_error("simulate: A is singular at the self-motion start");
        }
        // A pose Psi(s0, xi) has eta = s0 (Proposition 1), so the first projection starts at its answer.
        controller.reset(PandaExample::s0);

        const Pose gamma0 = example.path().pose(PandaExample::s0);
        const Pose g0 = PandaArm::output(q, nullptr);
        const Matrix3d Rerr = gamma0.block<3,3>(0,0).transpose() * g0.block<3,3>(0,0);
        run.startAttitudeDeg = logSO3(Rerr).norm() * 180.0 / EIGEN_PI;
        run.startPosition = (g0.block<3,1>(0,3) - gamma0.block<3,1>(0,3)).norm();

        // The joint limits, once: the loop below must not allocate (fixed-size Eigen values do not).
        const Vec7 dqMax = PandaArm::dqMax();
        const Vec7 ddqMax = PandaArm::ddqMax();

        run.log.reserve(static_cast<std::size_t>(steps) + 1);  // every row of the run, even if it stops early
        run.minJointMargin = run.minSigma = std::numeric_limits<double>::infinity();
        run.ok = true;
        for (long k = 0; k <= steps; ++k) {
            const double t = k * dt;
            Vec7 qdd = Vec7::Zero();
            Controller<PandaArm>::Status st;
            if (!controller.step(q, qd, t, dt, qdd, st)) {
                run.ok = false;
                run.failureTime = t;
                break;
            }

            PandaLogRow row;
            row.t = t;
            row.eta = st.eta;
            row.etaDot = st.etaDot;
            row.ref = st.ref;
            row.eTan = st.eTan;
            row.xi = st.xi;
            row.xiNorm = st.xi.norm();
            row.xiDesign = PandaExample::designedResponse(xi0, t);
            row.beta = st.beta;
            row.betaDesign = PandaExample::designedResponse(beta0, t);
            row.q = q;
            row.qd = qd;
            row.qdd = qdd;
            row.sigmaMin = st.sigmaMin;
            int marginJoint = -1;
            row.jointMargin = PandaArm::jointMargin(q, &marginJoint);
            row.iterations = st.iterations;
            run.log.push_back(row);  // within the reserved capacity: no allocation

            // Summary numbers; a strict comparison keeps the first time each extreme is reached.
            const double xiDeviation = std::abs(row.xiNorm - row.xiDesign);
            if (xiDeviation > run.maxXiDeviation) {
                run.maxXiDeviation = xiDeviation;
                run.maxXiDeviationTime = t;
            }
            run.maxXiNorm = std::max(run.maxXiNorm, row.xiNorm);
            if (std::abs(row.beta) > run.maxBeta) {
                run.maxBeta = std::abs(row.beta);
                run.maxBetaTime = t;
            }
            const double betaDeviation = std::abs(row.beta - row.betaDesign);
            if (betaDeviation > run.maxBetaDeviation) {
                run.maxBetaDeviation = betaDeviation;
                run.maxBetaDeviationTime = t;
            }
            if (row.jointMargin < run.minJointMargin) {
                run.minJointMargin = row.jointMargin;
                run.minMarginJoint = marginJoint;
                run.minMarginTime = t;
            }
            if (row.sigmaMin < run.minSigma) {
                run.minSigma = row.sigmaMin;
                run.minSigmaTime = t;
            }
            int joint = -1;
            const double velocityFraction = (qd.cwiseAbs().array() / dqMax.array()).maxCoeff(&joint);
            if (velocityFraction > run.maxVelocityFraction) {
                run.maxVelocityFraction = velocityFraction;
                run.maxVelocityJoint = joint;
                run.maxVelocityTime = t;
            }
            const double accelerationFraction = (qdd.cwiseAbs().array() / ddqMax.array()).maxCoeff(&joint);
            if (accelerationFraction > run.maxAccelerationFraction) {
                run.maxAccelerationFraction = accelerationFraction;
                run.maxAccelerationJoint = joint;
                run.maxAccelerationTime = t;
            }
            run.maxIterations = std::max(run.maxIterations, st.iterations);

            // The arm moves with qdd* held over [t, t + dt): exact for a constant acceleration.
            q += dt * qd + (0.5 * dt * dt) * qdd;
            qd += dt * qdd;
        }

        // A run that stopped early keeps the rows before the failed step, and the numbers below are over those rows.
        // The log is empty only if the very first step fails. That cannot happen with this example's two starts:
        // the first step does not depend on T or dt (the first projection starts from reset(s0)), and it succeeds
        // at both (tested). The guard keeps front() and back() off an empty log all the same.
        if (run.log.empty()) {  // no rows: no extremes to report
            run.minJointMargin = run.minSigma = std::numeric_limits<double>::quiet_NaN();
        } else {
            const PandaLogRow& first = run.log.front();
            const PandaLogRow& last = run.log.back();
            run.startXiNorm = first.xiNorm;
            run.startBeta = first.beta;
            run.finalXiNorm = last.xiNorm;
            run.finalEtaDot = last.etaDot;
            run.finalBeta = last.beta;
            Vec7 qStar = Vec7::Zero();
            example.lift().eval(last.eta, qStar);
            const Vec7 gap = last.q - qStar;
            run.finalLiftGap = gap.cwiseAbs().maxCoeff();
            run.finalLiftGapNorm = gap.norm();
        }
        return run;
    }

}  // namespace se3pf
