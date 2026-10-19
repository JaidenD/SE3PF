#include "se3pf/experiments.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <utility>

namespace se3pf {

    // ---------------------------------------------------------------------------------------------------------
    // SpeedProfile

    void SpeedProfile::eval(double t, double& nu, double& dnuDt) const
    {
        // r(x) = quinticRamp(x), r'(x) = 30 x^2 (1 - x)^2
        nu = 0.0;
        dnuDt = 0.0;
        if (!(t >= 0.0))  // also NaN: no motion
            return;
        if (t < rampUp) {
            const double x = t / rampUp;
            nu = nuMax * quinticRamp(x);
            dnuDt = nuMax * 30.0 * x * x * (1.0 - x) * (1.0 - x) / rampUp;
            return;
        }
        if (t < rampUp + cruise) {
            nu = nuMax;
            return;
        }
        if (t < end()) {
            const double x = (t - rampUp - cruise) / rampDown;
            nu = nuMax * (1.0 - quinticRamp(x));
            dnuDt = -nuMax * 30.0 * x * x * (1.0 - x) * (1.0 - x) / rampDown;
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    // FrictionFeedforward

    bool FrictionFeedforward::active() const
    {
        return (coulomb.array() != 0.0).any() || (viscous.array() != 0.0).any();
    }

    Vec7 FrictionFeedforward::torque(const Vec7& v) const noexcept
    {
        Vec7 tau = Vec7::Zero();
        for (int i = 0; i < 7; ++i)
            tau[i] = coulomb[i] * std::tanh(v[i] / velocityScale) + viscous[i] * v[i];
        return tau;
    }

    FrictionFeedforward FrictionFeedforward::scaled(double k) const
    {
        FrictionFeedforward f = *this;
        f.coulomb = k * coulomb;
        f.viscous = k * viscous;
        return f;
    }

    FrictionFeedforward FrictionFeedforward::gaz2019()
    {
        FrictionFeedforward f;
        // phi1 of the sigmoid, joints 1 to 7 (FI_11 ... FI_17 of the paper's code); Fc = phi1 / 2, the mean of
        // the sigmoid's two asymptotes (experiments.hpp)
        f.coulomb = 0.5 * (Vec7() << 0.54615, 0.87224, 0.64068, 1.2794, 0.83904, 0.30301, 0.56489).finished();
        f.viscous = Vec7::Zero();
        f.velocityScale = 0.01;
        return f;
    }

    // ---------------------------------------------------------------------------------------------------------
    // The catalog

    namespace {
        // Transverse and redundancy poles at -4 (the paper's lambda = 4: K_p = 16, K_d = 8), k_t = 5 for a velocity
        // reference, and tangential position poles at -4 too (k_p = 16, k_d = 8) for a position reference.
        Gains experimentGains()
        {
            Gains g = Gains::fromLambda(PandaExample::lambda);
            g.kt = PandaExample::kt;
            g.kpEta = 16.0;
            g.kdEta = 8.0;
            g.KpBeta = PandaExample::KpBeta;
            g.KdBeta = PandaExample::KdBeta;
            return g;
        }

        std::vector<ExperimentSpec> catalog()
        {
            std::vector<ExperimentSpec> v;
            ExperimentSpec base;
            base.gains = experimentGains();
            base.rampTime = 1.0;

            {
                ExperimentSpec e = base;
                e.name = "hold";
                e.title = "Hold one point of the path";
                e.purpose = "Starts on the path at s0 = 0.1 and regulates eta = s0, xi = 0 and beta = 0 for 10 s. "
                            "Checks the computed-torque loop and the steady-state errors before anything moves "
                            "along the path.";
                e.watch = "The arm should not move visibly. Pass: no safety stop, |xi| and |beta| below 2e-2 "
                          "throughout (an error of the robot's gravity model or friction leaves a small steady "
                          "error; MuJoCo with an exact model: 1e-14).";
                e.mode = TangentialReference::Mode::Position;
                e.duration = 10.0;
                // Nothing should move: a tighter velocity and flange-speed check stops the arm early if it does.
                e.limits.dqFrac = 0.2;
                e.limits.ddqFrac = 0.2;
                e.limits.flangeSpeed = 0.1;
                v.push_back(e);
            }
            {
                ExperimentSpec e = base;
                e.name = "loops";
                e.title = "Two loops at constant speed";
                e.purpose = "Starts on the path and moves along it: nu_d rises smoothly from 0 to 0.05 over 4 s, "
                            "cruises for 26.8 s and falls back to 0 over 4 s (two loops of L = 0.770 if eta_dot "
                            "followed nu_d exactly), then 3.2 s at rest. With beta = 0 the arm follows the closed "
                            "lift q*: the same configuration at the same point of the path in every loop.";
                e.watch = "Pass: no safety stop, |xi| below 2e-2, |beta| below 3e-2 rad, eta_dot near nu_d while "
                          "cruising. Friction: on 2026-10-19 the arm stalled here without compensation (eta_dot 1 % "
                          "of nu_d; |xi| 0.018 and beta 0.030 held by static friction: inside the acceleration "
                          "margin the wrist joints can deliver only 0.02 to 0.06 N m), so the friction feedforward "
                          "tau_f (the Coulomb levels of Gaz et al. 2019, joint 6 raised to 0.25 N m; half of them "
                          "was not enough) is on. If the arm still stalls, raise the joints that did not move. "
                          "An Acceleration stop on a wrist joint soon after the motion "
                          "starts means tau_f exceeds that joint's friction: lower it. If |beta| exceeds 3e-2 here, "
                          "expect 'paper' to stop at its 0.05 rad margin.";
                e.speed = {0.05, 4.0, 26.8, 4.0};
                e.duration = 38.0;
                e.limits.flangeSpeed = 0.2;
                // Full Gaz et al. 2019 levels, joint 6 raised: at half of them (run loops_2, 2026-10-19) the wrist
                // never moved, and the torques held without motion were 0.235 (joint 5), 0.148 (6) and 0.141 N m (7),
                // joint 6 at Gaz's full level 0.152 already.
                e.friction = FrictionFeedforward::gaz2019();
                e.friction.coulomb[5] = 0.25;
                v.push_back(e);
            }
            {
                ExperimentSpec e = base;
                e.name = "recover";
                e.title = "Converge to the path while moving (half the paper's start offset)";
                e.purpose = "Starts off the path at Psi(0.1, xi0 / 2) (|xi| = 0.0156: 7.5 deg and 8.4 mm) with "
                            "beta = 0, and moves along the path: nu_d rises from 0 to 0.08 over 2 s, cruises for 14 s "
                            "and falls back to 0 over the last 4 s of 20 s. xi decays as the paper's designed "
                            "response |xi0| (1 + 4t) e^{-4t}, after the 1 s ramp.";
                e.watch = "Pass: no safety stop; |xi| decays from 0.0156 without overshoot within about 2 s and "
                          "stays below 2e-2; the arm ends at rest. MuJoCo (1 N m s/rad of damping): |beta| up to "
                          "0.015 rad, eta_dot 0.059 instead of 0.08 while cruising, steady |xi| 1.1e-2.";
                e.xi0 = 0.5 * PandaExample::xi0();
                e.speed = {PandaExample::nuD, 2.0, 14.0, 4.0};
                e.duration = 20.0;
                e.limits.flangeSpeed = 0.2;
                v.push_back(e);
            }
            {
                ExperimentSpec e = base;
                e.name = "paper";
                e.title = "The paper's experiment";
                e.purpose = "The paper's Panda run: starts off the path at Psi(0.1, xi0) (|xi0| = 0.0312: 15.1 deg "
                            "and 16.8 mm) with beta = 0, poles at -4, k_t = 5, nu_d = 0.08. On the robot nu_d rises "
                            "over 2 s and falls back to 0 over the last 4 s of 20 s, so the arm starts and ends at "
                            "rest, and mu is ramped in over 1 s.";
                e.watch = "Pass: no safety stop; |xi| follows the designed response for about 2 s and then stays "
                          "below 2e-2; compare with the prediction (plot_robot_log.py --overlay). MuJoCo (1 N m s/rad "
                          "of damping): commanded |qdd*| 17 % of a limit, torque 25 %, |beta| up to 0.015 rad, "
                          "eta_dot 0.059 instead of 0.08 while cruising, steady |xi| 1.1e-2.";
                e.xi0 = PandaExample::xi0();
                e.speed = {PandaExample::nuD, 2.0, 14.0, 4.0};
                e.duration = 20.0;
                e.limits.flangeSpeed = 0.2;
                v.push_back(e);
            }
            {
                ExperimentSpec e = base;
                e.name = "self-motion";
                e.title = "Self-motion back to the lift";
                e.purpose = "Starts on the path at s0 = 0.1 with the arm moved off the lift along ker J_b by "
                            "beta = 0.04 rad, and holds eta = s0: beta decays to zero (the redundancy channel, "
                            "eq:beta-law) while xi stays zero, so the elbow moves and the flange stays.";
                e.watch = "Pass: no safety stop; the flange does not move visibly (|xi| below 2e-2), |beta| below "
                          "1e-2 rad after 3 s.";
                e.mode = TangentialReference::Mode::Position;
                e.beta0 = 0.04;  // inside the 0.05 rad margin of HardwareBounds::beta
                e.duration = 6.0;
                e.limits.dqFrac = 0.2;
                e.limits.flangeSpeed = 0.1;
                v.push_back(e);
            }
            return v;
        }
    }  // namespace

    const std::vector<ExperimentSpec>& experiments()
    {
        static const std::vector<ExperimentSpec> v = catalog();
        return v;
    }

    const ExperimentSpec* findExperiment(const std::string& name)
    {
        for (const ExperimentSpec& e : experiments())
            if (e.name == name)
                return &e;
        return nullptr;
    }

    Vec7 RunConstants::stopDamping()
    {
        return (Vec7() << 30.0, 30.0, 30.0, 30.0, 3.0, 3.0, 1.5).finished();
    }

    Vec7 RunConstants::stopStiffness()
    {
        return (Vec7() << 300.0, 300.0, 300.0, 300.0, 50.0, 50.0, 20.0).finished();
    }

    // ---------------------------------------------------------------------------------------------------------
    // SafetyStop

    void SafetyStop::begin(double t, const Vec7& q) noexcept
    {
        active_ = true;
        t0_ = t;
        restSince_ = -1.0;
        settled_ = false;
        qStop_ = q;
    }

    Vec7 SafetyStop::torque(const Vec7& q, const Vec7& dq, const Vec7& qStop) noexcept
    {
        const Vec7 K = RunConstants::stopStiffness();
        const Vec7 D = RunConstants::stopDamping();
        const Vec7 limit = RunConstants::stopTorqueFrac * PandaArm::tauMax();
        Vec7 tau = Vec7::Zero();
        for (int i = 0; i < 7; ++i) {
            double u = 0.0;
            if (std::isfinite(q[i]) && std::isfinite(qStop[i]))
                u -= K[i] * (q[i] - qStop[i]);
            if (std::isfinite(dq[i]))
                u -= D[i] * dq[i];
            tau[i] = std::clamp(u, -limit[i], limit[i]);  // u is finite here
        }
        return tau;
    }

    Vec7 SafetyStop::step(double t, const Vec7& q, const Vec7& dq, bool& done) noexcept
    {
        if (!active_) {
            done = true;
            return Vec7::Zero();
        }
        // A joint whose q was not finite when the stop began is held where it is first measured again
        for (int i = 0; i < 7; ++i)
            if (!std::isfinite(qStop_[i]) && std::isfinite(q[i]))
                qStop_[i] = q[i];
        const bool atRest = dq.allFinite() && dq.norm() < RunConstants::stopSpeed;
        if (!atRest)
            restSince_ = -1.0;
        else if (restSince_ < 0.0)
            restSince_ = t;
        settled_ = restSince_ >= 0.0 && t - restSince_ >= RunConstants::stopSettle;
        // A NaN time ends the stop too (the robot's own controller takes the arm)
        done = settled_ || !(t - t0_ < RunConstants::stopTime);
        return torque(q, dq, qStop_);
    }

    // ---------------------------------------------------------------------------------------------------------
    // ExperimentSetup

    namespace {
        const ExperimentSpec& validated(const ExperimentSpec& s)
        {
            const SpeedProfile& p = s.speed;
            if (!std::isfinite(s.duration) || !(s.duration > 0.0))
                throw std::invalid_argument("ExperimentSetup: the duration must be finite and positive");
            if (!(s.rampTime >= 0.0) || !std::isfinite(s.rampTime))
                throw std::invalid_argument("ExperimentSetup: the ramp time must be finite and >= 0");
            if (!s.xi0.allFinite() || !std::isfinite(s.beta0))
                throw std::invalid_argument("ExperimentSetup: the start must be finite");
            if (s.mode == TangentialReference::Mode::Velocity) {
                if (!(p.rampUp >= 0.0 && p.cruise >= 0.0 && p.rampDown >= 0.0) || !std::isfinite(p.nuMax)
                    || !std::isfinite(p.end()))
                    throw std::invalid_argument("ExperimentSetup: the speed profile needs finite times >= 0");
                if (!(p.end() <= s.duration))
                    throw std::invalid_argument("ExperimentSetup: the speed profile must end by the duration");
            }
            const FrictionFeedforward& f = s.friction;
            if (!f.coulomb.allFinite() || !f.viscous.allFinite() || (f.coulomb.array() < 0.0).any()
                || (f.viscous.array() < 0.0).any())
                throw std::invalid_argument("ExperimentSetup: the friction feedforward needs finite Fc and Fv >= 0");
            if (!std::isfinite(f.velocityScale) || !(f.velocityScale > 0.0))
                throw std::invalid_argument("ExperimentSetup: the friction feedforward needs a velocity scale > 0");
            return s;
        }
    }  // namespace

    ExperimentSetup::ExperimentSetup(const PandaExample& example, const ExperimentSpec& spec)
        : example_(example), spec_(validated(spec)),
          q0_(example.start(example.redundancy(), spec.xi0, spec.beta0))
    {
    }

    TangentialReference ExperimentSetup::reference() const
    {
        if (spec_.mode == TangentialReference::Mode::Position) {
            const double etaD = s0();
            return TangentialReference::fromPosition([etaD](double, double& e, double& de, double& dde) {
                e = etaD;
                de = 0.0;
                dde = 0.0;
            });
        }
        const SpeedProfile profile = spec_.speed;
        return TangentialReference::fromVelocity(
            [profile](double t, double, double& nu, double& dnuDt, double& dnuDeta) {
                profile.eval(t, nu, dnuDt);
                dnuDeta = 0.0;
            });
    }

    double ExperimentSetup::referenceSpeed(double t) const
    {
        if (spec_.mode == TangentialReference::Mode::Position)
            return 0.0;  // eta_d = s0 (reference()): the point does not move
        double nu = 0.0, dnuDt = 0.0;
        spec_.speed.eval(t, nu, dnuDt);
        return nu;
    }

    Controller<PandaArm> ExperimentSetup::controller() const
    {
        Controller<PandaArm> c(example_.tube(), spec_.gains, reference(), &example_.redundancy());
        c.setRampTime(spec_.rampTime);
        c.reset(s0());  // a pose Psi(s0, xi) has eta = s0 (Proposition 1)
        return c;
    }

    // ---------------------------------------------------------------------------------------------------------
    // ExperimentLoop

    ExperimentLoop::ExperimentLoop(const ExperimentSetup& setup, const SafetyMonitor& monitor)
        : setup_(setup), monitor_(monitor), controller_(setup.controller())
    {
    }

    void ExperimentLoop::follow(double t, const Vec7& q, const Vec7& dq, const ModelTerms* model, Cycle& out)
    {
        out.violation = Violation::None;
        out.finished = false;
        out.tau.setZero();
        out.tauG.setZero();
        out.tauF.setZero();
        out.qdd.setZero();

        if (!q.allFinite() || !dq.allFinite() || !std::isfinite(t)) {
            out.status = Status();
            out.status.t = t;
            out.violation = Violation::NonFinite;
            return;
        }
        // The controller's dt is only the prediction horizon of eta-hat (the next cycle comes after one period).
        Vec7 qdd = Vec7::Zero();
        if (!controller_.step(q, dq, t, RunConstants::dt, qdd, out.status)) {
            out.violation = Violation::Projection;
            return;
        }
        eta_ = out.status.eta;
        out.qdd = qdd;

        Violation v = monitor_.checkOutputs(out.status.xi, out.status.beta, out.status.sigmaMin);
        if (v == Violation::None)
            v = monitor_.checkState(q, dq);
        if (v == Violation::None)
            v = monitor_.checkAcceleration(qdd);
        // The friction feedforward of the reference motion: v = q*'(eta) eta_dot_ref (FrictionFeedforward)
        const FrictionFeedforward& friction = setup_.spec().friction;
        if (friction.active()) {
            Vec7 qStar = Vec7::Zero(), dqStar = Vec7::Zero();
            setup_.example().lift().eval(out.status.eta, qStar, &dqStar);
            out.tauF = friction.torque(Vec7(dqStar * setup_.referenceSpeed(t)));
            if (v == Violation::None && !out.tauF.allFinite())
                v = Violation::NonFinite;
        }
        if (model) {
            // eq:fb without G (the robot adds its own G, with the configured load), plus the friction feedforward
            out.tau = model->M * qdd + model->c + out.tauF;
            out.tauG = out.tau + model->G;
            if (v == Violation::None)
                v = monitor_.checkTorque(out.tauG, out.tau, tauPrev_, RunConstants::dt);
            tauPrev_ = out.tau;
        }
        if (v == Violation::None && t > setup_.timeout())
            v = Violation::Timeout;
        out.violation = v;
        out.finished = v == Violation::None && t >= setup_.duration() && dq.norm() < RunConstants::restSpeed;
    }

    bool ExperimentLoop::observe(const Vec7& q, const Vec7& dq, Status& st) const
    {
        st = Status();
        Controller<PandaArm>::OutputState out;
        if (!q.allFinite() || !dq.allFinite() || !controller_.outputs(q, eta_, out))
            return false;
        const Vec7 yDot = out.A * dq;
        st.eta = out.eta;
        st.etaDot = yDot[0];
        st.xi = out.xi;
        st.xiDot = yDot.segment<5>(1);
        st.beta = out.beta;
        st.betaDot = yDot[6];
        st.iterations = out.iterations;
        st.ok = true;
        return true;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Logs

    ExperimentLogRow makeLogRow(double t, int phase, const ExperimentLoop::Status& st, const Vec7& q, const Vec7& dq,
                                const Vec7& qdd, const Vec7& tau, const Vec7& tauG, const Vec7& tauF)
    {
        ExperimentLogRow r;
        r.t = t;
        r.phase = phase;
        r.eta = st.eta;
        r.etaDot = st.etaDot;
        r.ref = st.ref;
        r.eTan = st.eTan;
        r.xi = st.xi;
        r.xiNorm = st.xi.norm();
        r.beta = st.beta;
        r.betaDot = st.betaDot;
        r.sigmaMin = st.sigmaMin;
        r.ramp = st.ramp;
        r.newton = st.iterations;
        r.q = q;
        r.dq = dq;
        r.qdd = qdd;
        r.tau = tau;
        r.tauG = tauG;
        r.tauF = tauF;
        return r;
    }

    void RunPeaks::add(const ExperimentLogRow& r)
    {
        if (r.phase != 0)
            return;
        if (rows == 0) {
            minJointMargin = minSigma = std::numeric_limits<double>::infinity();
        }
        ++rows;
        int j = -1;
        const double margin = PandaArm::jointMargin(r.q, &j);
        if (margin < minJointMargin || std::isnan(margin)) {
            minJointMargin = margin;
            minMarginJoint = j;
        }
        const double dqFrac = r.dq.cwiseAbs().cwiseQuotient(PandaArm::dqMax()).maxCoeff(&j);
        if (dqFrac > maxDqFrac) {
            maxDqFrac = dqFrac;
            maxDqJoint = j;
        }
        const double ddqFrac = r.qdd.cwiseAbs().cwiseQuotient(PandaArm::ddqMax()).maxCoeff(&j);
        if (ddqFrac > maxDdqFrac) {
            maxDdqFrac = ddqFrac;
            maxDdqJoint = j;
        }
        if (torques) {
            const double tauFrac = r.tauG.cwiseAbs().cwiseQuotient(PandaArm::tauMax()).maxCoeff(&j);
            if (tauFrac > maxTauFrac) {
                maxTauFrac = tauFrac;
                maxTauJoint = j;
            }
            maxTauRate = std::max(maxTauRate, (r.tau - tauPrev_).cwiseAbs().maxCoeff() / RunConstants::dt);
            tauPrev_ = r.tau;
        }
        const double tauF = r.tauF.cwiseAbs().maxCoeff(&j);
        if (tauF > maxTauF) {
            maxTauF = tauF;
            maxTauFJoint = j;
        }
        Mat67 Jb = Mat67::Zero();
        PandaArm::output(r.q, &Jb);
        const Vec6 zeta = Jb * r.dq;
        maxFlangeSpeed = std::max(maxFlangeSpeed, zeta.tail<3>().norm());
        maxXi = std::max(maxXi, r.xiNorm);
        maxBeta = std::max(maxBeta, std::abs(r.beta));
        minSigma = std::min(minSigma, r.sigmaMin);
        maxNewton = std::max(maxNewton, r.newton);
        finalT = r.t;
        finalXi = r.xiNorm;
        finalBeta = r.beta;
        finalEtaDot = r.etaDot;
        finalETan = r.eTan;
    }

    std::string RunPeaks::summary(const SafetyLimits& b) const
    {
        char buf[900];
        char torque[200];
        if (torques)
            std::snprintf(torque, sizeof(torque), " | torque %.1f %% (joint %d) [%.0f %%] | torque rate %.1f N m/s [%.0f]",
                          100.0 * maxTauFrac, maxTauJoint + 1, 100.0 * b.tauFrac, maxTauRate,
                          b.tauRateFrac * PandaArm::tauRateMax);
        else
            std::snprintf(torque, sizeof(torque), " | torque not checked (no dynamics model)");
        std::snprintf(buf, sizeof(buf),
                      "joint margin %.3f rad (joint %d) [>= %.2f] | |qd| %.1f %% (joint %d) [%.0f %%] | |qdd*| %.1f %% "
                      "(joint %d) [%.0f %%]%s | tau_f <= %.2f N m (joint %d) | flange %.3f m/s [%.2f] | |xi| %.4f "
                      "[%.3f] | |beta| %.4f rad [%.3f] | sigma_min(A) %.4f [>= %.0e] | Newton <= %d | final |xi| "
                      "%.1e, beta %.1e",
                      minJointMargin, minMarginJoint + 1, b.jointMargin, 100.0 * maxDqFrac, maxDqJoint + 1,
                      100.0 * b.dqFrac, 100.0 * maxDdqFrac, maxDdqJoint + 1, 100.0 * b.ddqFrac, torque, maxTauF,
                      maxTauFJoint + 1, maxFlangeSpeed, b.flangeSpeed, maxXi, b.xi, maxBeta, b.beta, minSigma,
                      b.sigmaMin, maxNewton, finalXi, finalBeta);
        return buf;
    }

}  // namespace se3pf
