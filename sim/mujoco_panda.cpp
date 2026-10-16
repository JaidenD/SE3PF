#include "mujoco_panda.hpp"

#include <mujoco/mujoco.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

#ifndef SE3PF_MENAGERIE_PANDA_DIR
#error "SE3PF_MENAGERIE_PANDA_DIR must be set by cmake/mujoco.cmake (an absolute path to the Menagerie Panda)"
#endif

namespace se3pf {

    namespace {
        constexpr int kArm = MujocoPanda::kArm;

        // Turn the Menagerie position servos (force = kp (ctrl - q) - kd qdot) into torque motors, force = ctrl, for
        // the seven arm actuators. ctrl is clamped to the actuator's force range, which Menagerie sets to the Franka
        // torque limits. Returns false unless actuator i drives joint i (a hinge) through control i, force output i
        // and dof i, with unit gear, a force limit and no actuator dynamics, damping, armature or delay, so that
        // ctrl[i] = tau_i exactly. apps/mujoco_smoke.cpp has a copy with the same checks (it avoids se3pf on purpose):
        // keep the two in sync.
        bool useTorqueMotors(mjModel* m)
        {
            if (m->nactuator < kArm || m->njnt < kArm || m->nv < kArm || m->nu < kArm)
                return false;
            for (int i = 0; i < kArm; ++i) {
                const int c = m->actuator_ctrladr[i];
                if (m->jnt_type[i] != mjJNT_HINGE || m->jnt_qposadr[i] != i || m->jnt_dofadr[i] != i
                    || m->actuator_trntype[i] != mjTRN_JOINT || m->actuator_trnid[2 * i] != i || c != i
                    || m->actuator_ctrlnum[i] != 1 || m->actuator_outadr[i] != i || m->actuator_gear[6 * i] != 1.0
                    || m->actuator_dyntype[i] != mjDYN_NONE || m->actuator_damping[i] != 0.0
                    || m->actuator_armature[i] != 0.0 || m->actuator_delay[i] != 0.0
                    || !m->actuator_forcelimited[i])
                    return false;
                m->actuator_gaintype[i] = mjGAIN_FIXED;
                m->actuator_biastype[i] = mjBIAS_NONE;
                std::fill(m->actuator_gainprm + i * mjNGAIN, m->actuator_gainprm + (i + 1) * mjNGAIN, 0.0);
                std::fill(m->actuator_biasprm + i * mjNBIAS, m->actuator_biasprm + (i + 1) * mjNBIAS, 0.0);
                m->actuator_gainprm[i * mjNGAIN] = 1.0;
                m->actuator_ctrllimited[c] = 1;
                m->actuator_ctrlrange[2 * c] = m->actuator_forcerange[2 * i];
                m->actuator_ctrlrange[2 * c + 1] = m->actuator_forcerange[2 * i + 1];
            }
            return true;
        }

        // The keyframe "home" of the Menagerie models (fingers open at 0.04 m with ctrl = 255 in panda.xml). Its arm
        // controls are position targets of the servos, which would be torques now, so they are zeroed.
        void resetHome(const mjModel* m, mjData* d)
        {
            const int key = mj_name2id(m, mjOBJ_KEY, "home");
            if (key >= 0)
                mj_resetDataKeyframe(m, d, key);
            else
                mj_resetData(m, d);
            for (int i = 0; i < kArm; ++i)
                d->ctrl[i] = 0.0;
        }

        Vec7 armVector(const mjtNum* v)
        {
            Vec7 x = Vec7::Zero();
            for (int i = 0; i < kArm; ++i)
                x[i] = v[i];
            return x;
        }

        void setArm(mjtNum* dst, const Vec7& x)
        {
            for (int i = 0; i < kArm; ++i)
                dst[i] = x[i];
        }

        double maxAbs(const Vec7& v) { return v.cwiseAbs().maxCoeff(); }
    }  // namespace

    // ---------------------------------------------------------------------------------------------------------
    // MujocoPanda

    void MujocoPanda::ModelDeleter::operator()(mjModel* m) const { mj_deleteModel(m); }
    void MujocoPanda::DataDeleter::operator()(mjData* d) const { mj_deleteData(d); }

    MujocoPanda::MujocoPanda(const std::string& xmlPath, const Options& options) : options_(options)
    {
        if (!(options.timestep > 0.0))
            throw std::invalid_argument("MujocoPanda: the timestep must be positive");

        char error[1000] = "";
        m_.reset(mj_loadXML(xmlPath.c_str(), nullptr, error, sizeof(error)));
        if (!m_)
            throw std::runtime_error("MujocoPanda: cannot load " + xmlPath + ": " + error);
        mjModel* m = m_.get();
        link7_ = mj_name2id(m, mjOBJ_BODY, "link7");
        if (link7_ < 0 || !useTorqueMotors(m))
            throw std::runtime_error("MujocoPanda: " + xmlPath + " is not a Menagerie Panda (body link7, actuator i "
                                     "on joint i with unit gear)");

        m->opt.timestep = options.timestep;
        m->opt.integrator = options.integrator == Integrator::RK4 ? mjINT_RK4 : mjINT_EULER;
        if (!options.damping)
            m->opt.disableflags |= mjDSBL_DAMPER;
        // MuJoCo's Euler integrates joint damping implicitly by default: it integrates the acceleration
        // (M + h D)^-1 (tau - bias + passive), not the qacc = M^-1 (tau - bias + passive) it reports, an O(h)
        // difference (0.035 rad/s^2 at |qdot| ~ 0.3 rad/s). Explicit damping makes the reported qacc the integrated
        // one, so the exact step stays exact with the damping on. (RK4 is always explicit.)
        m->opt.disableflags |= mjDSBL_EULERDAMP;
        m->opt.disableflags |= mjDSBL_AUTORESET;  // a bad state stays visible (warning, NaN) instead of a reset

        d_.reset(mj_makeData(m));
        scratch_.reset(mj_makeData(m));
        if (!d_ || !scratch_)
            throw std::runtime_error("MujocoPanda: cannot allocate mjData");
        resetHome(m, d_.get());
        resetHome(m, scratch_.get());

        const std::size_t nv = static_cast<std::size_t>(m->nv);
        fullM_.assign(nv * nv, 0.0);
        jacp_.assign(3 * nv, 0.0);
        jacr_.assign(3 * nv, 0.0);
        vecIn_.assign(nv, 0.0);
        vecOut_.assign(nv, 0.0);
        mj_forward(m, d_.get());
    }

    bool MujocoPanda::hasHand() const { return m_->nq > kArm; }

    Vec7 MujocoPanda::torqueLimit() const
    {
        Vec7 limit = Vec7::Zero();
        for (int i = 0; i < kArm; ++i)
            limit[i] = m_->actuator_ctrlrange[2 * m_->actuator_ctrladr[i] + 1];
        return limit;
    }

    void MujocoPanda::setState(const Vec7& q, const Vec7& qd)
    {
        setArm(d_->qpos, q);
        setArm(d_->qvel, qd);
        mj_forward(m_.get(), d_.get());
    }

    Vec7 MujocoPanda::q() const { return armVector(d_->qpos); }
    Vec7 MujocoPanda::qd() const { return armVector(d_->qvel); }
    Vec7 MujocoPanda::qacc() const { return armVector(d_->qacc); }
    double MujocoPanda::time() const { return d_->time; }

    void MujocoPanda::dynamics(Mat<7, 7>& M, Vec7& bias, Vec7& passive)
    {
        mj_forward(m_.get(), d_.get());
        // MuJoCo keeps M sparse (there is no d->qM since 3.15); mj_fullM writes the dense nv x nv matrix, row-major.
        mj_fullM(m_.get(), d_.get(), fullM_.data());
        const int nv = m_->nv;
        for (int r = 0; r < kArm; ++r)
            for (int c = 0; c < kArm; ++c)
                M(r, c) = fullM_[static_cast<std::size_t>(r * nv + c)];
        bias = armVector(d_->qfrc_bias);       // C(q, qdot) qdot + G(q), by recursive Newton-Euler
        passive = armVector(d_->qfrc_passive); // -D qdot (zero with the damper switched off)
    }

    Vec7 MujocoPanda::gravity(const Vec7& q)
    {
        mjData* s = scratch_.get();
        setArm(s->qpos, q);
        for (int i = 0; i < m_->nv; ++i)
            s->qvel[i] = 0.0;
        mj_forward(m_.get(), s);
        return armVector(s->qfrc_bias);  // C qdot = 0 at rest, so the bias is G
    }

    Pose MujocoPanda::flange(const Vec7& q, Mat67* Jb)
    {
        mjData* s = scratch_.get();
        setArm(s->qpos, q);
        mj_kinematics(m_.get(), s);
        mj_comPos(m_.get(), s);  // mj_jac needs the frames of mj_comPos

        // xmat is row-major: R(r, c) = xmat[3 r + c]; the flange is 0.107 m along link7's z-axis (column 2).
        const mjtNum* xm = s->xmat + 9 * link7_;
        Matrix3d R;
        R << xm[0], xm[1], xm[2],
             xm[3], xm[4], xm[5],
             xm[6], xm[7], xm[8];
        const Vector3d p = Vector3d(s->xpos[3 * link7_], s->xpos[3 * link7_ + 1], s->xpos[3 * link7_ + 2])
                           + PandaArm::kFlangeOffset * R.col(2);
        Pose g = Pose::Identity();
        g.block<3,3>(0,0) = R;
        g.block<3,1>(0,3) = p;

        if (Jb) {
            // mj_jac: world-frame rotational and translational Jacobians of the point p fixed to link7, 3 x nv each,
            // row-major. In the flange frame: omega_b = R^T omega, v_b = R^T pdot.
            const mjtNum point[3] = {p[0], p[1], p[2]};
            mj_jac(m_.get(), s, jacp_.data(), jacr_.data(), point, link7_);
            const int nv = m_->nv;
            Eigen::Matrix<double, 3, 7> jr = Eigen::Matrix<double, 3, 7>::Zero();
            Eigen::Matrix<double, 3, 7> jp = Eigen::Matrix<double, 3, 7>::Zero();
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < kArm; ++c) {
                    jr(r, c) = jacr_[static_cast<std::size_t>(r * nv + c)];
                    jp(r, c) = jacp_[static_cast<std::size_t>(r * nv + c)];
                }
            Jb->topRows<3>() = R.transpose() * jr;
            Jb->bottomRows<3>() = R.transpose() * jp;
        }
        return g;
    }

    Vec7 MujocoPanda::exactStep(const Vec7& qdd)
    {
        mjModel* m = m_.get();
        mjData* d = d_.get();
        mj_step1(m, d);  // positions and velocities: M, qfrc_bias = C qdot + G, qfrc_passive

        // tau = M qdd + C qdot + G - passive (eq:fb, with MuJoCo's own model). The finger dofs of panda.xml (if
        // any) get zero acceleration in the product; the exact step is meant for panda_nohand.xml.
        std::fill(vecIn_.begin(), vecIn_.end(), 0.0);
        setArm(vecIn_.data(), qdd);
        mj_mulM(m, d, vecOut_.data(), vecIn_.data());
        Vec7 tau = Vec7::Zero();
        for (int i = 0; i < kArm; ++i) {
            tau[i] = vecOut_[static_cast<std::size_t>(i)] + d->qfrc_bias[i] - d->qfrc_passive[i];
            d->ctrl[i] = tau[i];
        }
        mj_step2(m, d);  // actuation, qacc, and one Euler step (mj_step1 / mj_step2 always integrate with Euler)
        return tau;
    }

    void MujocoPanda::holdStep(const Vec7& tau, int substeps)
    {
        setArm(d_->ctrl, tau);
        for (int i = 0; i < substeps; ++i)
            mj_step(m_.get(), d_.get());
    }

    void MujocoPanda::setAppliedTorque(const Vec7& tau) { setArm(d_->qfrc_applied, tau); }

    int MujocoPanda::contacts() const { return d_->ncon; }
    int MujocoPanda::constraintRows() const { return d_->nefc; }

    int MujocoPanda::warnings() const
    {
        int n = 0;
        for (int i = 0; i < mjNWARNING; ++i)
            n += d_->warning[i].number;
        return n;
    }

    bool MujocoPanda::withinTorqueLimits(const Vec7& tau) const
    {
        for (int i = 0; i < kArm; ++i) {
            const int c = m_->actuator_ctrladr[i];
            if (!(tau[i] >= m_->actuator_ctrlrange[2 * c] && tau[i] <= m_->actuator_ctrlrange[2 * c + 1]))
                return false;  // also NaN
        }
        return true;
    }

    std::string defaultPandaModelDir() { return SE3PF_MENAGERIE_PANDA_DIR; }

    // ---------------------------------------------------------------------------------------------------------
    // The hardware margins

    HardwareMargins hardwareMargins(const Vec7& q, const Vec7& qd, const Vec7& qdd, const Vec7& tau, const Vec7& G,
                                    const Vec7& tauPrev, const Vec7& GPrev, double dt)
    {
        HardwareMargins h;
        h.jointMargin = PandaArm::jointMargin(q, &h.marginJoint);
        // Each joint against its own limit (cwiseQuotient), never against one limit for all
        h.dqFrac = qd.cwiseAbs().cwiseQuotient(PandaArm::dqMax());
        h.ddqFrac = qdd.cwiseAbs().cwiseQuotient(PandaArm::ddqMax());
        h.tauFrac = tau.cwiseAbs().cwiseQuotient(PandaArm::tauMax());
        // libfranka is sent tau - G (it adds G itself), so the rate margin applies to that difference
        h.tauRateCmd = maxAbs((tau - G) - (tauPrev - GPrev)) / dt;
        h.tauRate = maxAbs(tau - tauPrev) / dt;
        Mat67 Jb = Mat67::Zero();
        PandaArm::output(q, &Jb);
        const Vec6 zeta = Jb * qd;  // body twist (omega, v) of the flange
        h.flangeSpeed = zeta.tail<3>().norm();
        return h;
    }

    // ---------------------------------------------------------------------------------------------------------
    // The closed loop

    Controller<PandaArm> exampleController(const PandaExample& example)
    {
        // The example's own gains and reference, so this runs the controller that the joint-acceleration-level
        // simulate() and its tests check.
        Controller<PandaArm> controller(example.tube(), PandaExample::gains(), PandaExample::reference(),
                                        &example.redundancy());
        controller.reset(PandaExample::s0);
        return controller;
    }

    namespace {
        // max_i |q_i - q*_i(eta)| and its 2-norm at q, with eta projected from the guess (the controller's
        // prediction). NaN if the projection fails.
        void liftError(const PandaExample& example, const Controller<PandaArm>& controller, const Vec7& q,
                       double& maxErr, double& norm2)
        {
            Controller<PandaArm>::OutputState out;
            maxErr = norm2 = std::nan("");
            if (!controller.outputs(q, controller.prediction(), out))
                return;
            Vec7 qStar = Vec7::Zero();
            example.lift().eval(out.eta, qStar);
            maxErr = (q - qStar).cwiseAbs().maxCoeff();
            norm2 = (q - qStar).norm();
        }
    }  // namespace

    MujocoRun simulateMujoco(const PandaExample& example, const MujocoRunOptions& options)
    {
        using Start = MujocoRunOptions::Start;
        using Torque = MujocoRunOptions::Torque;
        if (options.plant.empty() || !(options.dt > 0.0) || options.substeps < 1 || !(options.duration >= 0.0))
            throw std::invalid_argument("simulateMujoco: need a plant model, dt > 0, substeps >= 1, duration >= 0");
        // The control periods t_k = k dt <= T: T = duration for the paper's start, one loop L / nu_d for the
        // self-motion start (9627 periods at 1 ms). Also checks that T and dt are finite and T / dt is in range.
        const double T = options.start == MujocoRunOptions::Start::Paper
                             ? options.duration
                             : example.duration(PandaExample::RunKind::SelfMotion);
        const long periods = stepCount(T, options.dt, PandaExample::maxSteps) + 1;
        if (options.torque == Torque::Exact && !options.model.empty())
            throw std::invalid_argument("simulateMujoco: the exact model is the plant's own model");
        if (options.plantGravity && options.model.empty())
            throw std::invalid_argument("simulateMujoco: plantGravity needs a separate model file");

        MujocoRun run;
        run.options = options;

        // The plant, and the model of the computed torque (the plant itself unless another file is given)
        MujocoPanda::Options plantOptions;
        plantOptions.damping = options.damping;
        plantOptions.timestep = options.dt / options.substeps;
        // mj_step1 / mj_step2 of the exact step always integrate with Euler
        plantOptions.integrator = options.torque == Torque::Exact ? MujocoPanda::Integrator::Euler
                                                                  : options.integrator;
        MujocoPanda plant(options.plant, plantOptions);
        if (options.torque == Torque::Exact && plant.hasHand())
            throw std::invalid_argument("simulateMujoco: the exact step needs panda_nohand.xml (the finger dofs "
                                        "couple to the arm through M)");
        std::unique_ptr<MujocoPanda> separateModel;
        if (!options.model.empty())
            separateModel = std::make_unique<MujocoPanda>(options.model, plantOptions);
        MujocoPanda& model = separateModel ? *separateModel : plant;

        Controller<PandaArm> controller = exampleController(example);
        controller.setRampTime(options.rampTime);

        // Start (design doc, Panda simulation): at rest at the paper's start, or on the path with beta = 0.2 and
        // the velocity A(q0)^-1 (nu_d, 0, ..., 0), so that eta_dot = nu_d and xi_dot = beta_dot = 0 at t = 0, as the
        // designed response beta0 (1 + 4t) e^{-4t} assumes.
        Vec7 q0 = Vec7::Zero(), qd0 = Vec7::Zero();
        if (options.start == Start::Paper) {
            q0 = example.paperStart();
            run.xi0 = PandaExample::xi0().norm();
            run.beta0 = 0.0;
        } else {
            q0 = example.selfMotionStart();
            run.xi0 = 0.0;
            run.beta0 = PandaExample::selfMotionBeta0;
            Controller<PandaArm>::OutputState out;
            if (!controller.outputs(q0, PandaExample::s0, out))
                throw std::runtime_error("simulateMujoco: no outputs at the self-motion start");
            Vec7 y0 = Vec7::Zero();
            y0[0] = PandaExample::nuD;
            qd0 = out.A.partialPivLu().solve(y0);
        }
        plant.setState(q0, qd0);

        // The baseline: the same controller at the joint-acceleration level, from the same start
        Controller<PandaArm> kinematic = exampleController(example);
        kinematic.setRampTime(options.rampTime);
        Vec7 qKin = q0, qdKin = qd0;
        bool baselineOk = options.baseline;

        if (options.keepLog)
            run.log.reserve(static_cast<std::size_t>(periods));

        // Before t = 0 the arm rests at q0, held by gravity compensation: the delivered torque is G(q0) and the
        // robot's commanded torque tau - G is 0. The first torque-rate samples compare with these.
        const Vec7 G0 = options.plantGravity ? plant.gravity(q0) : model.gravity(q0);
        Vec7 tauPrev = G0, GPrev = G0;

        const Vec7 dqMax = PandaArm::dqMax(), ddqMax = PandaArm::ddqMax();
        const double dt = options.dt;
        const double steadyFrom = std::max(0.0, (periods - 1) * dt - 5.0);

        Mat<7, 7> M = Mat<7, 7>::Zero();
        Vec7 bias = Vec7::Zero(), passive = Vec7::Zero();
        Mat<7, 7> Mp = Mat<7, 7>::Zero();
        Vec7 biasP = Vec7::Zero(), passiveP = Vec7::Zero();

        const auto wallStart = std::chrono::steady_clock::now();
        run.ok = true;
        for (long k = 0; k < periods; ++k) {
            const double t = k * dt;
            const Vec7 q = plant.q();
            const Vec7 qd = plant.qd();

            Vec7 qdd = Vec7::Zero();
            Controller<PandaArm>::Status st;
            if (!q.allFinite() || !qd.allFinite() || !controller.step(q, qd, t, dt, qdd, st)) {
                run.ok = false;
                run.failureTime = t;
                break;
            }

            // The torque, and the plant moves on to t + dt
            // G is the gravity term inside tau: tau - G is what libfranka would be sent (it adds G itself).
            Vec7 tau = Vec7::Zero(), G = Vec7::Zero();
            int contacts = 0, rows = 0;
            if (options.torque == Torque::Exact) {
                G = plant.gravity(q);
                // qdd* held for the period; the torque follows the plant's state at every plant step
                for (int j = 0; j < options.substeps; ++j) {
                    const Vec7 tauJ = plant.exactStep(qdd);
                    if (j == 0) {
                        tau = tauJ;  // logged: the torque at the start of the period
                        contacts = plant.contacts();  // of mj_step1, at (q, qd)
                        rows = plant.constraintRows();
                    }
                    run.maxQaccError = std::max(run.maxQaccError, maxAbs(plant.qacc() - qdd));
                }
            } else {
                plant.dynamics(Mp, biasP, passiveP);  // evaluates the plant at (q, qd): contacts and rows
                contacts = plant.contacts();
                rows = plant.constraintRows();
                if (separateModel) {
                    model.setState(q, qd);
                    model.dynamics(M, bias, passive);
                    G = model.gravity(q);
                    if (options.plantGravity) {
                        const Vec7 Gplant = plant.gravity(q);
                        bias += Gplant - G;
                        G = Gplant;
                    }
                } else {
                    M = Mp;
                    bias = biasP;
                    G = plant.gravity(q);
                }
                tau = M * qdd + bias;  // eq:fb from the state at the start of the period; passive not compensated
                plant.holdStep(tau, options.substeps);
            }

            // Margins of the hardware section (design doc), at the start of the period
            const HardwareMargins hm = hardwareMargins(q, qd, qdd, tau, G, tauPrev, GPrev, dt);
            tauPrev = tau;
            GPrev = G;
            if (hm.jointMargin < run.minJointMargin) {
                run.minJointMargin = hm.jointMargin;
                run.minMarginJoint = hm.marginJoint;
            }
            run.peakDqFrac = run.peakDqFrac.cwiseMax(hm.dqFrac);
            run.peakDdqFrac = run.peakDdqFrac.cwiseMax(hm.ddqFrac);
            run.peakTauFrac = run.peakTauFrac.cwiseMax(hm.tauFrac);
            run.maxDqFrac = std::max(run.maxDqFrac, hm.dqFrac.maxCoeff());
            const double ddqFrac = hm.ddqFrac.maxCoeff();
            if (ddqFrac > run.maxDdqFrac) {
                run.maxDdqFrac = ddqFrac;
                run.maxDdqFracTime = t;
            }
            run.maxTauFrac = std::max(run.maxTauFrac, hm.tauFrac.maxCoeff());
            run.maxTauRateCmd = std::max(run.maxTauRateCmd, hm.tauRateCmd);
            run.maxTauRate = std::max(run.maxTauRate, hm.tauRate);
            run.maxFlangeSpeed = std::max(run.maxFlangeSpeed, hm.flangeSpeed);
            run.minSigma = std::min(run.minSigma, st.sigmaMin);
            run.maxIterations = std::max(run.maxIterations, st.iterations);
            run.maxContacts = std::max(run.maxContacts, contacts);
            run.maxConstraintRows = std::max(run.maxConstraintRows, rows);
            if (!plant.withinTorqueLimits(tau))
                ++run.clampedSteps;

            // Transverse and redundancy outputs against their designed responses
            const double xiNorm = st.xi.norm();
            const double xiDesign = PandaExample::designedResponse(run.xi0, t);
            const double betaDesign = PandaExample::designedResponse(run.beta0, t);
            const double deviation = std::abs(xiNorm - xiDesign);
            if (deviation > run.maxDeviation) {
                run.maxDeviation = deviation;
                run.maxDeviationTime = t;
            }
            run.maxXi = std::max(run.maxXi, xiNorm);
            run.maxBeta = std::max(run.maxBeta, std::abs(st.beta));
            run.maxBetaDeviation = std::max(run.maxBetaDeviation, std::abs(st.beta - betaDesign));
            if (t >= steadyFrom) {
                run.steadyXi = std::max(run.steadyXi, xiNorm);
                run.steadyBeta = std::max(run.steadyBeta, std::abs(st.beta));
            }
            run.finalXi = xiNorm;
            run.finalBeta = st.beta;
            run.finalEtaDot = st.etaDot;

            // Kinematics cross-check: the library's FK and J_b against MuJoCo's
            if (options.crossCheck) {
                Mat67 Jb = Mat67::Zero(), JbMj = Mat67::Zero();
                const Pose g = PandaArm::output(q, &Jb);
                const Pose gMj = plant.flange(q, &JbMj);
                run.maxPoseError = std::max(run.maxPoseError, (gMj - g).cwiseAbs().maxCoeff());
                run.maxJacobianError = std::max(run.maxJacobianError, (JbMj - Jb).cwiseAbs().maxCoeff());
            }

            // The baseline at the same instant, then its exact update for an acceleration held over the period
            double xiKin = 0.0, betaKin = 0.0, qErr = 0.0;
            if (baselineOk) {
                Vec7 qddKin = Vec7::Zero();
                Controller<PandaArm>::Status stKin;
                if (kinematic.step(qKin, qdKin, t, dt, qddKin, stKin)) {
                    xiKin = stKin.xi.norm();
                    betaKin = stKin.beta;
                    qErr = maxAbs(q - qKin);
                    run.maxBaselineQError = std::max(run.maxBaselineQError, qErr);
                    run.baselineMaxDeviation = std::max(run.baselineMaxDeviation, std::abs(xiKin - xiDesign));
                    run.baselineMaxXi = std::max(run.baselineMaxXi, xiKin);
                    run.baselineMaxBeta = std::max(run.baselineMaxBeta, std::abs(betaKin));
                    run.baselineMaxBetaDeviation = std::max(run.baselineMaxBetaDeviation,
                                                            std::abs(betaKin - betaDesign));
                    run.baselineFinalXi = xiKin;
                    run.baselineFinalEtaDot = stKin.etaDot;
                    run.baselineMaxDqFrac = std::max(run.baselineMaxDqFrac,
                                                     qdKin.cwiseAbs().cwiseQuotient(dqMax).maxCoeff());
                    run.baselineMaxDdqFrac = std::max(run.baselineMaxDdqFrac,
                                                      qddKin.cwiseAbs().cwiseQuotient(ddqMax).maxCoeff());
                    qKin += dt * qdKin + (0.5 * dt * dt) * qddKin;
                    qdKin += dt * qddKin;
                } else {
                    baselineOk = false;
                    run.baselineFinalXi = std::nan("");
                }
            }

            if (options.keepLog) {
                MujocoLogRow r;
                r.t = t;
                r.eta = st.eta;
                r.etaDot = st.etaDot;
                r.xi = st.xi;
                r.xiNorm = xiNorm;
                r.xiDesign = xiDesign;
                r.beta = st.beta;
                r.betaDesign = betaDesign;
                r.q = q;
                r.qd = qd;
                r.tau = tau;
                r.G = G;
                r.qdd = qdd;
                r.sigmaMin = st.sigmaMin;
                r.jointMargin = hm.jointMargin;
                r.flangeSpeed = hm.flangeSpeed;
                r.tauRateCmd = hm.tauRateCmd;
                r.tauRate = hm.tauRate;
                r.contacts = contacts;
                r.constraintRows = rows;
                r.ramp = st.ramp;
                r.xiNormKin = xiKin;
                r.betaKin = betaKin;
                r.qErrKin = qErr;
                run.log.push_back(r);
            }
            ++run.steps;
        }
        run.wallTime = std::chrono::duration<double>(std::chrono::steady_clock::now() - wallStart).count();
        run.warnings = plant.warnings();
        if (run.warnings > 0)
            run.ok = false;

        // Distance from the lift at the end, after the last period
        liftError(example, controller, plant.q(), run.finalLiftError, run.finalLiftError2);
        if (options.baseline) {
            if (baselineOk)
                liftError(example, kinematic, qKin, run.baselineFinalLiftError, run.baselineFinalLiftError2);
            else
                run.baselineFinalLiftError = run.baselineFinalLiftError2 = std::nan("");
        }
        return run;
    }

}  // namespace se3pf
