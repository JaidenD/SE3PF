// MuJoCo smoke test for the Panda simulation (an app, not part of the library).
//   1. Load a MuJoCo Menagerie Panda model through the C API and print what the closed-loop simulation relies on.
//   2. Step 1000 times with zero control and print the final qpos.
//   3. Torque mode (panda_nohand.xml, nv = 7): turn the position servos into torque motors and check the two facts
//      the computed-torque simulation needs:
//        a. tau = qfrc_bias holds the arm at rest at the rest posture q_r for 5 s (gravity compensation);
//        b. tau = M qdd* + qfrc_bias - qfrc_passive gives qacc = qdd* (exact model, eq:fb).
// Usage:  mujoco_smoke <path/to/mujoco_menagerie/franka_emika_panda/panda_nohand.xml>
#include <mujoco/mujoco.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {

    constexpr int kArmDof = 7;
    constexpr double kPi = 3.14159265358979323846;
    // Rest posture q_r = (0, -pi/4, 0, -3pi/4, 0, pi/2, pi/4) (design doc, Conventions).
    const double kRest[kArmDof] = {0.0, -kPi / 4, 0.0, -3 * kPi / 4, 0.0, kPi / 2, kPi / 4};

    const char* integratorName(int integrator)
    {
        switch (integrator) {
            case mjINT_EULER: return "Euler (semi-implicit)";
            case mjINT_RK4: return "RK4";
            case mjINT_IMPLICIT: return "implicit";
            case mjINT_IMPLICITFAST: return "implicitfast";
            default: return "other";
        }
    }

    const char* nameOr(const char* name) { return name ? name : "(unnamed)"; }

    double maxAbsDiff(const mjtNum* a, const double* b, int n)
    {
        double err = 0.0;
        for (int i = 0; i < n; ++i)
            err = std::max(err, std::abs(a[i] - b[i]));
        return err;
    }

    double maxAbs(const mjtNum* a, int n)
    {
        double v = 0.0;
        for (int i = 0; i < n; ++i)
            v = std::max(v, std::abs(a[i]));
        return v;
    }

    void printVector(const char* label, const mjtNum* v, int n)
    {
        std::printf("%s", label);
        for (int i = 0; i < n; ++i)
            std::printf(" %+.6f", v[i]);
        std::printf("\n");
    }

    void printModel(const mjModel* m)
    {
        std::printf("MuJoCo %s\n", mj_versionString());
        std::printf("nq = %d, nv = %d, nu = %d, nbody = %d, neq = %d\n", static_cast<int>(m->nq),
                    static_cast<int>(m->nv), static_cast<int>(m->nu), static_cast<int>(m->nbody),
                    static_cast<int>(m->neq));
        std::printf("timestep = %g s, integrator = %s, gravity = (%g, %g, %g)\n", m->opt.timestep,
                    integratorName(m->opt.integrator), m->opt.gravity[0], m->opt.gravity[1], m->opt.gravity[2]);
        for (int j = 0; j < m->njnt; ++j) {
            const int dof = m->jnt_dofadr[j];
            std::printf("  joint %d %-14s qposadr %d dofadr %d range [%+.4f, %+.4f] armature %g damping %g "
                        "frictionloss %g\n",
                        j, nameOr(mj_id2name(m, mjOBJ_JOINT, j)), m->jnt_qposadr[j], dof, m->jnt_range[2 * j],
                        m->jnt_range[2 * j + 1], m->dof_armature[dof], m->dof_damping[dof], m->dof_frictionloss[dof]);
        }
        for (int a = 0; a < m->nactuator; ++a) {
            const char* target = m->actuator_trntype[a] == mjTRN_JOINT
                                     ? nameOr(mj_id2name(m, mjOBJ_JOINT, m->actuator_trnid[2 * a]))
                                     : "(not a joint)";
            const int c = m->actuator_ctrladr[a];
            std::printf("  actuator %d %-10s on %-14s gain %g bias (%g, %g, %g) ctrlrange [%g, %g] "
                        "forcerange [%g, %g]\n",
                        a, nameOr(mj_id2name(m, mjOBJ_ACTUATOR, a)), target, m->actuator_gainprm[a * mjNGAIN],
                        m->actuator_biasprm[a * mjNBIAS], m->actuator_biasprm[a * mjNBIAS + 1],
                        m->actuator_biasprm[a * mjNBIAS + 2], m->actuator_ctrlrange[2 * c],
                        m->actuator_ctrlrange[2 * c + 1], m->actuator_forcerange[2 * a],
                        m->actuator_forcerange[2 * a + 1]);
        }
    }

    // Turn the Menagerie position servos (force = kp (ctrl - q) - kd qdot) into torque motors: force = ctrl.
    // The control is clamped to the actuator's force range, which Menagerie sets to the Franka torque limits.
    // Returns false unless actuator i drives joint i (a hinge) through control i, force output i and dof i, with unit
    // gear, a force limit and no actuator dynamics, damping, armature or delay, so that ctrl[i] = tau_i exactly.
    // Keep in sync with useTorqueMotors in sim/mujoco_panda.cpp (this app avoids se3pf on purpose, so it has a copy).
    bool useTorqueMotors(mjModel* m)
    {
        if (m->nactuator < kArmDof || m->njnt < kArmDof || m->nv < kArmDof || m->nu < kArmDof)
            return false;
        for (int i = 0; i < kArmDof; ++i) {
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

    // The flange (Franka frame F, the tool frame of the design doc) is not a body or site of the Menagerie model:
    // it is body link7 moved 0.107 m along its own z-axis, with the same orientation. (The attachment_site of
    // panda_nohand.xml sits at the flange but is turned 135 deg about z; the hand body of panda.xml -45 deg.)
    void printFlange(const mjModel* m, mjData* d, const char* label, const double* q)
    {
        const int link7 = mj_name2id(m, mjOBJ_BODY, "link7");
        if (link7 < 0)
            return;
        for (int i = 0; i < kArmDof; ++i)
            d->qpos[i] = q[i];
        mj_kinematics(m, d);
        const mjtNum* R = d->xmat + 9 * link7;  // row-major 3 x 3; column 2 is the z-axis
        double p[3] = {0.0, 0.0, 0.0};
        for (int r = 0; r < 3; ++r)
            p[r] = d->xpos[3 * link7 + r] + 0.107 * R[3 * r + 2];  // PandaArm::kFlangeOffset (not used here)
        std::printf("flange at %s: p = (%.6f, %.6f, %.6f), "
                    "R rows (%+.6f %+.6f %+.6f) (%+.6f %+.6f %+.6f) (%+.6f %+.6f %+.6f)\n",
                    label, p[0], p[1], p[2], R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8]);
    }

    // Check a: hold the rest posture with tau = qfrc_bias - qfrc_passive for `seconds`. With compensate = false the
    // torque is zero, which shows that the check can fail (the arm falls). Returns max_t |q(t) - q_r|.
    double holdRest(const mjModel* m, mjData* d, double seconds, bool compensate, double* finalSpeed)
    {
        mj_resetData(m, d);
        for (int i = 0; i < kArmDof; ++i)
            d->qpos[i] = kRest[i];
        const int steps = static_cast<int>(std::lround(seconds / m->opt.timestep));
        double drift = 0.0;
        for (int k = 0; k < steps; ++k) {
            mj_step1(m, d);  // position- and velocity-dependent terms: M, qfrc_bias, qfrc_passive
            for (int i = 0; i < kArmDof; ++i)
                d->ctrl[i] = compensate ? d->qfrc_bias[i] - d->qfrc_passive[i] : 0.0;
            mj_step2(m, d);  // actuation, acceleration, integration
            drift = std::max(drift, maxAbsDiff(d->qpos, kRest, kArmDof));
        }
        *finalSpeed = maxAbs(d->qvel, kArmDof);
        return drift;
    }

    // Check b: tau = M qdd* + qfrc_bias - qfrc_passive (eq:fb) for a constant qdd*, from a moving start.
    // Prints max |qacc - qdd*| after one mj_forward and over `steps` steps of mj_step1 / mj_step2; returns the
    // larger of the two.
    double constantAcceleration(const mjModel* m, mjData* d, int steps)
    {
        const double qdd[kArmDof] = {0.5, -0.4, 0.3, 0.6, -0.8, 1.0, -1.2};
        const double v0[kArmDof] = {0.2, -0.1, 0.15, 0.1, -0.2, 0.25, -0.3};
        double tau[kArmDof] = {0.0};
        mj_resetData(m, d);
        for (int i = 0; i < kArmDof; ++i) {
            d->qpos[i] = kRest[i];
            d->qvel[i] = v0[i];
        }

        // One evaluation: forward dynamics with the computed torque must return qdd* (no active constraints).
        mj_forward(m, d);
        mj_mulM(m, d, tau, qdd);  // M qdd*, with M including the armature
        for (int i = 0; i < kArmDof; ++i)
            d->ctrl[i] = tau[i] + d->qfrc_bias[i] - d->qfrc_passive[i];
        mj_forward(m, d);
        double accErr = maxAbsDiff(d->qacc, qdd, kArmDof);
        std::printf("    one mj_forward: max |qacc - qdd*| = %.3e (nefc = %d, ncon = %d)\n", accErr,
                    static_cast<int>(d->nefc), static_cast<int>(d->ncon));

        // Closed loop at the control rate: recompute tau every step, as the simulation will.
        double maxTau = 0.0;
        for (int k = 0; k < steps; ++k) {
            mj_step1(m, d);
            mj_mulM(m, d, tau, qdd);
            for (int i = 0; i < kArmDof; ++i)
                d->ctrl[i] = tau[i] + d->qfrc_bias[i] - d->qfrc_passive[i];
            maxTau = std::max(maxTau, maxAbs(d->ctrl, kArmDof));
            mj_step2(m, d);
            accErr = std::max(accErr, maxAbsDiff(d->qacc, qdd, kArmDof));
        }
        // Exact solution for a constant acceleration. MuJoCo's Euler and implicit integrators update qvel first and
        // then qpos with the new qvel, which adds h t qdd* / 2 to the position: an O(h) error, not a model error.
        const double t = steps * m->opt.timestep;
        double velErr = 0.0, posErr = 0.0, eulerBias = 0.0;
        for (int i = 0; i < kArmDof; ++i) {
            velErr = std::max(velErr, std::abs(d->qvel[i] - (v0[i] + t * qdd[i])));
            posErr = std::max(posErr, std::abs(d->qpos[i] - (kRest[i] + t * v0[i] + 0.5 * t * t * qdd[i])));
            eulerBias = std::max(eulerBias, std::abs(0.5 * m->opt.timestep * t * qdd[i]));
        }
        std::printf("    %d steps (%.3f s): max |qacc - qdd*| = %.3e, max |qvel - (v0 + t qdd*)| = %.3e\n", steps, t,
                    accErr, velErr);
        std::printf("    max |qpos - exact| = %.3e (Euler's h t |qdd*| / 2 = %.3e), max |tau| = %.2f N m\n", posErr,
                    eulerBias, maxTau);
        return accErr;
    }

}  // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <path/to/franka_emika_panda/panda_nohand.xml>\n", argv[0]);
        return 2;
    }
    char error[1000] = "";
    mjModel* m = mj_loadXML(argv[1], nullptr, error, sizeof(error));
    if (!m) {
        std::fprintf(stderr, "could not load %s: %s\n", argv[1], error);
        return 1;
    }
    mjData* d = mj_makeData(m);

    // 1. The model
    printModel(m);
    const double zero[kArmDof] = {0.0};
    printFlange(m, d, "q = 0", zero);  // design doc: (0.088, 0, 0.926), R = diag(1, -1, -1)
    printFlange(m, d, "q_r  ", kRest);  // design doc: (0.306891, 0, 0.590282)

    // 2. 1000 steps with zero control. The Menagerie actuators are position servos, so ctrl = 0 asks for q = 0,
    //    clamped to each actuator's ctrlrange (joint 4 cannot reach 0).
    mj_resetData(m, d);
    printVector("qpos(0)    :", d->qpos, static_cast<int>(m->nq));
    for (int k = 0; k < 1000; ++k)
        mj_step(m, d);
    std::printf("after 1000 steps, t = %.3f s:\n", d->time);
    printVector("qpos(t)    :", d->qpos, static_cast<int>(m->nq));

    // 3. Torque mode. The exit code reports the exact-model pass: 0 if it holds the rest posture, the control
    //    without torque falls, and qacc = qdd* to 1e-10.
    bool ok = true;
    if (m->nv != kArmDof || !useTorqueMotors(m)) {
        std::printf("torque checks skipped: they need panda_nohand.xml (nv = 7, actuator i on joint i)\n");
    } else {
        m->opt.timestep = 0.001;  // the 1 kHz control period
        mj_resetData(m, d);
        for (int i = 0; i < kArmDof; ++i)
            d->qpos[i] = kRest[i];
        mj_forward(m, d);
        printVector("G(q_r)     :", d->qfrc_bias, kArmDof);
        mjtNum M[kArmDof * kArmDof] = {0.0};
        mj_fullM(m, d, M);  // dense, row-major (MuJoCo 3.15 keeps M sparse in d->M; there is no d->qM any more)
        std::printf("diag M(q_r): ");
        for (int i = 0; i < kArmDof; ++i)
            std::printf(" %.6f", M[i * kArmDof + i]);
        std::printf("  (includes armature %g)\n", m->dof_armature[0]);

        for (int pass = 0; pass < 2; ++pass) {
            // Pass 0: the model's joint damping (1 N m s/rad), compensated through qfrc_passive. qacc matches qdd*,
            //         but implicitfast (and Euler) integrate damping implicitly, (M + h D) dv = h f, so the
            //         compensated motion is not exact: qvel is off by about 1e-2 rad/s after 1 s.
            // Pass 1: damping switched off (mjDSBL_DAMPER): the exact-model setting.
            if (pass == 1)
                m->opt.disableflags |= mjDSBL_DAMPER;
            std::printf("torque motors, timestep %g s, %s, joint damping %s:\n", m->opt.timestep,
                        integratorName(m->opt.integrator), pass == 0 ? "on (compensated)" : "off");
            double speed = 0.0;
            const double drift = holdRest(m, d, 5.0, true, &speed);
            std::printf("    hold q_r for 5 s with tau = qfrc_bias - qfrc_passive: max |q - q_r| = %.3e rad, "
                        "final |qvel| = %.3e\n",
                        drift, speed);
            const double fall = holdRest(m, d, 0.1, false, &speed);
            std::printf("    control, tau = 0 for 0.1 s: max |q - q_r| = %.3e rad\n", fall);
            const double accErr = constantAcceleration(m, d, 1000);
            if (pass == 1)
                ok = drift < 1e-10 && fall > 1e-3 && accErr < 1e-10;
        }
        std::printf("torque checks (damping off): %s\n", ok ? "PASS" : "FAIL");
    }

    mj_deleteData(d);
    mj_deleteModel(m);
    return ok ? 0 : 1;
}
