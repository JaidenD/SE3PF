"""Task 4: the passive, gravity-free Panda whose flange is dragged along gamma(nu t) by a weld to a mocap body.

Model: menagerie panda_nohand.xml with the actuators (position servos) and the keyframe removed, joint damping 0,
gravity 0, contacts disabled, armature 0.1 kept, joint limits kept (we check they never engage). A massless body
'flange' at link7 + 0.107 m along z (the DH flange; checked), a mocap body 'target', and a weld equality
flange <-> target with relpose = identity, solref = (tc, 1), solimp = (0.9999, 0.9999, 0.001).
MuJoCo treats the mocap body as having zero velocity inside the constraint law a = -b v - k r, so a target moving at
speed V is tracked with a steady lag r = (b/k) V = 2 tc V: a pure time delay of 2 tc. We either leave it
(lead = 0) or compensate it by feeding the mocap gamma(nu (t + lead)), lead = 2 tc; both are reported.
The mocap is held during a step, so it is set to the step's midpoint time.
"""
import json, os, re, sys, time
import numpy as np
import mujoco
from pandakin import fk, path_pose, path_vel, pose_error, L, DATA, MJ_XML, kernel_vector, load_rest, logSO3
from hl import spectral_deriv, trig_interp

NU = 0.08
HERE = os.path.dirname(os.path.abspath(__file__))


def build_model(tc=5e-4, dt=1e-4, integrator="RK4"):
    xml = open(MJ_XML).read()
    assets = os.path.join(os.path.dirname(MJ_XML), "assets")
    xml = xml.replace('meshdir="assets"', 'meshdir="%s"' % assets)
    xml = xml.replace('<option integrator="implicitfast" impratio="10"/>',
                      '<option integrator="%s" timestep="%g" gravity="0 0 0"><flag contact="disable"/></option>'
                      % (integrator, dt))
    xml = xml.replace('damping="1"', 'damping="0"')
    xml = re.sub(r"<actuator>.*?</actuator>", "", xml, flags=re.S)
    xml = re.sub(r"<keyframe>.*?</keyframe>", "", xml, flags=re.S)
    xml = xml.replace('<body name="attachment"', '<body name="flange" pos="0 0 0.107"/>\n<body name="attachment"')
    xml = xml.replace("<worldbody>", '<worldbody>\n<body name="target" mocap="true"/>')
    xml = xml.replace("</mujoco>", '<equality><weld name="drag" body1="flange" body2="target" relpose="0 0 0 1 0 0 0" '
                      'solref="%g 1" solimp="0.9999 0.9999 0.001 0.5 2"/></equality>\n</mujoco>' % tc)
    m = mujoco.MjModel.from_xml_string(xml)
    return m


def quat_of(R):
    qt = np.zeros(4)
    mujoco.mju_mat2Quat(qt, np.ascontiguousarray(R).reshape(9))
    return qt


def run(m, q0, qd0, T, lead=0.0, log_every=10):
    d = mujoco.MjData(m)
    d.qpos[:] = q0
    d.qvel[:] = qd0
    dt = m.opt.timestep
    nsteps = int(round(T / dt))
    tm = (np.arange(nsteps) + 0.5) * dt + lead
    Rm, pm = path_pose(NU * tm)
    fl = m.body("flange").id
    logs_t, logs_q, logs_qd, logs_ep, logs_ea = [], [], [], [], []
    limit_hits = 0
    for k in range(nsteps + 1):
        if k % log_every == 0 or k == nsteps:
            t = k * dt
            mujoco.mj_kinematics(m, d)
            Rg, pg = path_pose(NU * t)
            Rf = d.xmat[fl].reshape(3, 3)
            logs_t.append(t); logs_q.append(d.qpos.copy()); logs_qd.append(d.qvel.copy())
            logs_ep.append(np.linalg.norm(d.xpos[fl] - pg))
            logs_ea.append(np.linalg.norm(logSO3(Rf.T @ Rg)))
        if k == nsteps:
            break
        d.mocap_pos[0] = pm[k]
        d.mocap_quat[0] = quat_of(Rm[k])
        mujoco.mj_step(m, d)
        if k % 1000 == 0 and d.nefc > 6:  # more than the 6 weld rows: a joint limit is active
            limit_hits += 1
    return dict(t=np.array(logs_t), q=np.array(logs_q), qd=np.array(logs_qd), ep=np.array(logs_ep),
                ea=np.array(logs_ea), limit_hits=limit_hits)


def report(name, r, qref, qref_d, laps):
    """Deviation from a reference lift (trig-interpolated node values qref, derivative qref_d) along the run."""
    s = NU * r["t"]
    qr = trig_interp(qref, L, s)
    dev = np.linalg.norm(r["q"] - qr, axis=1)
    T1 = L / NU
    i1 = int(np.argmin(abs(r["t"] - T1)))
    out = dict(name=name, laps=laps, max_dev_lap1=float(dev[:i1 + 1].max()), dev_end_lap1=float(dev[i1]),
               max_dev_all=float(dev.max()), dev_end=float(dev[-1]),
               closure_q_lap1=float(np.linalg.norm(r["q"][i1] - r["q"][0])),
               closure_qd_lap1=float(np.linalg.norm(r["qd"][i1] - r["qd"][0])),
               flange_pos_err_max=float(r["ep"].max()), flange_pos_err_rms=float(np.sqrt(np.mean(r["ep"] ** 2))),
               flange_att_err_max=float(r["ea"].max()), limit_hits=int(r["limit_hits"]))
    if laps > 1:
        out["dev_end_of_each_lap"] = [float(dev[int(np.argmin(abs(r["t"] - j * T1)))]) for j in range(1, laps + 1)]
    return out, dev


if __name__ == "__main__":
    tc = float(sys.argv[1]) if len(sys.argv) > 1 else 5e-4
    res = {}
    HM = np.load(os.path.join(DATA, "harm_M_1000_u0.npz"))
    qh = HM["q"]; dqh = spectral_deriv(qh, L, 1)
    _, qrest, _ = load_rest(1000)
    dqr = spectral_deriv(qrest, L, 1)
    m = build_model(tc=tc)
    # sanity: at q = q_h(0) and target = gamma(0) the weld residual vanishes
    d = mujoco.MjData(m); d.qpos[:] = qh[0]
    R0, p0 = path_pose(0.0)
    d.mocap_pos[0] = p0; d.mocap_quat[0] = quat_of(R0)
    mujoco.mj_forward(m, d)
    res["weld_residual_at_start"] = float(np.abs(d.efc_pos[:6]).max())
    print("weld residual at a consistent start: %.2e (nefc %d)" % (res["weld_residual_at_start"], d.nefc))
    T1 = L / NU
    dae_path = os.path.join(DATA, "dae_orbit_M.npz")
    runs = [("harmonic_M_lead0", qh[0], NU * dqh[0], 0.0, 1, qh),
            ("harmonic_M_lead2tc", qh[0], NU * dqh[0], 2 * tc, 3, qh),
            ("rest_lift_lead2tc", qrest[0], NU * dqr[0], 2 * tc, 1, qrest)]
    _, _, J0 = fk(qh[0])
    xk = kernel_vector(J0)
    runs.append(("harmonic_M_perturbed_w+0.01_lead2tc", qh[0], NU * (dqh[0] + 0.01 * xk), 2 * tc, 3, qh))
    if os.path.exists(dae_path):
        Dd = np.load(dae_path)
        runs.append(("dae_periodic_orbit_IC_lead2tc", Dd["Q"][0], NU * Dd["V"][0], 2 * tc, 3, qh))
    devs = {}
    for name, q0, qd0, lead, laps, ref in runs:
        t0 = time.time()
        r = run(m, q0, qd0, laps * T1, lead=lead)
        out, dev = report(name, r, ref, None, laps)
        if name.startswith("rest"):
            o2, dev2 = report(name, r, qh, None, laps)
            out["max_dev_from_harmonic_lap1"] = o2["max_dev_lap1"]
            out["dev_from_harmonic_end"] = o2["dev_end"]
        out["wall_s"] = time.time() - t0
        res[name] = out
        devs[name] = (r["t"], dev, r["ep"], r["ea"])
        print(json.dumps(out), flush=True)
    res["tc"] = tc; res["dt"] = m.opt.timestep; res["integrator"] = "RK4"
    json.dump(res, open(os.path.join(HERE, "results_mujoco_tc%g.json" % tc), "w"), indent=1)
    np.savez(os.path.join(DATA, "mujoco_tc%g.npz" % tc), **{k + "_t": v[0] for k, v in devs.items()},
             **{k + "_dev": v[1] for k, v in devs.items()}, **{k + "_ep": v[2] for k, v in devs.items()},
             **{k + "_ea": v[3] for k, v in devs.items()})
