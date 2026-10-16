"""Referee: own stepping loop for the MuJoCo dragged arm (weld model built by num/mujoco_drag.build_model), started on
the stored K = M harmonic lift q_h(0), qdot = nu q_h'(0) (spectral derivative). Reports |q - q_h(nu t)| at the end of
laps 1..3 and the max over lap 1 (reference interpolated trigonometrically), with the 2 tc lag compensation."""
import sys, os, time
import numpy as np
sys.path.insert(0, "../num")
import mujoco
from mujoco_drag import build_model, quat_of, NU
from pandakin import path_pose, L, DATA
from hl import spectral_deriv, trig_interp
tc = 5e-4; m = build_model(tc=tc); dt = m.opt.timestep; lead = 2*tc
qh = np.load(os.path.join(DATA, "harm_M_1000_u0.npz"))["q"]; dqh = spectral_deriv(qh, L, 1)
qr = np.load("../rest_lift_1000.npy")
def lap_run(q0, qd0, laps, ref):
    d = mujoco.MjData(m); d.qpos[:] = q0; d.qvel[:] = qd0
    T = L/NU; n = int(round(T/dt)); out = []; maxdev = 0; t0 = time.time()
    for k in range(laps*n):
        tm = (k+0.5)*dt+lead; R, p = path_pose(np.array([NU*tm])); d.mocap_pos[0] = p[0]; d.mocap_quat[0] = quat_of(R[0])
        mujoco.mj_step(m, d)
        if (k+1) % n == 0: out.append(np.linalg.norm(d.qpos-ref[0]))
        if k < n and (k+1) % 200 == 0:
            s = (NU*(k+1)*dt) % L; maxdev = max(maxdev, np.linalg.norm(d.qpos-trig_interp(ref, L, [s])[0]))
    return out, maxdev, time.time()-t0
ends, md, w = lap_run(qh[0], NU*dqh[0], 3, qh)
print(f"K=M harmonic start: |q - q_h| at lap ends {['%.3e' % e for e in ends]}, ratios {[round(ends[i+1]/ends[i], 2) for i in range(len(ends)-1)]}, "
      f"max over lap 1 (every 200 steps) {md:.3e}  [{w:.0f}s]", flush=True)
ends, md, w = lap_run(qr[0], NU*spectral_deriv(qr, L, 1)[0], 1, qr)
print(f"rest-lift start: |q(T) - q_rest(0)| {ends[0]:.3e}, max deviation from rest lift over lap 1 {md:.3e}", flush=True)
