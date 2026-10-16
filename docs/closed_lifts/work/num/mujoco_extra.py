"""Task 4, supplementary: MuJoCo dragged arm started on the exact periodic orbit of the ideal DAE (dae.py), compared
with that orbit over 3 laps; decomposition of the velocity closure error into its self-motion (vertical) and
constraint (flange twist) parts; growth per lap of the deviation vs the DAE Floquet multiplier.
"""
import json, os
import numpy as np
import mujoco
from mujoco_drag import build_model, run, NU
from pandakin import fk, kernel_vector, path_vel, L, DATA
from hl import trig_interp, spectral_deriv

tc = 5e-4
m = build_model(tc=tc)
D = np.load(os.path.join(DATA, "dae_orbit_M.npz"))
Q, V = D["Q"][:-1], D["V"][:-1]          # 4000 samples of the exact periodic orbit (s-time, nu = 1)
qh = np.load(os.path.join(DATA, "harm_M_1000_u0.npz"))["q"]
T1 = L / NU
res = {}
for name, q0, qd0 in (("dae_orbit_IC", Q[0], NU * V[0]),
                      ("discrete_harmonic_IC", qh[0], NU * spectral_deriv(qh, L, 1)[0])):
    r = run(m, q0, qd0, 3 * T1, lead=2 * tc)
    s = NU * r["t"]
    dev_orbit = np.linalg.norm(r["q"] - trig_interp(Q, L, s), axis=1)
    ends = [int(np.argmin(abs(r["t"] - j * T1))) for j in (1, 2, 3)]
    out = dict(dev_from_dae_orbit_max_lap1=float(dev_orbit[:ends[0] + 1].max()),
               dev_from_dae_orbit_end_of_laps=[float(dev_orbit[i]) for i in ends],
               growth_per_lap=[float(dev_orbit[ends[1]] / dev_orbit[ends[0]]), float(dev_orbit[ends[2]] / dev_orbit[ends[1]])])
    # velocity closure after one lap, decomposed
    i1 = ends[0]
    qT, qdT = r["q"][i1], r["qd"][i1]
    _, _, J = fk(qT)
    x = kernel_vector(J)
    dqd = qdT - r["qd"][0]
    out["qd_closure_total"] = float(np.linalg.norm(dqd))
    out["qd_closure_selfmotion_component_x^T dqd"] = float(x @ dqd)
    out["flange_twist_error_|J qd - nu zeta|_at_T"] = float(np.linalg.norm(J @ qdT - NU * path_vel(NU * r["t"][i1])))
    # the same quantity sampled along the lap (weld ringing)
    tw = []
    for j in range(0, i1 + 1, 50):
        _, _, Jj = fk(r["q"][j])
        tw.append(np.linalg.norm(Jj @ r["qd"][j] - NU * path_vel(NU * r["t"][j])))
    out["flange_twist_error_max_along_lap1"] = float(max(tw))
    out["flange_pos_err_max"] = float(r["ep"].max()); out["flange_att_err_max"] = float(r["ea"].max())
    res[name] = out
    print(name, json.dumps(out, indent=1), flush=True)
    np.savez(os.path.join(DATA, "mujoco_extra_%s.npz" % name), t=r["t"], dev=dev_orbit, ep=r["ep"], ea=r["ea"])
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_mujoco_extra.json"), "w"), indent=1)
