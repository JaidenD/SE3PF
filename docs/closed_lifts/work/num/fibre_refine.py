"""Refine the homology-class change of the T^7 fibre between s = 0.2310 and 0.2503: fine step h = 0.002,
fibre traced from the rest-lift configuration over gamma(s) (by constrained IK from the neighbouring node)."""
import json, os
import numpy as np
from pandakin import fk, kernel_vector, path_pose, QMIN, QMAX, L, DATA, load_rest, constrained_ik
from fibre import wrap, project
s_all, q_all, x_all = load_rest(2000)
res = []
for k0 in range(240, 262, 1):          # s = k0 L / 2000: 0.0924 ... ; we want 0.231..0.250 -> k0 = 600..650
    pass
for k0 in range(600, 652, 3):
    s0, qrest = s_all[k0], q_all[k0]
    Rt, pt = path_pose(s0)
    q = qrest.copy(); x_prev = kernel_vector(fk(q)[2]); h = 0.002; tau = 0.0; P = [q.copy()]
    for i in range(20000):
        def f(qq):
            x = kernel_vector(fk(qq)[2]); return x if x @ x_prev >= 0 else -x
        k1 = f(q); k2 = f(q + 0.5*h*k1); k3 = f(q + 0.5*h*k2); k4 = f(q + h*k3)
        q = project(q + (h/6)*(k1 + 2*k2 + 2*k3 + k4), Rt, pt); x_prev = f(q); tau += h; P.append(q.copy())
        if tau > 0.5 and np.linalg.norm(wrap(q - qrest)) < 0.6*h: break
    P = np.array(P); sv = np.linalg.svd(fk(P)[2], compute_uv=False)[:, -1]
    j = int(np.argmin(sv))
    row = dict(s=float(s0), length=float(tau), winding=np.round((P[-1]-P[0])/(2*np.pi)).astype(int).tolist(),
               sigmin_fibre=float(sv.min()), q_at_sigmin=np.round(wrap(P[j]), 3).tolist())
    res.append(row); print(json.dumps(row), flush=True)
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_fibre_refine.json"), "w"), indent=1)
