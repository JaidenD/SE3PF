"""Finer scan of the fibre over gamma(s) through the rest-lift node: homology class in T^7 (winding of each joint),
smallest sigma_min(J_b) on the whole fibre and on its admissible arc (inside the joint limits), arc length."""
import json, os
import numpy as np
from pandakin import fk, kernel_vector, path_pose, pose_error, QMIN, QMAX, L, DATA, load_rest
from fibre import wrap, project
res = []
s_all, q_all, _ = load_rest(1000)
for k0 in range(0, 1000, 25):
    s0, qrest = s_all[k0], q_all[k0]
    Rt, pt = path_pose(s0)
    q = qrest.copy(); x_prev = kernel_vector(fk(q)[2]); h = 0.01; tau = 0.0; P = [q.copy()]
    for i in range(4000):
        def f(qq):
            x = kernel_vector(fk(qq)[2]); return x if x @ x_prev >= 0 else -x
        k1 = f(q); k2 = f(q + 0.5*h*k1); k3 = f(q + 0.5*h*k2); k4 = f(q + h*k3)
        q = project(q + (h/6)*(k1 + 2*k2 + 2*k3 + k4), Rt, pt); x_prev = f(q); tau += h; P.append(q.copy())
        if tau > 0.5 and np.linalg.norm(wrap(q - qrest)) < 0.6*h: break
    P = np.array(P); sv = np.linalg.svd(fk(P)[2], compute_uv=False)[:, -1]
    ins = np.all((wrap(P) >= QMIN) & (wrap(P) <= QMAX), axis=1)
    n = len(P) - 1; fwd = 0
    while fwd < n and ins[fwd + 1]: fwd += 1
    bwd = 0
    while bwd < n and ins[(-bwd - 1) % n]: bwd += 1
    arc = np.r_[np.arange(0, fwd + 1), np.arange(n - bwd, n)]
    row = dict(s=float(s0), length=float(tau), winding=np.round((P[-1] - P[0]) / (2*np.pi)).astype(int).tolist(),
               sigmin_fibre=float(sv.min()), sigmin_admissible_arc=float(sv[arc].min()), arc_length=float((fwd + bwd)*h))
    res.append(row); print(json.dumps(row), flush=True)
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_fibre_scan.json"), "w"), indent=1)
