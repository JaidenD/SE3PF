"""Joint-limit-aware minimization: E + (kappa/2) h sum max(0, violation of the 0.10 rad safety band)^2, continuation
kappa = 1e2 ... 1e6, N = 500, started in each far basin (from the far critical points of analyze_far.py) and in the
near basin. Reports the energy E (penalty excluded) and the smallest joint margin of each limit-respecting lift."""
import json, os
import numpy as np
from hl import Problem
from pandakin import MjDyn, DATA, joint_margin
dyn = MjDyn()
res = {}
for met in ("I", "M"):
    starts = {"near": np.load(os.path.join(DATA, "harm_%s_1000_u0.npz" % met))["q"][::2]}
    for c in (-2.0, 2.0):
        starts["far%+g" % c] = np.load(os.path.join(DATA, "far_%s_start%+g.npz" % (met, c)))["q"][::2]
    for name, q0 in starts.items():
        P0 = Problem(500, met, dyn=dyn if met == "M" else None)
        u = P0.u_of_q(q0)
        for kap in (1e2, 1e3, 1e4, 1e5, 1e6):
            P = Problem(500, met, dyn=P0.dyn, lim_kappa=kap)
            try:
                u, q, Ek, hist = P.solve(u, tol=1e-7, maxit=120)
            except RuntimeError as e:
                print(met, name, kap, "failed", e); break
            E = P0.energy(q)
            mar = joint_margin(q)[0].min()
            print("%s %-7s kappa=%.0e  E=%.6f  E+pen=%.6f  min margin=%.4f  it=%d  g/h=%.1e" % (
                met, name, kap, E, Ek, mar, len(hist) - 1, hist[-1][2]), flush=True)
            res["%s_%s_k%.0e" % (met, name, kap)] = dict(E=E, Epen=Ek, min_margin=float(mar), iterations=len(hist) - 1)
        np.savez(os.path.join(DATA, "limits_%s_%s.npz" % (met, name)), u=u, q=q, s=P0.s)
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_limits.json"), "w"), indent=1)
