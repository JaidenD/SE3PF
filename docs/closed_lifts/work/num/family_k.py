"""Task 5: minimizers of E_k = E + (k/2) oint |q - q_r|^2 ds (Euclidean penalty) for K = I and K = M.

Continuation in k upward from the harmonic lift (k = 0) and downward from the rest-posture lift (start u = 0 at the
largest k), N = 500; the k = 100 and 1e4 cases are repeated at N = 1000.
"""
import json, os
import numpy as np
from hl import Problem
from pandakin import MjDyn, DATA, L

dyn = MjDyn()
ks = [0.0, 1.0, 10.0, 100.0, 1e3, 1e4]
res = {}
for met in ("I", "M"):
    sol_up = {}
    for direction in ("up", "down"):
        order = ks if direction == "up" else ks[::-1]
        u = None
        for k in order:
            P = Problem(500, met, kpen=k, dyn=dyn if met == "M" else None)
            if u is None:
                u = np.zeros(P.N) if direction == "down" else np.load(os.path.join(DATA, "harm_%s_500_u0.npz" % met))["u"]
            u, q, Ek, hist = P.solve(u, tol=8e-9, maxit=60)
            P0 = Problem(500, met, kpen=0.0, dyn=P.dyn)
            E = P0.energy(q)
            harm = np.load(os.path.join(DATA, "harm_%s_500_u0.npz" % met))["q"]
            H = P.hessian_u(u, q)
            lam = np.linalg.eigvalsh(H)[0]
            row = dict(k=k, Ek=Ek, E=E, penalty=Ek - E, umin=float(u.min()), umax=float(u.max()),
                       max_dist_rest=float(np.linalg.norm(q - P.qrest, axis=1).max()),
                       max_dist_harm=float(np.linalg.norm(q - harm, axis=1).max()),
                       lam_min_plain=float(lam), iterations=len(hist) - 1)
            key = "%s_N500_k%g_%s" % (met, k, direction)
            res[key] = row
            if direction == "up":
                sol_up[k] = u.copy()
                np.savez(os.path.join(DATA, "family_%s_k%g.npz" % (met, k)), u=u, q=q, s=P.s)
            else:
                row["max_|u_down - u_up|"] = float(np.abs(u - sol_up[k]).max())
            print("%s %-4s k=%-6g E=%.9f pen=%.3e u in [%.4f, %.4f] dist(rest)=%.3e dist(harm)=%.3e lam_min=%.3e it=%d %s" % (
                met, direction, k, E, Ek - E, u.min(), u.max(), row["max_dist_rest"], row["max_dist_harm"], lam,
                row["iterations"], ("|du vs up| = %.1e" % row["max_|u_down - u_up|"]) if direction == "down" else ""),
                flush=True)
    for k in (100.0, 1e4):
        P = Problem(1000, met, kpen=k, dyn=dyn if met == "M" else None)
        u, q, Ek, hist = P.solve(np.zeros(1000), tol=3e-8, maxit=60)
        res["%s_N1000_k%g" % (met, k)] = dict(k=k, max_dist_rest=float(np.linalg.norm(q - P.qrest, axis=1).max()),
                                             umin=float(u.min()), umax=float(u.max()))
        print("%s N=1000 k=%g dist(rest)=%.4e u in [%.4f, %.4f]" % (met, k, res["%s_N1000_k%g" % (met, k)]["max_dist_rest"],
                                                                     u.min(), u.max()), flush=True)
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_family.json"), "w"), indent=1)
