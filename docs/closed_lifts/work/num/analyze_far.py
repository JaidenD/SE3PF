"""The other critical points found from far starts (probe_starts.py): are they inside the joint limits, are they
local minima, are they grid-converged? N = 500 from u0 = +-2 (const), then refined at N = 1000 from the trig-
interpolated N = 500 solution."""
import json, os
import numpy as np
from hl import Problem, trig_interp
from pandakin import MjDyn, DATA, L, QMIN, QMAX, fk, joint_margin

dyn = MjDyn()
res = {}


def analyse(P, u, q, met):
    t, _ = P.tangents(q)
    Kt = t if met == "I" else np.einsum("kij,kj->ki", dyn.Mbatch(q), t)
    D = P.h * np.einsum("ki,ki->k", t, Kt)
    Dm = 1 / np.sqrt(D)
    H = P.hessian_u(u, q)
    w = np.linalg.eigvalsh(Dm[:, None] * H * Dm[None, :])
    mar, jj = joint_margin(q)
    per_joint = np.minimum(q - QMIN, QMAX - q).min(0)
    r, Fn, _ = P.el_residual(q)
    g, _ = P.grad_u(u)
    return dict(E=float(P.energy(q)), lam_gen=w[:3].tolist(), n_negative=int((w < 0).sum()),
                min_margin=float(mar.min()), per_joint_min_margin=per_joint.tolist(),
                fraction_nodes_inside_limits=float((mar >= 0).mean()),
                chart_quality=float(P.chart_quality(q)), EL_max=float(np.abs(r).max()),
                gradmax_over_h=float(np.abs(g).max() / P.h), umin=float(u.min()), umax=float(u.max()),
                minSigma=float(np.linalg.svd(fk(q)[2], compute_uv=False)[:, -1].min()))


for met in ("I", "M"):
    for c in (-2.0, 2.0):
        P = Problem(500, met, dyn=dyn if met == "M" else None)
        u, q, E, hist = P.solve(c + 0 * P.s, tol=8e-9, maxit=200)
        a500 = analyse(P, u, q, met)
        P2 = Problem(1000, met, dyn=dyn if met == "M" else None)
        q2g = trig_interp(q, L, P2.s)
        u2, q2, E2, hist2 = P2.solve(P2.u_of_q(q2g), tol=3e-8, maxit=200)
        a1000 = analyse(P2, u2, q2, met)
        a1000["max_dq_vs_N500_at_common_nodes"] = float(np.abs(q2[::2] - q).max())
        key = "%s_start%+g" % (met, c)
        res[key] = dict(N500=a500, N1000=a1000, iterations500=len(hist) - 1, iterations1000=len(hist2) - 1)
        np.savez(os.path.join(DATA, "far_%s.npz" % key), s=P2.s, u=u2, q=q2)
        print(key, json.dumps(res[key], indent=1), flush=True)
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_far.json"), "w"), indent=1)
