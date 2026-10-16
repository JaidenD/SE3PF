"""Floquet multipliers of the far critical points (all with positive second variation): the theory predicts
real positive multipliers (lambda, 1/lambda), trace > 2, for every nondegenerate minimizer."""
import json, os
import numpy as np
from dae import poincare, monodromy, uw_from_state
from hl import spectral_deriv
from pandakin import L, DATA
res = {}
for met, c in (("I", -2.0), ("I", 2.0), ("M", -2.0), ("M", 2.0)):
    F = np.load(os.path.join(DATA, "far_%s_start%+g.npz" % (met, c)))
    q = F["q"]; dq = spectral_deriv(q, L, 1)
    uw = uw_from_state(q[0], dq[0])
    n = 1000
    for it in range(5):
        Fv = poincare(uw, met, n, q[0]) - uw
        Mo = monodromy(uw, met, n, q[0])
        if np.abs(Fv).max() < 1e-10:
            break
        uw = uw - np.linalg.solve(Mo - np.eye(2), Fv)
    ev = np.linalg.eigvals(Mo)
    res["%s_%+g" % (met, c)] = dict(shooting_residual=float(np.abs(Fv).max()), multipliers=[str(e) for e in ev],
                                    trace=float(np.trace(Mo)), det=float(np.linalg.det(Mo)))
    print(met, c, res["%s_%+g" % (met, c)], flush=True)
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_dae_far.json"), "w"), indent=1)
