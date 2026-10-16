"""Uniqueness probe with far starts u0 = const in [-2, 2] (the linear chart is monotone on the whole admissible
fibre arc over gamma(0), u in [-2.32, 2.39]); N = 500."""
import json, os
import numpy as np
from hl import Problem
from pandakin import MjDyn, DATA
dyn = MjDyn()
res = {}
for met in ("I", "M"):
    P = Problem(500, met, dyn=dyn if met == "M" else None)
    ref = np.load(os.path.join(DATA, "harm_%s_500_u0.npz" % met))
    for c in (-2.0, -1.5, -1.0, 1.0, 1.5, 2.0):
        try:
            u0 = c + 0 * P.s
            q0 = P.chart(u0)
            cq0 = float(P.chart_quality(q0))
            u, q, E, hist = P.solve(u0, tol=8e-9, maxit=80)
            out = dict(start=c, chart_quality_at_start=cq0, E=E, iterations=len(hist) - 1,
                       max_du_vs_reference=float(np.abs(u - ref["u"]).max()), gradmax_over_h=hist[-1][2])
        except Exception as e:
            out = dict(start=c, failed=str(e))
        res["%s_%g" % (met, c)] = out
        print(met, json.dumps(out), flush=True)
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_probe_starts.json"), "w"), indent=1)
