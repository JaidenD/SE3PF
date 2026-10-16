"""Search for admissible self-motion components over gamma(s) other than the arc through the rest lift."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import rk
here = os.path.dirname(os.path.abspath(__file__))
arcs = np.load(os.path.join(here, "arcs_m50_d0.01.npy"), allow_pickle=True)
rng = np.random.default_rng(1)
h = rk.L/50
for j in [0, 12, 25, 37]:
    s = j*h; arc = arcs[j]["Q"]
    found = []; n_conv = 0
    for trial in range(400):
        q = rk.QMIN + (rk.QMAX-rk.QMIN)*rng.random(7)
        for it in range(60):
            e = rk.pose_err(q, s)
            if np.linalg.norm(e) < 1e-12: break
            J = rk.jb(q)
            dq = -np.linalg.lstsq(J, e, rcond=None)[0]
            nd = np.linalg.norm(dq)
            if nd > 0.3: dq *= 0.3/nd
            q = q + dq
        if np.linalg.norm(rk.pose_err(q, s)) > 1e-10: continue
        n_conv += 1
        if not (np.all(q > rk.QMIN) and np.all(q < rk.QMAX)): continue
        dist = np.min(np.linalg.norm(arc - q, axis=1))
        if dist > 0.03:
            found.append(q)
    print(f"node {j}: converged {n_conv}/400; admissible IK solutions NOT on the rest-lift arc: {len(found)}")
    if found:
        F = np.array(found)
        # cluster
        reps = []
        for q in F:
            if all(np.linalg.norm(q-r) > 0.5 for r in reps): reps.append(q)
        for r in reps[:6]:
            print("   rep", np.round(r, 3), "sigmin", round(np.linalg.svd(rk.jb(r), compute_uv=False)[-1], 3))
        np.save(os.path.join(here, f"othercomp_{j}.npy"), F)
