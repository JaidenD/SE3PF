"""Trace the admissible self-motion arc through the rest-posture lift at m nodes (arclength step d),
and search for other admissible fibre components by random IK starts."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import rk

HL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
m = int(sys.argv[1]) if len(sys.argv) > 1 else 50
d = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
rest = np.load(os.path.join(HL, "rest_lift_1000.npy"))
h = rk.L/m
idx = (np.arange(m)*1000)//m
assert np.allclose(idx*rk.L/1000, np.arange(m)*h)


def inbox(q):
    return np.all(q > rk.QMIN) and np.all(q < rk.QMAX)


def trace_dir(q0, s, sign, xref):
    pts = []; q = q0.copy(); x = rk.kernel(rk.jb(q)); x = x*np.sign(x@xref)*sign
    for it in range(4000):
        qp = q + d*x
        qn = rk.project(qp, s, x=x)
        xn = rk.kernel(rk.jb(qn)); xn = xn*np.sign(xn@x)
        if not inbox(qn):
            # bisect exit point roughly
            return pts, qn
        pts.append(qn); q = qn; x = xn
    return pts, None


data = {}
for j in range(m):
    s = j*h; q0 = rest[idx[j]]
    xref = rk.kernel(rk.jb(q0))
    if j > 0:
        xref = xref*np.sign(xref@data[j-1]["x0"])
    fwd, ef = trace_dir(q0, s, +1, xref)
    bwd, eb = trace_dir(q0, s, -1, xref)
    Q = np.array(bwd[::-1] + [q0] + fwd)
    k0 = len(bwd)
    data[j] = dict(Q=Q, k0=k0, x0=xref, exitf=ef, exitb=eb)
    if j % 10 == 0:
        jf = np.argmin(np.minimum(ef-rk.QMIN, rk.QMAX-ef)) if ef is not None else -1
        jbk = np.argmin(np.minimum(eb-rk.QMIN, rk.QMAX-eb)) if eb is not None else -1
        print(f"node {j}: arc t in [{-k0*d:.2f}, {len(fwd)*d:.2f}], exits at joints {jbk+1},{jf+1}; "
              f"sigmin on arc {min(np.linalg.svd(rk.jb(q),compute_uv=False)[-1] for q in Q[::20]):.3f}", flush=True)
np.save(os.path.join(os.path.dirname(__file__), f"arcs_m{m}_d{d}.npy"), np.array([data[j] for j in range(m)], dtype=object), allow_pickle=True)
