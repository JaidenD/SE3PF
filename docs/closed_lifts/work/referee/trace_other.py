"""Referee: continue the two other admissible self-motion arcs found over gamma(0) around the loop (m nodes) and save
them in the same format as arcs_m50_d0.01.npy (k0 = index of a reference point), so dp2.py can run on them."""
import sys
import numpy as np, rk
m, d, which = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
h = rk.L/m
def inbox(q): return np.all(q > rk.QMIN) and np.all(q < rk.QMAX)
def trace_dir(q0, s, sign):
    pts = []; q = q0.copy(); x = rk.kernel(rk.jb(q))*sign
    for it in range(3000):
        qn = rk.project(q + d*x, s, x=x); xn = rk.kernel(rk.jb(qn)); xn = xn*np.sign(xn@x)
        if not inbox(qn): return pts, qn
        pts.append(qn); q = qn; x = xn
    raise RuntimeError("no exit")
arcs0 = np.load("othercomp_arcs_s0.npy", allow_pickle=True)
A0 = arcs0[which]; qref = A0[len(A0)//2]; xref = rk.kernel(rk.jb(qref))
out = []
for j in range(m):
    s = j*h
    q0 = rk.project(qref, s)          # least-squares Newton onto the new fibre from the previous reference point
    x0 = rk.kernel(rk.jb(q0)); x0 = x0*np.sign(x0@xref)
    f, ef = trace_dir(q0, s, +1); b, eb = trace_dir(q0, s, -1)
    # orient so that +1 is along x0
    if f and (f[0]-q0)@x0 < 0: f, b, ef, eb = b, f, eb, ef
    Q = np.array(b[::-1]+[q0]+f); k0 = len(b)
    # re-centre the reference at the arc midpoint so the next node starts in the middle of the arc
    qref = Q[len(Q)//2]; xref = x0
    out.append(dict(Q=Q, k0=k0, x0=x0, exitf=ef, exitb=eb))
    if j % 10 == 0:
        mg = lambda e: np.argmin(np.minimum(e-rk.QMIN, rk.QMAX-e))+1
        print(f"node {j}: arc length {d*(len(Q)-1):.2f}, ends at joints {mg(eb)},{mg(ef)}, min|q-q_r| {np.min(np.linalg.norm(Q-rk.QR,axis=1)):.2f}", flush=True)
# closure check: does the arc at node m (= node 0) coincide with the arc at node 0?
qe = rk.project(qref, 0.0); A = out[0]["Q"]
print("closure: distance from continued reference to the node-0 arc:", np.min(np.linalg.norm(A-qe, axis=1)))
np.save(f"arcs_other{which}_m{m}_d{d}.npy", np.array(out, dtype=object), allow_pickle=True)
