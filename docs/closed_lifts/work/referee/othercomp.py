"""Referee: the other admissible self-motion components over gamma(0) (IK solutions found by random restarts that are
not on the arc through the rest lift). Trace their admissible arcs and report ends, lengths, and whether they are
distinct from the main arc."""
import numpy as np, rk
d = 0.01
def inbox(q): return np.all(q > rk.QMIN) and np.all(q < rk.QMAX)
def trace(q0, s, sign):
    pts = []; q = q0.copy(); x = rk.kernel(rk.jb(q))*sign
    for it in range(3000):
        qn = rk.project(q + d*x, s, x=x)
        xn = rk.kernel(rk.jb(qn)); xn = xn*np.sign(xn@x)
        if not inbox(qn): return pts, qn
        if len(pts) > 50 and np.linalg.norm(qn-q0) < 0.6*d: return pts, None  # closed circle
        pts.append(qn); q = qn; x = xn
    return pts, None
main = np.load("arcs_m50_d0.01.npy", allow_pickle=True)[0]["Q"]
F = np.load("othercomp_0.npy")
arcs = []
for q in F:
    if any(np.min(np.linalg.norm(A-q, axis=1)) < 0.05 for A in arcs): continue
    f, ef = trace(q, 0.0, +1); b, eb = trace(q, 0.0, -1)
    A = np.array(b[::-1]+[q]+f); arcs.append(A)
    def lim(e):
        if e is None: return "closed"
        mg = np.minimum(e-rk.QMIN, rk.QMAX-e); j = np.argmin(mg); return f"j{j+1}{'-' if e[j]-rk.QMIN[j] < rk.QMAX[j]-e[j] else '+'}"
    sig = min(np.linalg.svd(rk.jb(p), compute_uv=False)[-1] for p in A[::10])
    print(f"arc {len(arcs)}: length {d*(len(A)-1):.2f} rad, ends {lim(eb)} / {lim(ef)}, min dist to main arc "
          f"{min(np.min(np.linalg.norm(main-p, axis=1)) for p in A[::5]):.3f}, sigma_min>={sig:.3f}, "
          f"q1 range [{A[:,0].min():.2f},{A[:,0].max():.2f}] q2 range [{A[:,1].min():.2f},{A[:,1].max():.2f}]")
    # energy-relevant: minimum of |q - q_r| on arc
    print("    min |q-q_r| on arc:", round(np.min(np.linalg.norm(A-rk.QR, axis=1)), 3), " (main arc:", round(np.min(np.linalg.norm(main-rk.QR, axis=1)), 3), ")")
print("covered:", sum(any(np.min(np.linalg.norm(A-q, axis=1)) < 0.05 for A in arcs) for q in F), "of", len(F))
np.save("othercomp_arcs_s0.npy", np.array(arcs, dtype=object), allow_pickle=True)
