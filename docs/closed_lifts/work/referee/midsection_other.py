"""Referee: midsection of the OTHER admissible components (traced by trace_other.py, m = 50, step 0.01): bounding joints
at the arc ends per node, and the s-derivative of the midsection's arclength position (kinks at face switches).
Arclength position is measured from a fixed reference: we track the midsection point's distance along the arc from
the arc point closest to the previous node's midsection."""
import numpy as np, rk
def fd(y): return np.roll(y, -1)-y
for w in (0, 1):
    arcs = np.load(f"arcs_other{w}_m50_d0.01.npy", allow_pickle=True); m = len(arcs); h = rk.L/m
    JB, JF, QM, Ls, A, B = [], [], [], [], [], []
    for a in arcs:
        Q = a["Q"]; full = np.vstack([a["exitb"], Q, a["exitf"]])
        def end(qin, qout):
            best, jb = 1.0, -1
            for j in range(7):
                for lim, sg in ((rk.QMIN[j], -1), (rk.QMAX[j], 1)):
                    if (qout[j]-lim)*sg > 0 and (qin[j]-lim)*sg < 0:
                        f = (lim-qin[j])/(qout[j]-qin[j])
                        if f < best: best, jb = f, j
            return qin+best*(qout-qin), jb
        qb, jb = end(Q[0], a["exitb"]); qf, jf = end(Q[-1], a["exitf"])
        P = np.vstack([qb, Q, qf]); c = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
        half = c[-1]/2; i = np.searchsorted(c, half)-1; f = (half-c[i])/(c[i+1]-c[i]); qm = P[i]+f*(P[i+1]-P[i])
        JB.append(jb+1); JF.append(jf+1); QM.append(qm); Ls.append(c[-1])
    QM = np.array(QM); Ls = np.array(Ls)
    # arclength speed of the midsection along the self-motion: project node-to-node displacement on the kernel vector
    X = np.array([rk.kernel(rk.jb(q)) for q in QM])
    for j in range(1, m):
        if X[j]@X[j-1] < 0: X[j] = -X[j]
    v = np.array([X[j]@(QM[(j+1) % m]-QM[j]) for j in range(m)])/h   # self-motion speed of the midsection
    sw = [j for j in range(m) if JB[j] != JB[(j+1) % m] or JF[j] != JF[(j+1) % m]]
    print(f"component {w}: lower-end joints {sorted(set(JB))}, upper-end joints {sorted(set(JF))}, switches after nodes {sw}; "
          f"arc length [{Ls.min():.2f},{Ls.max():.2f}]")
    print("   self-motion speed of the midsection x.dq/ds per node:", np.round(v, 2).tolist())
    print(f"   max |second difference of that speed|/h: {np.abs(fd(v)).max()/h:.1f} (median {np.median(np.abs(fd(v)))/h:.1f})")
