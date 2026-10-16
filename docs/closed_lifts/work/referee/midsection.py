"""Referee: the 'midsection lift' of alternatives.md on the traced admissible arcs (m = 100 nodes, step 0.005 rad).
(1) bounding joints at both ends for every node and how close the ends are to a switch of the bounding face;
(2) smoothness of the end points and of the midsection in s (periodic finite differences);
(3) dependence on the joint-space metric used for the fibre measure (Euclidean, a diagonal re-weighting, K = M);
(4) a modified joint box in which the bounding face switches: the midsection acquires a kink (jump in d/ds)."""
import sys
import numpy as np
import rk
arcs = np.load("arcs_m100_d0.005.npy", allow_pickle=True)
m = len(arcs); h = rk.L/m; d = 0.005

def ends_and_arc(a, lo_box, hi_box):
    """cut the traced arc at the first exit from the box [lo_box, hi_box] in each direction from the rest-lift node;
    returns (points incl. interpolated end points, cumulative Euclidean arclength, bounding joints)"""
    Q = list(a["Q"]); k0 = a["k0"]
    full = np.array([a["exitb"]] + Q + [a["exitf"]])   # exit points are outside the original box
    k0 += 1
    inside = np.all((full > lo_box) & (full < hi_box), axis=1)
    # walk from k0 outwards
    def walk(step):
        k = k0
        while inside[k+step]: k += step
        q_in, q_out = full[k], full[k+step]
        # fraction along segment to the first violated face
        best, jb = 1.0, -1
        for j in range(7):
            for lim, sgn in ((lo_box[j], -1), (hi_box[j], +1)):
                if (q_out[j]-lim)*sgn > 0 and (q_in[j]-lim)*sgn < 0:
                    f = (lim-q_in[j])/(q_out[j]-q_in[j])
                    if f < best: best, jb = f, j
        return k, q_in + best*(q_out-q_in), jb
    kb, qb, jbk = walk(-1); kf, qf, jfw = walk(+1)
    P = np.vstack([qb, full[kb:kf+1], qf])
    return P, jbk, jfw, full[kb], full[kf]

def midpoint(P, metric=None):
    dP = np.diff(P, axis=0)
    if metric is None: ds = np.linalg.norm(dP, axis=1)
    else:
        Mm = metric(0.5*(P[1:]+P[:-1])); ds = np.sqrt(np.einsum("ki,kij,kj->k", dP, Mm, dP))
    c = np.concatenate([[0], np.cumsum(ds)]); half = c[-1]/2
    i = np.searchsorted(c, half)-1; f = (half-c[i])/(c[i+1]-c[i])
    return P[i] + f*(P[i+1]-P[i]), c[-1]

def run(lo_box, hi_box, metric=None, label=""):
    QM, A, B, JB, JF, M2 = [], [], [], [], [], []
    for a in arcs:
        P, jbk, jfw, qbin, qfin = ends_and_arc(a, lo_box, hi_box)
        qm, Ltot = midpoint(P, metric)
        QM.append(qm); JB.append(jbk); JF.append(jfw)
        # Euclidean arclength coordinate of ends relative to the rest-lift node (along the traced arc)
        dist = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
        rest = a["Q"][a["k0"]]; ir = np.argmin(np.linalg.norm(P-rest, axis=1))
        A.append(-dist[ir]); B.append(dist[-1]-dist[ir])
        # second-closest face at each end point (how far from a switch of the bounding face)
        sec = []
        for q, jj in ((P[0], jbk), (P[-1], jfw)):
            mg = np.minimum(q-lo_box, hi_box-q); mg[jj] = np.inf; sec.append(mg.min())
        M2.append(min(sec))
    QM = np.array(QM); A = np.array(A); B = np.array(B)
    return QM, A, B, np.array(JB), np.array(JF), np.array(M2)

def fd(y, k):
    for _ in range(k): y = np.roll(y, -1)-y
    return y

# (1)-(2) original box, Euclidean arclength
QM, A, B, JB, JF, M2 = run(rk.QMIN, rk.QMAX)
print("bounding joints (lower end / upper end) over all", m, "nodes:", set(JB+1), set(JF+1))
print(f"arc ends a(s) in [{A.min():.3f},{A.max():.3f}], b(s) in [{B.min():.3f},{B.max():.3f}]; smallest margin of the "
      f"second-closest face at an end point: {M2.min():.3f} rad")
tm = 0.5*(A+B)
print(f"midsection offset from rest lift along the arc: [{tm.min():.3f},{tm.max():.3f}] rad")
for k in (1, 2, 3):
    print(f"  max |Delta^{k} t_mid|/h^{k} = {np.abs(fd(tm, k)).max()/h**k:.3g},  max |Delta^{k} a|/h^{k} = {np.abs(fd(A, k)).max()/h**k:.3g}")
mg = np.minimum(QM-rk.QMIN, rk.QMAX-QM); print(f"midsection min joint margin {mg.min():.4f} (joint {np.unravel_index(np.argmin(mg), mg.shape)[1]+1}),"
      f" E_I = {np.sum((np.roll(QM,-1,0)-QM)**2)/(2*h):.3f} (m={m})")
adv = np.load("../adv/midpoint_500.npz")["Qm"]
print("  vs alternatives' midsection (500 nodes, every 5th):", np.abs(adv[::5]-QM).max())
np.save("midsection_I_m100.npy", QM)

# (3) metric dependence
wdiag = np.array([1, 1, 1, 1, 0.25, 0.25, 0.25])
QMw, *_ = run(rk.QMIN, rk.QMAX, metric=lambda X: np.broadcast_to(np.diag(wdiag), (len(X), 7, 7)))
print(f"metric diag(1,1,1,1,.25,.25,.25): midsection moves by up to {np.linalg.norm(QMw-QM, axis=1).max():.3f} rad (2-norm)")
try:
    sys.path.insert(0, "../num"); from pandakin import MjDyn
    dyn = MjDyn()
    QMm, *_ = run(rk.QMIN, rk.QMAX, metric=lambda X: dyn.Mbatch(X))
    print(f"metric K = M (MuJoCo Panda): midsection moves by up to {np.linalg.norm(QMm-QM, axis=1).max():.3f} rad (2-norm); "
          f"min margin {np.minimum(QMm-rk.QMIN, rk.QMAX-QMm).min():.4f}")
    np.save("midsection_M_m100.npy", QMm)
except Exception as e:
    print("K = M skipped:", e)

# (4) a box in which the bounding face switches: tighten joint 1 symmetrically to |q1| <= c
q1ends = []
for a in arcs:
    P, *_ = ends_and_arc(a, rk.QMIN, rk.QMAX); q1ends.append((abs(P[0][0]), abs(P[-1][0])))
q1ends = np.array(q1ends); print("|q1| at the two arc ends over s: lower end [%.3f,%.3f], upper end [%.3f,%.3f]" %
      (q1ends[:,0].min(), q1ends[:,0].max(), q1ends[:,1].min(), q1ends[:,1].max()))
c = 0.5*(q1ends[:, 0].min()+q1ends[:, 0].max())
lo2, hi2 = rk.QMIN.copy(), rk.QMAX.copy(); lo2[0], hi2[0] = -c, c
QM2, A2, B2, JB2, JF2, _ = run(lo2, hi2)
sw = np.where(np.roll(JB2, -1) != JB2)[0]
print(f"modified box |q1| <= {c:.3f}: lower-end bounding joint switches at nodes {sw.tolist()} "
      f"(joints {sorted(set((JB2+1).tolist()))}); upper-end joints {sorted(set((JF2+1).tolist()))}")
t2 = 0.5*(A2+B2); d1 = fd(t2, 1)/h; d2 = fd(t2, 2)/h**2
print("  d t_mid/ds around the switches:", [np.round(d1[max(i-3,0):i+4], 3).tolist() for i in sw])
print(f"  max |Delta^2 t_mid|/h^2: original box {np.abs(fd(tm,2)).max()/h**2:.3g}, modified box {np.abs(d2).max():.3g}")
np.savez("midsection_switch.npz", A=A2, B=B2, JB=JB2, JF=JF2, A0=A, B0=B, c=c)
