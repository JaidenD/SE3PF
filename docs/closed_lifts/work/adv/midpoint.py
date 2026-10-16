"""E5: the admissible self-motion arc [a(s), b(s)] (Euclidean arc length from the rest-posture lift, joint limits),
which joint bounds each end, and the arc-midpoint lift."""
import numpy as np
from harmonic import Fibres, load_rest, energy, report
from kinb import *

def flow_vec(F, Q, c, step):
    """move node k by arc length c[k] (vector) along its self-motion."""
    n = max(1, int(np.ceil(np.abs(c).max() / step))); hs = c / n
    for _ in range(n):
        X1 = F.kern(Q); Xm = F.kern(F.proj(Q + 0.5 * hs[:, None] * X1))
        Q = F.move(Q, hs, Xm)
    return Q

def slack(Q):
    return np.minimum(Q - QMIN, QMAX - Q)   # (N,7), negative = violated

def exit_arc(F, Q0, direction, coarse=0.02, fine=0.001, maxlen=6.0):
    """arc length to the first joint-limit violation along +-x, and the violating joint."""
    N = F.N; Q = Q0.copy(); th = np.zeros(N); done = np.zeros(N, bool)
    Qlast = Q.copy(); thlast = th.copy()
    while not done.all() and th.max() < maxlen:
        step = np.where(done, 0.0, direction * coarse)
        Qn = flow_vec(F, Q, step, coarse)
        bad = (slack(Qn).min(1) < 0) & ~done
        Qlast[~done & ~bad] = Qn[~done & ~bad]; thlast[~done & ~bad] += coarse
        done |= bad
        Q = np.where(done[:, None], Q, Qn)
    # fine pass from the last feasible point
    Q = Qlast.copy(); th = thlast.copy(); res = np.full(N, np.nan); jnt = np.zeros(N, int)
    for i in range(int(round(coarse / fine)) + 2):
        Qn = flow_vec(F, Q, np.full(N, direction * fine), fine)
        s0 = slack(Q).min(1); s1 = slack(Qn).min(1)
        hit = (s1 < 0) & np.isnan(res)
        frac = s0[hit] / (s0[hit] - s1[hit])
        res[hit] = th[hit] + frac * fine
        jnt[hit] = np.argmin(slack(Qn)[hit], 1) + 1
        Q = Qn; th = th + fine
    return res, jnt

if __name__ == '__main__':
    ss, Qrest = load_rest(500); F = Fibres(ss)
    b, jb = exit_arc(F, Qrest, +1)
    a, ja = exit_arc(F, Qrest, -1)
    a = -a
    print('upper end b(s): min %.4f max %.4f, bounding joints %s' % (np.nanmin(b), np.nanmax(b), np.unique(jb)))
    print('lower end a(s): min %.4f max %.4f, bounding joints %s' % (np.nanmin(a), np.nanmax(a), np.unique(ja)))
    print('any nan', np.isnan(a).any(), np.isnan(b).any())
    # smoothness: second differences
    d2a = np.abs(np.roll(a, -1) - 2 * a + np.roll(a, 1)) / F.h**2; d2b = np.abs(np.roll(b, -1) - 2 * b + np.roll(b, 1)) / F.h**2
    print('max |a\'\'| %.2f  max |b\'\'| %.2f (finite differences)' % (d2a.max(), d2b.max()))
    mid = 0.5 * (a + b)
    print('arc length of admissible arc: min %.4f max %.4f; midpoint offset from rest lift: min %.4f max %.4f' % ((b - a).min(), (b - a).max(), mid.min(), mid.max()))
    Qm = flow_vec(F, Qrest.copy(), mid, 0.01)
    E = energy(F, Qm)
    print('midpoint lift E=%.4f (rest %.4f)' % (E, energy(F, Qrest)))
    report('midpoint lift', F, Qm, Qrest, E, 0.0, np.array([np.nan]))
    np.savez('midpoint_500.npz', a=a, b=b, ja=ja, jb=jb, Qm=Qm)
