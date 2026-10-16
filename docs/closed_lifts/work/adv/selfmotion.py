"""Trace the full self-motion circle (joints on T^7, no limits) through a configuration on the fibre over gamma(s)."""
import numpy as np
from kin import *

def wrap(d):
    return (d + np.pi) % (2 * np.pi) - np.pi

def trace(q0, s, h=0.005, maxlen=60.0):
    """Returns arrays theta (Euclidean arc length from q0 along +x), Q (configs), smin (sigma_min J_b).
    Stops when the curve returns to q0 on T^7."""
    q = q0.copy(); th = 0.0
    TH, QQ, SM = [0.0], [q.copy()], []
    _, _, J = fk(q, True); SM.append(np.linalg.svd(J, compute_uv=False)[-1])
    while th < maxlen:
        _, _, J = fk(q, True); x1, _ = kernel(J)
        qm = q + 0.5 * h * x1
        _, _, Jm = fk(qm, True); xm, _ = kernel(Jm)
        qp = q + h * xm
        qn = project(qp, s, c=xm, cval=xm @ qp)
        th += np.linalg.norm(qn - q)  # chord ~ arc
        q = qn
        _, _, J = fk(q, True)
        TH.append(th); QQ.append(q.copy()); SM.append(np.linalg.svd(J, compute_uv=False)[-1])
        if th > 10 * h and np.linalg.norm(wrap(q - q0)) < 0.75 * h:
            break
    return np.array(TH), np.array(QQ), np.array(SM)

if __name__ == '__main__':
    rest = np.genfromtxt('../num/data/rest_500.csv', delimiter=',', names=True)
    qs = np.stack([rest['q%d' % i] for i in range(7)], 1); ss = rest['s']
    out = []
    for k in [0, 50, 100, 150, 200, 238, 250, 300, 350, 400, 450]:
        TH, QQ, SM = trace(qs[k], ss[k])
        ell = TH[-1] + np.linalg.norm(wrap(QQ[-1] - qs[k]))
        # admissible arc containing theta=0: inside limits (interpreting angles in (-pi, pi] then limits)
        Qw = wrap(QQ)
        ok = np.all((Qw > QMIN) & (Qw < QMAX), axis=1)
        # walk forward and backward from index 0
        n = len(TH)
        i = 0
        while i + 1 < n and ok[i + 1]: i += 1
        j = 0
        while ok[(j - 1) % n] and (j - 1) % n != i: j -= 1
        b = TH[i]; a = TH[j % n] - ell if j < 0 else 0.0
        def which(qq):
            d = np.minimum(qq - QMIN, QMAX - qq); return int(np.argmin(d)) + 1
        nb = which(Qw[min(i + 1, n - 1)]); na = which(Qw[(j - 1) % n])
        nseg = np.sum(ok[1:] != ok[:-1])
        print(f"k={k:3d} s={ss[k]:.4f}: circle length {ell:.4f} rad, min sigma_min on circle {SM.min():.4f}, "
              f"admissible arc around q*: [{a:+.4f}, {b:+.4f}] (length {b-a:.4f}), ends set by joints {na}/{nb}, "
              f"#limit crossings on circle {nseg}, margin(q*) {margin(qs[k]):.4f}")
        np.savez(f'sm_{k}.npz', TH=TH, QQ=QQ, SM=SM, ell=ell, a=a, b=b)
