"""E3/E4: closed lifts minimizing E_k(q) = sum_k |q_{k+1}-q_k|^2/(2h) + (kp/2) h sum_k |q_k - q_r|^2 over closed
sections (q_k on the fibre over gamma(s_k)), Euclidean joint metric, by Riemannian Newton on the product of fibres
(moving slice x_k^T q = const, x_k the unit kernel vector at the current q_k)."""
import numpy as np
from kinb import *
from kin import QR

def load_rest(N=500):
    rest = np.genfromtxt('../num/data/rest_%d.csv' % N, delimiter=',', names=True)
    return rest['s'], np.stack([rest['q%d' % i] for i in range(7)], 1)

class Fibres:
    def __init__(self, ss):
        self.ss = ss; self.N = len(ss); self.h = L / self.N
        self.Rg, self.pg, self.Z = path_b(ss)
    def proj(self, Q, C=None, cval=None):
        return project_b(Q, self.Rg, self.pg, C, cval)
    def kern(self, Q):
        return kernel_b(fkb(Q)[2])
    def move(self, Q, dth, X=None):
        """move each q_k along its fibre by slice amount dth_k (slice normal X_k)."""
        if X is None: X = self.kern(Q)
        Qp = Q + dth[:, None] * X
        return self.proj(Qp, X, np.einsum('mi,mi->m', X, Qp))
    def flow(self, Q, c, step=0.02):
        """move each q_k by Euclidean arc length c along its self-motion (midpoint rule + projection)."""
        n = max(1, int(np.ceil(abs(c) / step))); hs = c / n
        for _ in range(n):
            X1 = self.kern(Q); Xm = self.kern(self.proj(Q + 0.5 * hs * X1))
            Q = self.move(Q, np.full(self.N, hs), Xm)
        return Q

def energy(F, Q, kp=0.0, qr=QR):
    d = np.roll(Q, -1, 0) - Q
    return 0.5 * np.sum(d * d) / F.h + 0.5 * kp * F.h * np.sum((Q - qr) ** 2)

def grad_hess(F, Q, kp=0.0, qr=QR, eps=1e-4):
    X = F.kern(Q)
    lap = 2 * Q - np.roll(Q, -1, 0) - np.roll(Q, 1, 0)          # -(second difference)
    G = lap / F.h + kp * F.h * (Q - qr)                           # dE/dq_k
    g = np.einsum('mi,mi->m', X, G)
    # curvature of the fibre in the slice chart
    Qp = F.move(Q, np.full(F.N, eps), X); Qm = F.move(Q, np.full(F.N, -eps), X)
    kap = (Qp - 2 * Q + Qm) / eps**2
    diag = 2 / F.h + kp * F.h + np.einsum('mi,mi->m', kap, G)
    off = -np.einsum('mi,mi->m', X, np.roll(X, -1, 0)) / F.h       # H[k, k+1]
    H = np.diag(diag) + np.diag(off[:-1], 1) + np.diag(off[:-1], -1)
    H[-1, 0] = H[0, -1] = off[-1]
    return g, H, X

def solve(F, Q, kp=0.0, qr=QR, tol=1e-10, itmax=200, verbose=False):
    E = energy(F, Q, kp, qr)
    for it in range(itmax):
        g, H, X = grad_hess(F, Q, kp, qr)
        if np.abs(g).max() * F.h < tol: break   # g_k ~ h * (EL residual)
        w = np.linalg.eigvalsh(H)
        if w[0] <= 1e-8:
            H = H + (1e-6 - w[0]) * np.eye(F.N) * 1.0
        d = -np.linalg.solve(H, g)
        t = 1.0
        while True:
            Qn = F.move(Q, t * d, X); En = energy(F, Qn, kp, qr)
            if En <= E + 1e-4 * t * (g @ d) or t < 1e-6: break
            t *= 0.5
        Q, E = Qn, En
        if verbose: print(it, E, np.abs(g).max(), t)
    g, H, X = grad_hess(F, Q, kp, qr)
    w = np.linalg.eigvalsh(H)
    return Q, E, np.abs(g).max() / F.h, w

def report(name, F, Q, Qref, E, el, w):
    dist = np.linalg.norm(Q - Qref, axis=1)
    mar = np.minimum(Q - QMIN, QMAX - Q); m = mar.min(axis=1)
    sm = np.linalg.svd(fkb(Q)[2], compute_uv=False)[:, -1]
    print(f"{name}: E={E:.6f}  EL-resid={el:.1e}  Hess eig min/max {w[0]:.3e}/{w[-1]:.3e}  "
          f"max|q-q*_rest|={dist.max():.4f} (mean {dist.mean():.4f})  min margin {m.min():.4f} (joint {np.argmin(mar.min(0))+1})  "
          f"min sigma_min {sm.min():.4f}  feasible {bool(feasible(Q).all())}")
