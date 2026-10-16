"""Multi-start check of the kinetic-energy (K = M, MuJoCo menagerie Panda, armature 0.1) harmonic lift:
E_M = sum_k d_k^T M(m_k) d_k / (2h), d_k = q_{k+1}-q_k, m_k = midpoint. Riemannian Newton in a moving slice chart,
Hessian (cyclic tridiagonal) by central differences of the chart gradient with 5-colouring."""
import sys, time
import numpy as np
sys.path.insert(0, '../num')
from pandakin import MjDyn
from harmonic import Fibres, load_rest, report
from kinb import fkb, QMIN, QMAX, feasible, L

dyn = MjDyn()
def energy(F, Q):
    d = np.roll(Q, -1, 0) - Q; m = Q + 0.5 * d
    return 0.5 * sum(dk @ dyn.M(mk) @ dk for dk, mk in zip(d, m)) / F.h

def grad_q(F, Q, eps=1e-4):
    d = np.roll(Q, -1, 0) - Q; m = Q + 0.5 * d
    Md = np.zeros_like(Q); G = np.zeros_like(Q)
    for k in range(F.N):
        dk, mk = d[k], m[k]; Mk = dyn.M(mk); Md[k] = Mk @ dk
        nd = np.linalg.norm(dk); dh = dk / nd
        Mdot = (-dyn.M(mk + 2 * eps * dh) + 8 * dyn.M(mk + eps * dh) - 8 * dyn.M(mk - eps * dh) + dyn.M(mk - 2 * eps * dh)) / (12 * eps) * nd
        G[k] = Mdot @ dk - dyn.Cv(mk, dk)
    return (np.roll(Md, 1, 0) - Md) / F.h + (G + np.roll(G, 1, 0)) / (2 * F.h)

def chart_grad(F, Q, X, th):
    Qt = F.move(Q, th, X) if np.any(th != 0) else Q
    # tangent of the slice chart at Qt: J dq = 0, X^T dq = 1
    J = fkb(Qt)[2]; A = np.concatenate([J, X[:, None, :]], 1)
    rhs = np.zeros((F.N, 7, 1)); rhs[:, 6, 0] = 1
    T = np.linalg.solve(A, rhs)[..., 0]
    return np.einsum('ki,ki->k', T, grad_q(F, Qt)), Qt

def newton(F, Q, itmax=40, tol=1e-9, ncol=5, dth=1e-5):
    E = energy(F, Q)
    for it in range(itmax):
        X = F.kern(Q)
        g, _ = chart_grad(F, Q, X, np.zeros(F.N))
        if np.abs(g).max() / F.h < 1e-7: break
        H = np.zeros((F.N, F.N)); idx = np.arange(F.N)
        for c in range(ncol):
            e = (idx % ncol == c).astype(float)
            gp, _ = chart_grad(F, Q, X, dth * e); gm, _ = chart_grad(F, Q, X, -dth * e)
            col = (gp - gm) / (2 * dth)
            for k in np.where(e > 0)[0]:
                for j in (k - 1, k, (k + 1) % F.N):
                    H[j, k] = col[j]
        H = 0.5 * (H + H.T)
        w = np.linalg.eigvalsh(H)
        Hs = H if w[0] > 1e-8 else H + (1e-6 - w[0]) * np.eye(F.N)
        dd = -np.linalg.solve(Hs, g); t = 1.0
        while True:
            Qn = F.move(Q, t * dd, X); En = energy(F, Qn)
            if En <= E + 1e-4 * t * (g @ dd) or t < 1e-6: break
            t *= 0.5
        Q, E = Qn, En
    return Q, E, np.abs(g).max() / F.h, w

if __name__ == '__main__':
    ss, Qrest = load_rest(500); F = Fibres(ss)
    print('rest lift E_M = %.6f' % energy(F, Qrest))
    for c in [0.0, -3.0, -2.0, 2.0, 3.0]:
        t0 = time.time()
        Q0 = F.flow(Qrest.copy(), c) if c != 0 else Qrest.copy()
        Q, E, el, w = newton(F, Q0)
        report('K=M harmonic, start shifted %+.1f' % c, F, Q, Qrest, E, el, w)
        print('   (%.0f s)' % (time.time() - t0)); sys.stdout.flush()
        np.save('harmM_500_start%+.0f.npy' % c, Q)
