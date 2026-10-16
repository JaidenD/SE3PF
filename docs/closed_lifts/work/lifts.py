"""Lifts of the coning-circle path to the Panda's joint space; harmonic sections of the self-motion bundle.

Stress tests for theory.md. Everything uses the Euclidean joint metric K = I unless stated.
"""
import numpy as np
from pandageo import fk, kernel_vector, pinv, pose_error, ConingCircle, QR, QMIN, QMAX, expSO3, logSO3


# ----------------------------------------------------------------------------------------------------------
# Fibre retraction: q on the fibre over gamma(s) with x^T (q - q0) = u (Newton on the 7x7 system [J_b; x^T])
def retract(path, s, q0, x, u, tol=1e-13, maxit=30):
    """Vectorized over nodes: s (m,), q0 (m,7), x (m,7), u (m,). Returns q (m,7)."""
    q = q0 + u[:, None] * x
    for _ in range(maxit):
        e, J = pose_error(q, path, s)
        r = u - np.einsum('ma,ma->m', x, q - q0)
        A = np.concatenate([J, x[:, None, :]], axis=1)
        b = np.concatenate([e, r[:, None]], axis=1)
        dq = np.linalg.solve(A, b[..., None])[..., 0]
        q = q + dq
        if np.abs(dq).max() < tol:
            break
    return q


def rest_lift(path, m, qr=QR):
    """Pointwise rest-posture lift by projected Newton + continuation (the paper's Implementation section)."""
    s = np.arange(m + 1) * path.L / m
    Q = np.zeros((m + 1, 7))
    q = qr.copy()
    for j in range(m + 1):
        for it in range(200):
            e, J = pose_error(q, path, s[j])
            Jp = pinv(J)
            dq = Jp @ e + (np.eye(7) - Jp @ J) @ (qr - q)
            nrm = np.linalg.norm(dq)
            if nrm > 0.2:
                dq *= 0.2 / nrm
            q = q + dq
            if nrm < 1e-13:
                break
        Q[j] = q
    return s[:m], Q[:m], np.linalg.norm(Q[m] - Q[0])


# ----------------------------------------------------------------------------------------------------------
# Horizontal transport for the K-orthogonal (pseudoinverse) connection, with a weak pose feedback to stay on
# the fibre: q' = J^+ (zeta_gamma(s) + kfb e)
def transport(path, q0, s0, s1, n=4000, kfb=50.0):
    h = (s1 - s0) / n
    q = q0.copy()

    def f(q, s):
        e, J = pose_error(q, path, s)
        return pinv(J) @ (path.velocity(s) + kfb * e)

    s = s0
    for i in range(n):
        k1 = f(q, s); k2 = f(q + 0.5 * h * k1, s + 0.5 * h)
        k3 = f(q + 0.5 * h * k2, s + 0.5 * h); k4 = f(q + h * k3, s + h)
        q = q + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        s += h
    return q


def self_motion_curve(path, s, q0, ds=1e-3, nmax=20000, margin=0.0):
    """Trace the self-motion curve over gamma(s) through q0 in both directions (unit kernel field, RK4 with
    pose correction) until a joint limit (minus margin) or nmax steps. Returns arclength t and points."""
    out = {}
    for sign in (+1, -1):
        q = q0.copy(); ts = [0.0]; qs = [q.copy()]
        xref = kernel_vector(fk(q0)[2])

        def f(q):
            e, J = pose_error(q, path, s)
            x = kernel_vector(J)
            return sign * np.sign(x @ xref) * x + pinv(J) @ (20.0 * e)

        for i in range(nmax):
            k1 = f(q); k2 = f(q + 0.5 * ds * k1); k3 = f(q + 0.5 * ds * k2); k4 = f(q + ds * k3)
            qn = q + (ds / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            if np.any(qn < QMIN + margin) or np.any(qn > QMAX - margin):
                break
            if i > 10 and np.linalg.norm(qn - q0) < 0.6 * ds:  # closed self-motion circle
                out['closed'] = True
                break
            q = qn; ts.append(sign * ds * (i + 1)); qs.append(q.copy())
        out[sign] = (np.array(ts), np.array(qs))
    t = np.concatenate([out[-1][0][::-1], out[1][0][1:]])
    Q = np.concatenate([out[-1][1][::-1], out[1][1][1:]])
    return t, Q, out.get('closed', False)


# ----------------------------------------------------------------------------------------------------------
# Discrete harmonic sections. Unknowns: one fibre coordinate per node, in a chart re-centred at the current
# iterate (u_j = x_j^T (q_j^new - q_j)). Energies are sums of pair terms H_j(q_j, q_{j+1}).
def fibre_jets(path, s, q):
    """x = unit kernel vector, kappa = second derivative of the chart curve = curvature vector of the
    self-motion curve (Euclidean), at the nodes."""
    _, _, J = fk(q)
    x = kernel_vector(J)
    eps = 1e-6
    _, _, Jp = fk(q + eps * x); _, _, Jm = fk(q - eps * x)
    DJx = np.einsum('mab,mb->ma', (Jp - Jm) / (2 * eps), x)
    kappa = -np.einsum('mab,mb->ma', pinv(J), DJx)
    return x, kappa, J


def energy_full(q, h, k=0.0, qr=QR):
    d = np.roll(q, -1, axis=0) - q
    return 0.5 / h * np.sum(d * d) + 0.5 * k * h * np.sum((q - qr) ** 2)


def energy_vertical(q, h):
    qm = 0.5 * (q + np.roll(q, -1, axis=0))
    xm = kernel_vector(fk(qm)[2])
    d = np.roll(q, -1, axis=0) - q
    return 0.5 / h * np.sum(np.einsum('ma,ma->m', xm, d) ** 2)


def grad_hess_full(q, x, kap, h, k=0.0, qr=QR):
    lap = 2 * q - np.roll(q, 1, axis=0) - np.roll(q, -1, axis=0)
    g = np.einsum('ma,ma->m', lap, x) / h + k * h * np.einsum('ma,ma->m', q - qr, x)
    a = (2 * np.einsum('ma,ma->m', x, x) + np.einsum('ma,ma->m', lap, kap)) / h \
        + k * h * (np.einsum('ma,ma->m', x, x) + np.einsum('ma,ma->m', q - qr, kap))
    b = -np.einsum('ma,ma->m', x, np.roll(x, -1, axis=0)) / h  # H[j, j+1]
    return g, a, b


def grad_hess_pairs(q, x, kap, h, pair_energy, delta=1e-4):
    """Gradient and cyclic-tridiagonal Hessian of sum_j H(Q_j(u_j), Q_{j+1}(u_{j+1})) at u = 0, by central
    differences on the 2-jets Q_j(u) = q_j + u x_j + u^2/2 kappa_j (exact to second order)."""
    m = q.shape[0]
    qn, xn, kn = np.roll(q, -1, 0), np.roll(x, -1, 0), np.roll(kap, -1, 0)
    def Hval(a, b):
        A = q + a * x + 0.5 * a * a * kap
        B = qn + b * xn + 0.5 * b * b * kn
        return pair_energy(A, B)
    dl = delta
    H = {}
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            H[(a, b)] = Hval(a * dl, b * dl)
    d1 = (H[(1, 0)] - H[(-1, 0)]) / (2 * dl)
    d2 = (H[(0, 1)] - H[(0, -1)]) / (2 * dl)
    d11 = (H[(1, 0)] - 2 * H[(0, 0)] + H[(-1, 0)]) / dl**2
    d22 = (H[(0, 1)] - 2 * H[(0, 0)] + H[(0, -1)]) / dl**2
    d12 = (H[(1, 1)] - H[(1, -1)] - H[(-1, 1)] + H[(-1, -1)]) / (4 * dl**2)
    g = d1 + np.roll(d2, 1)
    a = d11 + np.roll(d22, 1)
    return g, a, d12


def cyclic_tridiag_matrix(a, b):
    m = a.shape[0]
    H = np.diag(a)
    idx = np.arange(m)
    H[idx, (idx + 1) % m] += b
    H[(idx + 1) % m, idx] += b
    return H


def harmonic_section(path, s, q, h, kind='full', k=0.0, tol=1e-10, maxit=60, verbose=False):
    """Damped Newton on the discrete Euler-Lagrange equations (cyclic tridiagonal Hessian)."""
    def E(q):
        return energy_full(q, h, k) if kind == 'full' else energy_vertical(q, h)

    if kind == 'vertical':
        def pair(A, B):
            xm = kernel_vector(fk(0.5 * (A + B))[2])
            return 0.5 / h * np.einsum('ma,ma->m', xm, B - A) ** 2
    hist = []
    for it in range(maxit):
        x, kap, J = fibre_jets(path, s, q)
        if kind == 'full':
            g, a, b = grad_hess_full(q, x, kap, h, k)
        else:
            g, a, b = grad_hess_pairs(q, x, kap, h, pair)
        H = cyclic_tridiag_matrix(a, b)
        gn = np.abs(g).max()
        hist.append((it, E(q), gn))
        if verbose:
            print(f'   it {it:2d}  E = {E(q):.12f}  |grad|_inf = {gn:.3e}')
        if gn < tol:
            break
        lam = np.linalg.eigvalsh(H)
        shift = 0.0 if lam[0] > 0 else -lam[0] + 1e-3 * np.abs(lam).max()
        du = -np.linalg.solve(H + shift * np.eye(len(a)), g)
        t = 1.0; E0 = E(q)
        while True:
            qn = retract(path, s, q, x, t * du)
            if E(qn) <= E0 + 1e-4 * t * (g @ du) or t < 1e-6:
                break
            t *= 0.5
        q = qn
    x, kap, J = fibre_jets(path, s, q)
    if kind == 'full':
        g, a, b = grad_hess_full(q, x, kap, h, k)
    else:
        g, a, b = grad_hess_pairs(q, x, kap, h, pair)
    return q, (a, b), hist


def monodromy_trace(a, b):
    """Discrete Floquet test (discrete Hill operator): transfer matrices of the periodic Jacobi equation
    b_j phi_{j+1} + a_j phi_j + b_{j-1} phi_{j-1} = 0. Returns trace of the monodromy (det = 1)."""
    m = a.shape[0]
    M = np.eye(2)
    for j in range(m):
        T = np.array([[-a[j] / b[j], -b[j - 1] / b[j]], [1.0, 0.0]])
        M = T @ M
        M /= 1.0  # no rescaling needed for moderate m
    return np.trace(M), np.linalg.det(M)
