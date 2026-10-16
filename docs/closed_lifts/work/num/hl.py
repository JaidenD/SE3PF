"""Discrete harmonic lifts of the Panda path: sections of the self-motion bundle over the loop R/LZ.

A closed lift is a periodic grid function q_k, k = 0..N-1, s_k = k L / N, with FK(q_k) = gamma(s_k). Chart: u_k is
the controller's beta at node k, x_k^T (q_k - q*_rest(s_k)) = u_k, x_k the unit kernel vector of J_b at the
rest-posture lift node. Discrete energy (periodic, q_N = q_0, d_k = q_{k+1} - q_k, h = L / N):
    K = I : E = sum_k |d_k|^2 / (2h)
    K = M : E = sum_k d_k^T M(m_k) d_k / (2h),  m_k = (q_k + q_{k+1}) / 2   (midpoint rule, MuJoCo M with armature)
plus an optional penalty (kpen / 2) h sum_k |q_k - q_r|^2 (the family E_k).
Gradient in u is exact up to roundoff: dE/du_k = t_k . dE/dq_k, t_k = dq_k/du_k from [J_b(q_k); x_k^T] t = (0, 1);
for K = M, grad_q (v^T M(q) v) = 2 (Mdot_v v - C(q, v) v) with C v from MuJoCo's RNE (gravity off) and
Mdot_v = directional derivative of M along v (fourth-order central difference, step 1e-4 rad).
Hessian in u: the energy couples u_k only to u_{k-1}, u_{k+1}, so the Hessian is periodic tridiagonal and is obtained
from 10 gradient evaluations (5 colours, central differences, du = 1e-5).
"""
import numpy as np
from pandakin import (fk, kernel_vector, constrained_ik, path_pose, path_vel, load_rest, joint_margin, MjDyn,
                      L, QR, pose_error)


def spectral_deriv(y, Lp, order=1):
    """Derivative of periodic samples y[k] (axis 0), k = 0..N-1, period Lp, by FFT."""
    N = y.shape[0]
    w = 2 * np.pi * np.fft.fftfreq(N, d=Lp / N)
    Y = np.fft.fft(y, axis=0)
    f = (1j * w) ** order
    if order % 2 == 1 and N % 2 == 0:
        f[N // 2] = 0.0
    return np.real(np.fft.ifft(f[:, None] * Y if y.ndim > 1 else f * Y, axis=0))


def trig_interp(y, Lp, s):
    """Trigonometric interpolation of periodic samples y (N x d) at points s."""
    N = y.shape[0]
    Y = np.fft.fft(y, axis=0) / N
    k = np.fft.fftfreq(N, d=1.0 / N)
    if N % 2 == 0:
        Y[N // 2] *= 0.5
        Y = np.concatenate([Y, Y[N // 2:N // 2 + 1]], 0)
        k = np.concatenate([k, [N // 2]])
    E = np.exp(2j * np.pi * np.outer(np.asarray(s), k) / Lp)
    return np.real(E @ Y)


class Problem:
    def __init__(self, N, metric="I", kpen=0.0, dyn=None, qref=QR, lim_kappa=0.0, lim_delta=0.10):
        self.N = N
        self.h = L / N
        self.s, self.qrest, self.x = load_rest(N)
        self.Rt, self.pt = path_pose(self.s)
        self.b0 = np.einsum("ki,ki->k", self.x, self.qrest)
        self.metric = metric
        self.kpen = kpen
        self.qref = np.asarray(qref)
        # optional quadratic penalty (lim_kappa / 2) h sum max(0, violation of [qmin + delta, qmax - delta])^2
        self.lim_kappa, self.lim_delta = lim_kappa, lim_delta
        self.dyn = dyn if dyn is not None else (MjDyn() if metric == "M" else None)
        self.qwarm = self.qrest.copy()
        self.ngrad = 0

    # -------------------------------------------------------------------------------------------- chart
    def chart(self, u, guess=None):
        q, ok = constrained_ik(self.Rt, self.pt, self.x, self.b0 + u, self.qwarm if guess is None else guess,
                               tol=1e-13, maxit=80, extra=1)
        if not ok:
            raise RuntimeError("chart: constrained IK failed")
        return q

    def u_of_q(self, q):
        return np.einsum("ki,ki->k", self.x, q) - self.b0

    def tangents(self, q):
        _, _, J = fk(q)
        A = np.concatenate([J, self.x[:, None, :]], 1)
        rhs = np.zeros((self.N, 7, 1)); rhs[:, 6, 0] = 1.0
        return np.linalg.solve(A, rhs)[..., 0], J

    # -------------------------------------------------------------------------------------------- energy
    def _M_mid(self, q):
        m = 0.5 * (q + np.roll(q, -1, 0))
        return m, self.dyn.Mbatch(m)

    def energy(self, q):
        d = np.roll(q, -1, 0) - q
        if self.metric == "I":
            E = 0.5 * np.sum(d * d) / self.h
        else:
            _, M = self._M_mid(q)
            E = 0.5 * np.einsum("ki,kij,kj->", d, M, d) / self.h
        if self.kpen:
            E += 0.5 * self.kpen * self.h * np.sum((q - self.qref) ** 2)
        if self.lim_kappa:
            E += 0.5 * self.lim_kappa * self.h * np.sum(self._viol(q) ** 2)
        return E

    def _viol(self, q):
        from pandakin import QMIN, QMAX
        return np.maximum(0.0, q - (QMAX - self.lim_delta)) - np.maximum(0.0, (QMIN + self.lim_delta) - q)

    def grad_q(self, q):
        self.ngrad += 1
        d = np.roll(q, -1, 0) - q
        if self.metric == "I":
            g = (q - np.roll(q, 1, 0)) / self.h - d / self.h
        else:
            m, M = self._M_mid(q)
            Md = np.einsum("kij,kj->ki", M, d)
            G = np.zeros_like(q)
            eps = 1e-4  # fourth-order central difference of M along d_k: truncation ~eps^4, roundoff ~1e-12
            for k in range(self.N):
                dk = d[k]
                nd = np.linalg.norm(dk)
                if nd == 0:
                    continue
                dh = dk / nd
                Mf = self.dyn.M
                Mdot = (-Mf(m[k] + 2 * eps * dh) + 8 * Mf(m[k] + eps * dh) - 8 * Mf(m[k] - eps * dh)
                        + Mf(m[k] - 2 * eps * dh)) / (12 * eps) * nd
                G[k] = Mdot @ dk - self.dyn.Cv(m[k], dk)
            g = (np.roll(Md, 1, 0) - Md) / self.h + (G + np.roll(G, 1, 0)) / (2 * self.h)
        if self.kpen:
            g = g + self.kpen * self.h * (q - self.qref)
        if self.lim_kappa:
            g = g + self.lim_kappa * self.h * self._viol(q)
        return g

    def grad_u(self, u, guess=None):
        q = self.chart(u, guess)
        t, _ = self.tangents(q)
        return np.einsum("ki,ki->k", t, self.grad_q(q)), q

    def hessian_u(self, u, q, du=1e-5, ncol=5):
        N = self.N
        assert N % ncol == 0
        H = np.zeros((N, N))
        idx = np.arange(N)
        for c in range(ncol):
            e = (idx % ncol == c).astype(float)
            gp, _ = self.grad_u(u + du * e, guess=q)
            gm, _ = self.grad_u(u - du * e, guess=q)
            col = (gp - gm) / (2 * du)
            for j in idx[idx % ncol == c]:
                for i in ((j - 1) % N, j, (j + 1) % N):
                    H[i, j] = col[i]
        return 0.5 * (H + H.T)

    # -------------------------------------------------------------------------------------------- Newton
    def solve(self, u0, tol=1e-9, maxit=40, verbose=False):
        """Newton with backtracking on E; stops when max |dE/du_k| / h < tol."""
        u = np.array(u0, dtype=float)
        q = self.chart(u)
        self.qwarm = q.copy()
        E = self.energy(q)
        hist = []
        for it in range(maxit):
            g, q = self.grad_u(u)
            gn = np.abs(g).max() / self.h
            hist.append((it, E, gn))
            if verbose:
                print("   it %2d  E = %.15f  max|g|/h = %.3e" % (it, E, gn), flush=True)
            if gn < tol:
                break
            H = self.hessian_u(u, q)
            w = np.linalg.eigvalsh(H)
            if w[0] <= 0:  # not convex here: shift to a descent direction
                H = H + (abs(w[0]) + 1e-6 * abs(w).max()) * np.eye(self.N)
            step = -np.linalg.solve(H, g)
            a = 1.0
            while True:
                un = u + a * step
                try:
                    qn = self.chart(un, guess=q)
                    En = self.energy(qn)
                except RuntimeError:
                    En = np.inf
                if En <= E + 1e-4 * a * g @ step or a < 1e-8:
                    break
                a *= 0.5
            if En > E + 1e-12 * abs(E) and gn > 1e3 * tol:
                print("   line search failed", flush=True)
            u, q, E = un, qn, En
            self.qwarm = q.copy()
        return u, q, E, hist

    # -------------------------------------------------------------------------------------------- diagnostics
    def el_residual(self, q):
        """x(q)^T (K q'' + C(q, q') q') along the lift, derivatives by FFT of the node values (spectral)."""
        dq = spectral_deriv(q, L, 1)
        ddq = spectral_deriv(q, L, 2)
        _, _, J = fk(q)
        x = kernel_vector(J)
        if self.metric == "I":
            F = ddq
        else:
            F = np.array([self.dyn.M(q[k]) @ ddq[k] + self.dyn.Cv(q[k], dq[k]) for k in range(self.N)])
        r = np.einsum("ki,ki->k", x, F)
        return r, np.linalg.norm(F, axis=1), dq

    def el_residual_fd4(self, q):
        """The same with fourth-order central differences (q' and q'')."""
        h = self.h
        r = lambda j: np.roll(q, -j, 0)
        dq = (r(-2) - 8 * r(-1) + 8 * r(1) - r(2)) / (12 * h)
        ddq = (-r(-2) + 16 * r(-1) - 30 * q + 16 * r(1) - r(2)) / (12 * h * h)
        _, _, J = fk(q)
        x = kernel_vector(J)
        if self.metric == "I":
            F = ddq
        else:
            F = np.array([self.dyn.M(q[k]) @ ddq[k] + self.dyn.Cv(q[k], dq[k]) for k in range(self.N)])
        return np.einsum("ki,ki->k", x, F)

    def summary(self, q):
        _, _, J = fk(q)
        sig = np.linalg.svd(J, compute_uv=False)[:, -1]
        mar, jj = joint_margin(q)
        k = int(np.argmin(mar))
        dq = spectral_deriv(q, L, 1)
        return dict(minMargin=float(mar[k]), marginJoint=int(jj[k]) + 1, marginS=float(self.s[k]),
                    minSigma=float(sig.min()), sigmaS=float(self.s[np.argmin(sig)]),
                    maxDq=float(np.linalg.norm(dq, axis=1).max()),
                    maxDistRest=float(np.linalg.norm(q - self.qrest, axis=1).max()),
                    maxDistRestInf=float(np.abs(q - self.qrest).max()))

    def chart_quality(self, q):
        """cos of the angle between the chart vector x_k and the fibre tangent x(q_k): the chart u is valid while > 0."""
        _, _, J = fk(q)
        return np.abs(np.einsum("ki,ki->k", self.x, kernel_vector(J))).min()
