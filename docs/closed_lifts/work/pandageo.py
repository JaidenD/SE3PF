"""Panda kinematics, the coning-circle path and SO(3) helpers (vectorized numpy).

Independent re-implementation used only to stress-test the claims of theory.md numerically.
Conventions as in the paper: twists (omega, v), body Jacobian, W = diag(wR I, wp I).
"""
import numpy as np

# Modified (Craig) DH: a_{i-1}, d_i, alpha_{i-1}; flange offset 0.107 m
A_DH = np.array([0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088])
D_DH = np.array([0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0])
AL_DH = np.array([0.0, -np.pi / 2, np.pi / 2, np.pi / 2, -np.pi / 2, np.pi / 2, np.pi / 2])
D_FLANGE = 0.107
QMIN = np.array([-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973])
QMAX = np.array([2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973])
QR = np.array([0.0, -np.pi / 4, 0.0, -3 * np.pi / 4, 0.0, np.pi / 2, np.pi / 4])


def hat(w):
    w = np.asarray(w, dtype=float)
    H = np.zeros(w.shape[:-1] + (3, 3))
    H[..., 0, 1] = -w[..., 2]; H[..., 0, 2] = w[..., 1]
    H[..., 1, 0] = w[..., 2];  H[..., 1, 2] = -w[..., 0]
    H[..., 2, 0] = -w[..., 1]; H[..., 2, 1] = w[..., 0]
    return H


def expSO3(w):
    w = np.asarray(w, dtype=float)
    th = np.linalg.norm(w, axis=-1)[..., None, None]
    W = hat(w)
    small = th < 1e-8
    ths = np.where(small, 1.0, th)
    a = np.where(small, 1.0 - th**2 / 6, np.sin(ths) / ths)
    b = np.where(small, 0.5 - th**2 / 24, (1 - np.cos(ths)) / ths**2)
    return np.eye(3) + a * W + b * (W @ W)


def logSO3(R):
    tr = np.trace(R, axis1=-2, axis2=-1)
    c = np.clip((tr - 1) / 2, -1.0, 1.0)
    th = np.arccos(c)
    v = np.stack([R[..., 2, 1] - R[..., 1, 2], R[..., 0, 2] - R[..., 2, 0], R[..., 1, 0] - R[..., 0, 1]], axis=-1)
    small = th < 1e-8
    f = np.where(small, 0.5 + th**2 / 12, th / (2 * np.sin(np.where(small, 1.0, th))))
    return f[..., None] * v


def fk(q):
    """q (...,7) -> R (...,3,3), p (...,3), body Jacobian Jb (...,6,7), rows (omega, v)."""
    q = np.asarray(q, dtype=float)
    shp = q.shape[:-1]
    R = np.broadcast_to(np.eye(3), shp + (3, 3)).copy()
    p = np.zeros(shp + (3,))
    Z = np.zeros(shp + (7, 3)); P = np.zeros(shp + (7, 3))
    for i in range(7):
        c, s = np.cos(q[..., i]), np.sin(q[..., i])
        ca, sa = np.cos(AL_DH[i]), np.sin(AL_DH[i])
        ca = 0.0 if abs(ca) < 1e-15 else ca
        sa = 0.0 if abs(sa) < 1e-15 else sa
        Ri = np.zeros(shp + (3, 3))
        Ri[..., 0, 0] = c; Ri[..., 0, 1] = -s
        Ri[..., 1, 0] = s * ca; Ri[..., 1, 1] = c * ca; Ri[..., 1, 2] = -sa
        Ri[..., 2, 0] = s * sa; Ri[..., 2, 1] = c * sa; Ri[..., 2, 2] = ca
        ti = np.array([A_DH[i], -D_DH[i] * sa, D_DH[i] * ca])
        p = p + np.einsum('...ab,b->...a', R, ti)
        R = R @ Ri
        Z[..., i, :] = R[..., :, 2]
        P[..., i, :] = p
    p = p + D_FLANGE * R[..., :, 2]
    Rt = np.swapaxes(R, -1, -2)
    w = Z
    v = np.cross(Z, p[..., None, :] - P)
    wb = np.einsum('...ab,...ib->...ai', Rt, w)
    vb = np.einsum('...ab,...ib->...ai', Rt, v)
    J = np.concatenate([wb, vb], axis=-2)
    return R, p, J


def kernel_vector(J):
    """Unit generalized cross product of the six rows of J (...,6,7): spans ker J, continuous orientation."""
    n = np.zeros(J.shape[:-2] + (7,))
    for j in range(7):
        cols = [k for k in range(7) if k != j]
        n[..., j] = (-1) ** j * np.linalg.det(J[..., :, cols])
    return n / np.linalg.norm(n, axis=-1, keepdims=True)


def pinv(J):
    JJt = J @ np.swapaxes(J, -1, -2)
    return np.swapaxes(J, -1, -2) @ np.linalg.inv(JJt)


class ConingCircle:
    """R = R0 Rz(th) A Rz(-th), p = c + r (cos th, sin th, 0); arc length in W = diag(wR I, wp I)."""

    def __init__(self, wR=0.01, wp=1.0, center=(0.40, 0.0, 0.45), radius=0.12, cone=0.25):
        self.wR, self.wp = wR, wp
        self.Wd = np.array([wR] * 3 + [wp] * 3)
        self.R0 = fk(QR)[0]
        self.c = np.array(center); self.r = radius
        self.A = expSO3(np.array([cone, 0.0, 0.0]))
        e3 = np.array([0.0, 0.0, 1.0])
        self.w0 = self.A.T @ e3 - e3
        self.sigma = np.sqrt(4 * wR * np.sin(cone / 2) ** 2 + wp * radius**2)
        self.L = 2 * np.pi * self.sigma

    @staticmethod
    def Rz(th):
        th = np.asarray(th, dtype=float)
        z = np.zeros(th.shape + (3,)); z[..., 2] = th
        return expSO3(z)

    def pose(self, s):
        th = np.asarray(s, dtype=float) / self.sigma
        Rz = self.Rz(th)
        R = self.R0 @ Rz @ self.A @ np.swapaxes(Rz, -1, -2)
        p = self.c + self.r * np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], axis=-1)
        return R, p

    def velocity(self, s):
        th = np.asarray(s, dtype=float) / self.sigma
        R, _ = self.pose(s)
        pd = self.r * np.stack([-np.sin(th), np.cos(th), np.zeros_like(th)], axis=-1)
        w = np.einsum('...ab,b->...a', self.Rz(th), self.w0)
        v = np.einsum('...ba,...b->...a', R, pd)
        return np.concatenate([w, v], axis=-1) / self.sigma

    def acceleration(self, s):
        th = np.asarray(s, dtype=float) / self.sigma
        R, _ = self.pose(s)
        Rz = self.Rz(th)
        w = np.einsum('...ab,b->...a', Rz, self.w0)
        dw = np.einsum('...ab,b->...a', Rz, np.cross(np.array([0, 0, 1.0]), self.w0))
        pd = self.r * np.stack([-np.sin(th), np.cos(th), np.zeros_like(th)], axis=-1)
        pdd = -self.r * np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], axis=-1)
        v = np.einsum('...ba,...b->...a', R, pd)
        dv = np.einsum('...ba,...b->...a', R, pdd) - np.cross(w, v)
        return np.concatenate([dw, dv], axis=-1) / self.sigma**2


def pose_error(q, path, s):
    """Body-frame error e = Log_I(Pi(q)^{-1} gamma(s)) of the product metric, and J_b(q)."""
    R, p, J = fk(q)
    Rg, pg = path.pose(s)
    Rt = np.swapaxes(R, -1, -2)
    e = np.concatenate([logSO3(Rt @ Rg), np.einsum('...ab,...b->...a', Rt, pg - p)], axis=-1)
    return e, J
