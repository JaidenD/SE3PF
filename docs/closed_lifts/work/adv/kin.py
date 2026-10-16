"""Independent Panda kinematics + coning-circle path for the adversarial checks (alternatives.md).

Conventions as in the paper: twists (omega, v), body Jacobian J_b (6x7), product metric W = diag(wR I, wp I).
Validated against the project library's dump (../num/data/fk_check.csv, rest_*.csv) in check_kin.py.
"""
import numpy as np

# modified DH (Craig): a_{i-1}, d_i, alpha_{i-1}
A_ = np.array([0, 0, 0, 0.0825, -0.0825, 0, 0.088])
D_ = np.array([0.333, 0, 0.316, 0, 0.384, 0, 0])
AL = np.array([0, -np.pi / 2, np.pi / 2, np.pi / 2, -np.pi / 2, np.pi / 2, np.pi / 2])
DF = 0.107
QMIN = np.array([-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973])
QMAX = np.array([2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973])
QR = np.array([0, -np.pi / 4, 0, -3 * np.pi / 4, 0, np.pi / 2, np.pi / 4])
WR, WP = 0.01, 1.0
Wd = np.array([WR] * 3 + [WP] * 3)


def hat(w):
    return np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])


def expSO3(w):
    th = np.linalg.norm(w)
    K = hat(w)
    if th < 1e-12:
        return np.eye(3) + K
    return np.eye(3) + np.sin(th) / th * K + (1 - np.cos(th)) / th**2 * K @ K


def logSO3(R):
    c = np.clip((np.trace(R) - 1) / 2, -1, 1)
    th = np.arccos(c)
    if th < 1e-8:
        return 0.5 * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    return th / (2 * np.sin(th)) * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])


def _T(a, d, al, th):
    ca, sa, ct, st = np.cos(al), np.sin(al), np.cos(th), np.sin(th)
    return np.array([[ct, -st, 0, a], [st * ca, ct * ca, -sa, -d * sa], [st * sa, ct * sa, ca, d * ca], [0, 0, 0, 1]])


def fk(q, jac=False):
    """Flange pose (R, p); with jac=True also the body Jacobian J_b (rows omega, v)."""
    T = np.eye(4)
    zs, ps = [], []
    for i in range(7):
        T = T @ _T(A_[i], D_[i], AL[i], q[i])
        zs.append(T[:3, 2].copy())
        ps.append(T[:3, 3].copy())
    T = T @ np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, DF], [0, 0, 0, 1.0]])
    R, p = T[:3, :3], T[:3, 3]
    if not jac:
        return R, p
    Js = np.zeros((6, 7))
    for i in range(7):
        Js[:3, i] = zs[i]
        Js[3:, i] = np.cross(ps[i], zs[i])
    # J_b = Ad_{g^-1} J_s, Ad_{g^-1} = [[R^T, 0], [-R^T p^, R^T]]
    Jb = np.zeros((6, 7))
    Jb[:3] = R.T @ Js[:3]
    Jb[3:] = -R.T @ hat(p) @ Js[:3] + R.T @ Js[3:]
    return R, p, Jb


def kernel(J):
    """Generalized cross product of the six rows (n_j = (-1)^j det(J without column j)), normalized."""
    n = np.array([(-1) ** j * np.linalg.det(np.delete(J, j, axis=1)) for j in range(7)])
    return n / np.linalg.norm(n), np.linalg.norm(n)


# ---- coning circle path (Panda example) ----
R0 = expSO3(np.zeros(3))
_R, _p = fk(QR)
R0 = _R.copy()
CEN = np.array([0.40, 0.0, 0.45])
RAD, CONE = 0.12, 0.25
ACONE = expSO3(np.array([CONE, 0, 0]))
W0 = ACONE.T @ np.array([0, 0, 1.0]) - np.array([0, 0, 1.0])
SIG = np.sqrt(4 * WR * np.sin(CONE / 2) ** 2 + WP * RAD**2)
L = 2 * np.pi * SIG


def Rz(t):
    return expSO3(np.array([0, 0, t]))


def path(s):
    """gamma(s) = (R, p) and zeta_gamma(s) (unit W-speed body velocity)."""
    th = s / SIG
    R = R0 @ Rz(th) @ ACONE @ Rz(-th)
    p = CEN + RAD * np.array([np.cos(th), np.sin(th), 0])
    pd = RAD * np.array([-np.sin(th), np.cos(th), 0])
    z = np.concatenate([Rz(th) @ W0, R.T @ pd]) / SIG
    return R, p, z


def pose_err(q, s):
    """Product-metric Log of FK(q)^-1 gamma(s), in the flange frame (omega, v)."""
    R, p, J = fk(q, True)
    Rg, pg, _ = path(s)
    return np.concatenate([logSO3(R.T @ Rg), R.T @ (pg - p)]), J


def margin(q):
    return np.min(np.minimum(q - QMIN, QMAX - q))


def project(q, s, c=None, cval=0.0, tol=1e-13, itmax=30):
    """Newton onto the fibre over gamma(s); with c (7-vector) also impose c^T q = cval (square 7x7 system),
    otherwise minimum-norm (pseudoinverse) steps."""
    q = q.copy()
    for _ in range(itmax):
        e, J = pose_err(q, s)
        if c is None:
            dq = J.T @ np.linalg.solve(J @ J.T, e)
            if np.linalg.norm(e) < tol:
                return q
        else:
            r = cval - c @ q
            if np.linalg.norm(e) < tol and abs(r) < tol:
                return q
            dq = np.linalg.solve(np.vstack([J, c]), np.concatenate([e, [r]]))
        q = q + dq
    return q
