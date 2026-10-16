"""Panda kinematics, the paper's coning-circle path, and the MuJoCo dynamics (M, C qdot), in numpy.

Conventions as in the project (se3pf): twists (omega, v), body Jacobian J_b in the flange frame, modified DH
(Craig), flange = joint-7 frame moved by 0.107 m along its z axis. Validated against the C++ library in
check_kinematics.py (data/fk_check.csv) and against MuJoCo's link frames.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

# modified DH: a_{i-1}, d_i, alpha_{i-1}
DH = np.array([[0.0, 0.333, 0.0],
               [0.0, 0.0, -np.pi / 2],
               [0.0, 0.316, np.pi / 2],
               [0.0825, 0.0, np.pi / 2],
               [-0.0825, 0.384, -np.pi / 2],
               [0.0, 0.0, np.pi / 2],
               [0.088, 0.0, np.pi / 2]])
FLANGE = 0.107
QMIN = np.array([-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973])
QMAX = np.array([2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973])
QR = np.array([0.0, -np.pi / 4, 0.0, -3 * np.pi / 4, 0.0, np.pi / 2, np.pi / 4])
W6 = np.diag([0.01, 0.01, 0.01, 1.0, 1.0, 1.0])


def hat(w):
    w = np.asarray(w)
    Z = np.zeros(w.shape[:-1] + (3, 3))
    Z[..., 0, 1] = -w[..., 2]; Z[..., 0, 2] = w[..., 1]
    Z[..., 1, 0] = w[..., 2]; Z[..., 1, 2] = -w[..., 0]
    Z[..., 2, 0] = -w[..., 1]; Z[..., 2, 1] = w[..., 0]
    return Z


def expSO3(w):
    w = np.asarray(w, dtype=float)
    th = np.linalg.norm(w, axis=-1)[..., None, None]
    K = hat(w)
    small = th < 1e-8
    ths = np.where(small, 1.0, th)
    A = np.where(small, 1.0 - th ** 2 / 6, np.sin(ths) / ths)
    B = np.where(small, 0.5 - th ** 2 / 24, (1 - np.cos(ths)) / ths ** 2)
    return np.eye(3) + A * K + B * (K @ K)


def logSO3(R):
    """Batched log on SO(3), angles < pi (robust near 0; fine for the small errors used here)."""
    tr = np.trace(R, axis1=-2, axis2=-1)
    c = np.clip((tr - 1) / 2, -1.0, 1.0)
    th = np.arccos(c)
    v = np.stack([R[..., 2, 1] - R[..., 1, 2], R[..., 0, 2] - R[..., 2, 0], R[..., 1, 0] - R[..., 0, 1]], -1)
    small = th < 1e-6
    s = np.where(small, 1.0, np.sin(th))
    f = np.where(small, 0.5 + th ** 2 / 12, th / (2 * s))
    return f[..., None] * v


def _tdh(a, d, al, th):
    """Modified DH transform for a batch of angles th (shape B) -> (B, 4, 4)."""
    ct, st = np.cos(th), np.sin(th)
    ca, sa = np.cos(al), np.sin(al)
    T = np.zeros(th.shape + (4, 4))
    T[..., 0, 0] = ct; T[..., 0, 1] = -st; T[..., 0, 3] = a
    T[..., 1, 0] = st * ca; T[..., 1, 1] = ct * ca; T[..., 1, 2] = -sa; T[..., 1, 3] = -d * sa
    T[..., 2, 0] = st * sa; T[..., 2, 1] = ct * sa; T[..., 2, 2] = ca; T[..., 2, 3] = d * ca
    T[..., 3, 3] = 1.0
    return T


def fk(q, jac=True):
    """Flange pose (R, p) and body Jacobian J_b (rows omega, v in the flange frame) for q of shape (..., 7)."""
    q = np.asarray(q, dtype=float)
    shp = q.shape[:-1]
    qf = q.reshape(-1, 7)
    B = qf.shape[0]
    T = np.tile(np.eye(4), (B, 1, 1))
    zs, ps = [], []
    for i in range(7):
        T = T @ _tdh(DH[i, 0], DH[i, 1], DH[i, 2], qf[:, i])
        zs.append(T[:, :3, 2].copy()); ps.append(T[:, :3, 3].copy())
    R = T[:, :3, :3]
    p = T[:, :3, 3] + FLANGE * T[:, :3, 2]
    if not jac:
        return R.reshape(shp + (3, 3)), p.reshape(shp + (3,))
    Z = np.stack(zs, -1)                      # (B, 3, 7)
    P = np.stack(ps, -1)
    Rt = np.transpose(R, (0, 2, 1))
    Jw = Rt @ Z
    Jv = Rt @ np.cross(Z, p[:, :, None] - P, axis=1)
    J = np.concatenate([Jw, Jv], axis=1)
    return R.reshape(shp + (3, 3)), p.reshape(shp + (3,)), J.reshape(shp + (6, 7))


def kernel_vector(J):
    """Unit kernel vector of a 6x7 J (batched), oriented by det([x^T; J]) > 0 (the library's kernelVector)."""
    J = np.asarray(J)
    n = np.stack([(-1) ** j * np.linalg.det(np.delete(J, j, axis=-1)) for j in range(7)], -1)
    return n / np.linalg.norm(n, axis=-1, keepdims=True)


def pose_error(R, p, Rt, pt):
    """e = (log(R^T Rt), R^T (pt - p)): product-metric Log_I of FK(q)^-1 target, in the flange frame."""
    RtR = np.swapaxes(R, -1, -2)
    w = logSO3(RtR @ Rt)
    v = np.einsum('...ij,...j->...i', RtR, pt - p)
    return np.concatenate([w, v], -1)


def constrained_ik(Rt, pt, a, b, q0, tol=1e-13, maxit=60, cap=0.2, extra=1):
    """Batched Newton on FK(q) = (Rt, pt), a^T q = b, from q0 (the library's constrainedIK). Returns q, ok.
    After the residual is below tol, `extra` more updates are taken (to reach roundoff level)."""
    q = np.array(q0, dtype=float, copy=True)
    done = 0
    for it in range(maxit):
        R, p, J = fk(q)
        e = pose_error(R, p, Rt, pt)
        r = b - np.einsum('...i,...i->...', a, q)
        err = np.maximum(np.abs(e).max(-1), np.abs(r))
        if np.all(err < tol):
            if done >= extra:
                return q, True
            done += 1
        A = np.concatenate([J, a[..., None, :]], -2)
        rhs = np.concatenate([e, r[..., None]], -1)
        dq = np.linalg.solve(A, rhs[..., None])[..., 0]
        nrm = np.linalg.norm(dq, axis=-1, keepdims=True)
        dq = dq * np.minimum(1.0, cap / np.maximum(nrm, 1e-300))
        q = q + dq
    R, p, J = fk(q)
    e = pose_error(R, p, Rt, pt)
    r = b - np.einsum('...i,...i->...', a, q)
    err = np.maximum(np.abs(e).max(-1), np.abs(r))
    return q, bool(np.all(err < 1e-10))


def joint_margin(q):
    q = np.asarray(q)
    m = np.minimum(q - QMIN, QMAX - q)
    return m.min(-1), m.argmin(-1)


# ---------------------------------------------------------------------------------------------------------------
# The paper's coning-circle path (closed form, as src/path.cpp)

def _meta():
    d = {}
    with open(os.path.join(DATA, "meta.csv")) as f:
        next(f)
        for line in f:
            k, v = line.strip().split(",")
            d[k] = float(v)
    return d


META = _meta()
R0 = np.array([[META["R0_%d%d" % (i, j)] for j in range(3)] for i in range(3)])
P0 = np.array([0.40, 0.0, 0.45])
RAD = 0.12
ALPHA = 0.25
A_CONE = expSO3(np.array([ALPHA, 0.0, 0.0]))
W0 = A_CONE.T @ np.array([0, 0, 1.0]) - np.array([0, 0, 1.0])
SIGMA = np.sqrt(4 * 0.01 * np.sin(ALPHA / 2) ** 2 + RAD ** 2)
L = 2 * np.pi * SIGMA


def Rz(th):
    th = np.asarray(th, dtype=float)
    c, s = np.cos(th), np.sin(th)
    R = np.zeros(th.shape + (3, 3))
    R[..., 0, 0] = c; R[..., 0, 1] = -s; R[..., 1, 0] = s; R[..., 1, 1] = c; R[..., 2, 2] = 1
    return R


def path_pose(s):
    th = np.asarray(s, dtype=float) / SIGMA
    R = R0 @ Rz(th) @ A_CONE @ Rz(-th)
    p = P0 + RAD * np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], -1)
    return R, p


def path_vel(s):
    th = np.asarray(s, dtype=float) / SIGMA
    R, _ = path_pose(s)
    w = Rz(th) @ W0
    pd = RAD * np.stack([-np.sin(th), np.cos(th), np.zeros_like(th)], -1)
    v = np.einsum('...ji,...j->...i', R, pd)
    return np.concatenate([w, v], -1) / SIGMA


def path_acc(s):
    th = np.asarray(s, dtype=float) / SIGMA
    R, _ = path_pose(s)
    w = Rz(th) @ W0
    dw = Rz(th) @ np.cross(np.array([0, 0, 1.0]), W0)
    pd = RAD * np.stack([-np.sin(th), np.cos(th), np.zeros_like(th)], -1)
    pdd = -RAD * np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], -1)
    v = np.einsum('...ji,...j->...i', R, pd)
    dv = np.einsum('...ji,...j->...i', R, pdd) - np.cross(w, v)
    return np.concatenate([dw, dv], -1) / SIGMA ** 2


def load_rest(N):
    """Rest-posture lift nodes (Lift(path, q_r, N) of the library) and the chart vectors x_k."""
    a = np.loadtxt(os.path.join(DATA, "rest_%d.csv" % N), delimiter=",", skiprows=1)
    s = a[:, 1]
    q = a[:, 26:33]
    x = a[:, 33:40]
    return s, q, x


def load_frame(N):
    a = np.loadtxt(os.path.join(DATA, "frame_%d.csv" % N), delimiter=",", skiprows=1)
    s = a[:, 1]
    Nf = a[:, 2:32].reshape(-1, 5, 6).transpose(0, 2, 1)    # column-major 6x5
    dN = a[:, 32:62].reshape(-1, 5, 6).transpose(0, 2, 1)
    return s, Nf, dN


# ---------------------------------------------------------------------------------------------------------------
# MuJoCo dynamics: M(q) (armature included) and C(q, v) v with gravity off

MJ_XML = "/home/jaiden/Desktop/SE3_pathFollowing/third_party/mujoco_menagerie/franka_emika_panda/panda_nohand.xml"


class MjDyn:
    def __init__(self):
        import mujoco
        self.mj = mujoco
        self.m = mujoco.MjModel.from_xml_path(MJ_XML)
        self.m.opt.gravity[:] = 0.0
        self.d = mujoco.MjData(self.m)
        self.Mbuf = np.zeros((7, 7))

    def M(self, q):
        mj, m, d = self.mj, self.m, self.d
        d.qpos[:] = q
        d.qvel[:] = 0.0
        mj.mj_kinematics(m, d)
        mj.mj_comPos(m, d)
        mj.mj_crb(m, d)
        mj.mj_fullM(m, d, self.Mbuf)
        return self.Mbuf.copy()

    def Mbatch(self, Q):
        Q = np.atleast_2d(Q)
        return np.array([self.M(q) for q in Q])

    def Cv(self, q, v):
        """C(q, v) v = qfrc_bias with gravity zero (Christoffel terms of M, armature has none)."""
        mj, m, d = self.mj, self.m, self.d
        d.qpos[:] = q
        d.qvel[:] = v
        mj.mj_kinematics(m, d)
        mj.mj_comPos(m, d)
        mj.mj_comVel(m, d)
        mj.mj_rne(m, d, 0, d.qfrc_bias)
        return d.qfrc_bias.copy()

    def Cv_batch(self, Q, V):
        return np.array([self.Cv(q, v) for q, v in zip(Q, V)])
