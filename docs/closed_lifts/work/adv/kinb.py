"""Batched (over configurations) Panda FK / body Jacobian / kernel vector / pose error."""
import numpy as np
from kin import A_, D_, AL, DF, QMIN, QMAX, QR, path, SIG, L, WR, WP

def _Tb(i, th):
    a, d, al = A_[i], D_[i], AL[i]
    ca, sa = np.cos(al), np.sin(al)
    ct, st = np.cos(th), np.sin(th)
    M = th.shape[0]
    T = np.zeros((M, 4, 4))
    T[:, 0, 0] = ct; T[:, 0, 1] = -st; T[:, 0, 3] = a
    T[:, 1, 0] = st * ca; T[:, 1, 1] = ct * ca; T[:, 1, 2] = -sa; T[:, 1, 3] = -d * sa
    T[:, 2, 0] = st * sa; T[:, 2, 1] = ct * sa; T[:, 2, 2] = ca; T[:, 2, 3] = d * ca
    T[:, 3, 3] = 1
    return T

def hatb(w):
    M = w.shape[0]; H = np.zeros((M, 3, 3))
    H[:, 0, 1] = -w[:, 2]; H[:, 0, 2] = w[:, 1]; H[:, 1, 0] = w[:, 2]
    H[:, 1, 2] = -w[:, 0]; H[:, 2, 0] = -w[:, 1]; H[:, 2, 1] = w[:, 0]
    return H

def fkb(Q):
    M = Q.shape[0]
    T = np.tile(np.eye(4), (M, 1, 1))
    zs = np.zeros((M, 7, 3)); ps = np.zeros((M, 7, 3))
    for i in range(7):
        T = T @ _Tb(i, Q[:, i])
        zs[:, i] = T[:, :3, 2]; ps[:, i] = T[:, :3, 3]
    R = T[:, :3, :3]; p = T[:, :3, 3] + DF * T[:, :3, 2]
    Jw = np.transpose(zs, (0, 2, 1))                      # M,3,7
    Jv = np.transpose(np.cross(ps, zs), (0, 2, 1))
    RT = np.transpose(R, (0, 2, 1))
    Jb = np.zeros((M, 6, 7))
    Jb[:, :3] = RT @ Jw
    Jb[:, 3:] = -RT @ hatb(p) @ Jw + RT @ Jv
    return R, p, Jb

def logSO3b(R):
    c = np.clip((np.trace(R, axis1=1, axis2=2) - 1) / 2, -1, 1)
    th = np.arccos(c)
    v = np.stack([R[:, 2, 1] - R[:, 1, 2], R[:, 0, 2] - R[:, 2, 0], R[:, 1, 0] - R[:, 0, 1]], 1)
    f = np.where(th < 1e-8, 0.5, th / (2 * np.sin(np.maximum(th, 1e-300))))
    return f[:, None] * v

def path_b(svec):
    out = [path(s) for s in svec]
    return (np.stack([o[0] for o in out]), np.stack([o[1] for o in out]), np.stack([o[2] for o in out]))

def pose_err_b(Q, Rg, pg):
    R, p, J = fkb(Q)
    RT = np.transpose(R, (0, 2, 1))
    e = np.concatenate([logSO3b(RT @ Rg), np.einsum('mij,mj->mi', RT, pg - p)], 1)
    return e, J

def kernel_b(J):
    n = np.stack([(-1) ** j * np.linalg.det(np.delete(J, j, axis=2)) for j in range(7)], 1)
    return n / np.linalg.norm(n, axis=1, keepdims=True)

def pinv_apply(J, v):
    return np.einsum('mji,mj->mi', J, np.linalg.solve(J @ np.transpose(J, (0, 2, 1)), v[..., None])[..., 0])

def project_b(Q, Rg, pg, C=None, cval=None, tol=1e-13, itmax=40):
    """Batched Newton onto fibres; with C (M,7) and cval (M,) also impose C q = cval."""
    Q = Q.copy()
    for _ in range(itmax):
        e, J = pose_err_b(Q, Rg, pg)
        if C is None:
            if np.abs(e).max() < tol: break
            Q = Q + pinv_apply(J, e)
        else:
            r = cval - np.einsum('mi,mi->m', C, Q)
            if np.abs(e).max() < tol and np.abs(r).max() < tol: break
            A = np.concatenate([J, C[:, None, :]], 1)
            Q = Q + np.linalg.solve(A, np.concatenate([e, r[:, None]], 1)[..., None])[..., 0]
    return Q

def feasible(Q):
    return np.all((Q > QMIN) & (Q < QMAX), axis=-1)
