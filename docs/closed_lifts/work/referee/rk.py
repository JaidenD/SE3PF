"""Referee's independent Panda kinematics + coning circle (written from the paper / Franka DH table, not from
the agents' code). Twists (omega, v); body Jacobian in the flange frame."""
import numpy as np

A_ = [0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088]
D_ = [0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0]
AL = [0.0, -np.pi/2, np.pi/2, np.pi/2, -np.pi/2, np.pi/2, np.pi/2]
FL = 0.107
QMIN = np.array([-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973])
QMAX = np.array([2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973])
QR = np.array([0.0, -np.pi/4, 0.0, -3*np.pi/4, 0.0, np.pi/2, np.pi/4])
WR, WP = 0.01, 1.0


def rx(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def hat(w):
    return np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])


def expso3(w):
    th = np.linalg.norm(w)
    K = hat(w)
    if th < 1e-12:
        return np.eye(3) + K
    return np.eye(3) + np.sin(th)/th*K + (1-np.cos(th))/th**2*K@K


def logso3(R):
    c = np.clip((np.trace(R)-1)/2, -1, 1)
    th = np.arccos(c)
    w = np.array([R[2, 1]-R[1, 2], R[0, 2]-R[2, 0], R[1, 0]-R[0, 1]])
    if th < 1e-8:
        return 0.5*w
    return th/(2*np.sin(th))*w


def fk_all(q):
    """returns flange (R,p) and list of joint z-axes and origins in base frame"""
    R = np.eye(3); p = np.zeros(3)
    zs, os_ = [], []
    for i in range(7):
        # T_i = Rx(al) Tx(a) Rz(q) Tz(d)
        R1 = R @ rx(AL[i]); p1 = p + R1 @ np.array([A_[i], 0, 0])
        R = R1 @ rz(q[i]); p = p1 + R @ np.array([0, 0, D_[i]])
        zs.append(R[:, 2].copy()); os_.append(p.copy())
    pf = p + R @ np.array([0, 0, FL])
    return R, pf, zs, os_


def fk(q):
    R, p, _, _ = fk_all(q)
    return R, p


def jb(q):
    R, p, zs, os_ = fk_all(q)
    J = np.zeros((6, 7))
    for i in range(7):
        J[:3, i] = R.T @ zs[i]
        J[3:, i] = R.T @ np.cross(zs[i], p - os_[i])
    return J


def kernel(J):
    n = np.array([(-1)**j*np.linalg.det(np.delete(J, j, axis=1)) for j in range(7)])
    return n/np.linalg.norm(n)


# coning circle
R0 = fk(QR)[0]
CEN = np.array([0.40, 0.0, 0.45]); RAD = 0.12; CONE = 0.25
Acone = expso3(np.array([CONE, 0, 0]))
w0 = Acone.T @ np.array([0, 0, 1.0]) - np.array([0, 0, 1.0])
SIG = np.sqrt(WR*w0@w0 + WP*RAD**2)
L = 2*np.pi*SIG


def gamma(s):
    th = s/SIG
    return R0 @ rz(th) @ Acone @ rz(-th), CEN + RAD*np.array([np.cos(th), np.sin(th), 0])


def zeta(s):
    th = s/SIG
    R = gamma(s)[0]
    pd = RAD*np.array([-np.sin(th), np.cos(th), 0])
    return np.concatenate([rz(th) @ w0, R.T @ pd])/SIG


def pose_err(q, s):
    """body error log(gamma(s)^{-1} Pi(q)) as (omega, v-ish) 6-vector (rotation log, R_g^T dp)"""
    R, p = fk(q); Rg, pg = gamma(s)
    return np.concatenate([logso3(Rg.T @ R), Rg.T @ (p - pg)])


def project(q, s, x=None, iters=20, tol=1e-13):
    """Newton onto the fibre over gamma(s) keeping x^T(q - q0) = 0 (hyperplane slice through q0)"""
    q0 = q.copy(); q = q.copy()
    for _ in range(iters):
        e = pose_err(q, s)
        if np.linalg.norm(e) < tol:
            break
        J = jb(q)
        if x is None:
            dq = -np.linalg.lstsq(J, e, rcond=None)[0]
        else:
            Aaug = np.vstack([J, x]); rhs = -np.concatenate([e, [x @ (q - q0)]])
            dq = np.linalg.solve(Aaug, rhs)
        q = q + dq
    return q
