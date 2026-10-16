"""Referee: K = M (MuJoCo menagerie Panda, armature 0.1) midpoint-rule energies of stored lifts, and the M-energy of
the sigma-images (the shoulder involution is NOT an M-isometry)."""
import sys, numpy as np
sys.path.insert(0, "../num"); from pandakin import MjDyn, L
import rk
dyn = MjDyn()
def EM(Q):
    h = L/len(Q); D = np.roll(Q, -1, 0)-Q; Mm = dyn.Mbatch(0.5*(Q+np.roll(Q, -1, 0)))
    return np.einsum("ki,kij,kj->", D, Mm, D)/(2*h)
def wrap(a): return (a+np.pi) % (2*np.pi)-np.pi
def sigma(Q): Q = Q.copy(); Q[:, 0] = wrap(Q[:, 0]+np.pi); Q[:, 1] = -Q[:, 1]; Q[:, 2] = wrap(Q[:, 2]+np.pi); return Q
for name, Q in [("rest (1000)", np.load("../rest_lift_1000.npy")), ("K=M near-rest (1000)", np.load("../num/data/harm_M_1000_u0.npz")["q"]),
                ("K=M far start-2 (1000)", np.load("../num/data/far_M_start-2.npz")["q"]), ("K=I A (1000)", np.load("../harmonic_full_1000.npy")),
                ("K=I B (1000)", np.load("../num/data/far_I_start-2.npz")["q"]), ("midsection (100)", np.load("midsection_I_m100.npy"))]:
    print(f"{name:24s} E_M = {EM(Q):.5f}   E_M(sigma q) = {EM(sigma(Q)):.5f}")
