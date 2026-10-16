"""Referee: the Panda shoulder involution sigma(q) = (q1 + pi, -q2, q3 + pi, q4, ..., q7) (mod 2 pi) leaves the forward
kinematics unchanged (joints 1-3 form a ZYZ spherical shoulder; the d3 offset is along z3). It is an isometry of the
Euclidean joint metric, so it maps closed lifts to closed lifts with the same K = I energy. Check it on random q, and
apply it to the stored lifts (A, B, C, rest lift, midsection): is the image admissible, and what is its energy?"""
import numpy as np, rk
def wrap(a): return (a+np.pi) % (2*np.pi) - np.pi
def sigma(Q):
    Q = np.array(Q, float).copy(); Q[..., 0] = wrap(Q[..., 0]+np.pi); Q[..., 1] = -Q[..., 1]; Q[..., 2] = wrap(Q[..., 2]+np.pi); return Q
rng = np.random.default_rng(0); err = 0
for _ in range(20):
    q = rk.QMIN+(rk.QMAX-rk.QMIN)*rng.random(7); R1, p1 = rk.fk(q); R2, p2 = rk.fk(sigma(q))
    err = max(err, np.abs(R1-R2).max(), np.abs(p1-p2).max())
print("max |FK(sigma q) - FK(q)| over 20 random q:", err)
def E(Q): h = rk.L/len(Q); return np.sum((np.roll(Q, -1, 0)-Q)**2)/(2*h)
for name, Q in [("A", np.load("../harmonic_full_1000.npy")), ("B", np.load("../num/data/far_I_start-2.npz")["q"]),
                ("C", np.load("../num/data/far_I_start+2.npz")["q"]), ("rest", np.load("../rest_lift_1000.npy")),
                ("midsection", np.load("midsection_I_m100.npy"))]:
    S = sigma(Q); jumps = np.abs(np.diff(np.vstack([S, S[:1]]), axis=0)).max()
    mg = np.minimum(S-rk.QMIN, rk.QMAX-S).min(axis=1)
    print(f"{name:10s}: E(q) = {E(Q):.6f}, E(sigma q) = {E(S):.6f} (max node jump {jumps:.3f}), sigma(q) admissible at "
          f"{np.mean(mg > 0)*100:.0f}% of nodes, min margin of sigma(q) {mg.min():.4f}, q1 range [{Q[:,0].min():.2f},{Q[:,0].max():.2f}], "
          f"q3 range [{Q[:,2].min():.2f},{Q[:,2].max():.2f}]")
