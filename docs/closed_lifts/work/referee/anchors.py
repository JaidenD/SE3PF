"""Referee: check the referee kinematics rk.py against the design-doc anchors and finite differences."""
import numpy as np, rk
R, p = rk.fk(np.zeros(7)); print("FK(0) p =", p, " R diag =", np.diag(R), " offdiag max", np.abs(R-np.diag(np.diag(R))).max())
R, p = rk.fk(rk.QR); print("FK(q_r) p =", p, "R=\n", R)
print("L =", rk.L, " sigma =", rk.SIG, " gamma(0).p =", rk.gamma(0.0)[1])
# unit W-speed and zeta vs finite difference
W = np.diag([rk.WR]*3+[rk.WP]*3)
rng = np.random.default_rng(0)
err = 0
for s in rng.random(5)*rk.L:
    z = rk.zeta(s); hh = 1e-6
    Rp, pp = rk.gamma(s+hh); Rm, pm = rk.gamma(s-hh); Rg, pg = rk.gamma(s)
    zfd = np.concatenate([rk.logso3(Rm.T@Rp)/(2*hh), Rg.T@(pp-pm)/(2*hh)])
    err = max(err, np.abs(z-zfd).max()); assert abs(z@W@z-1) < 1e-12
print("zeta vs FD:", err)
# J_b vs FD
err = 0
for t in range(3):
    q = rk.QMIN + (rk.QMAX-rk.QMIN)*rng.random(7); J = rk.jb(q); R, p = rk.fk(q)
    for j in range(7):
        e = np.zeros(7); e[j] = 1e-6
        Rp, pp = rk.fk(q+e); Rm, pm = rk.fk(q-e)
        col = np.concatenate([rk.logso3(Rm.T@Rp)/2e-6, R.T@(pp-pm)/2e-6])
        err = max(err, np.abs(col-J[:, j]).max())
print("J_b vs FD:", err)
rest = np.load("../rest_lift_1000.npy"); print("rest lift q*(0) =", np.round(rest[0], 5))
m = len(rest); h = rk.L/m
print("rest lift FK residual max:", max(np.linalg.norm(rk.pose_err(rest[j], j*h)) for j in range(0, m, 7)))
