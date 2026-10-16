"""Validate the numpy kinematics/path/dynamics against the C++ library dumps and MuJoCo."""
import numpy as np
import mujoco
from pandakin import *

a = np.loadtxt(os.path.join(DATA, "fk_check.csv"), delimiter=",")
q = a[:, :7]
G = a[:, 7:19].reshape(-1, 3, 4)
J = a[:, 19:61].reshape(-1, 6, 7)
R, p, Jp = fk(q)
print("FK vs C++: R %.2e  p %.2e  J_b %.2e" % (abs(R - G[:, :, :3]).max(), abs(p - G[:, :, 3]).max(), abs(Jp - J).max()))

for N in (250, 2000):
    arr = np.loadtxt(os.path.join(DATA, "rest_%d.csv" % N), delimiter=",", skiprows=1)
    s = arr[:, 1]
    Rg, pg = path_pose(s)
    print("path N=%d: R %.2e p %.2e zeta %.2e zeta' %.2e" % (
        N, abs(Rg.reshape(-1, 9) - arr[:, 2:11]).max(), abs(pg - arr[:, 11:14]).max(),
        abs(path_vel(s) - arr[:, 14:20]).max(), abs(path_acc(s) - arr[:, 20:26]).max()))
    qs = arr[:, 26:33]
    Rq, pq, Jq = fk(qs)
    e = pose_error(Rq, pq, Rg, pg)
    x = kernel_vector(Jq)
    print("   rest lift: |FK(q*) - gamma| %.2e, |x_py - x_cpp| %.2e" % (abs(e).max(), abs(x - arr[:, 33:40]).max()))
print("L python %.15f  C++ %.15f" % (L, META["L"]))

# MuJoCo frames: link7 + 0.107 z must equal the DH flange
m = mujoco.MjModel.from_xml_path(MJ_XML)
d = mujoco.MjData(m)
b7 = m.body("link7").id
err = 0.0
for qq in q:
    d.qpos[:] = qq
    mujoco.mj_kinematics(m, d)
    Rm = d.xmat[b7].reshape(3, 3)
    pm = d.xpos[b7] + FLANGE * Rm[:, 2]
    Rd, pd = fk(qq, jac=False)
    err = max(err, abs(Rm - Rd).max(), abs(pm - pd).max())
print("MuJoCo link7+0.107z vs DH flange: %.2e" % err)
print("armature", m.dof_armature, "damping", m.dof_damping, "frictionloss", m.dof_frictionloss)

# C(q,v)v from RNE vs the Christoffel identity Cv = Mdot v - 1/2 grad(v^T M v)
dyn = MjDyn()
rng = np.random.default_rng(0)
worst = 0
for t in range(5):
    qq = QR + 0.5 * rng.standard_normal(7)
    v = rng.standard_normal(7)
    h = 1e-5
    Md = (dyn.M(qq + h * v) - dyn.M(qq - h * v)) / (2 * h)
    gr = np.array([(v @ dyn.M(qq + h * e) @ v - v @ dyn.M(qq - h * e) @ v) / (2 * h) for e in np.eye(7)])
    Cv_fd = Md @ v - 0.5 * gr
    Cv = dyn.Cv(qq, v)
    worst = max(worst, abs(Cv - Cv_fd).max() / abs(Cv).max())
print("RNE C(q,v)v vs Christoffel identity (FD of M): rel err %.2e" % worst)
d.qpos[:] = QR
mujoco.mj_forward(m, d)
Mf = np.zeros((7, 7)); mujoco.mj_fullM(m, d, Mf)
print("M(q_r) via MjDyn vs mj_forward: %.2e ; diag M(q_r) =" % abs(Mf - dyn.M(QR)).max(), np.round(np.diag(Mf), 4))
