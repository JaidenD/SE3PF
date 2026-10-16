"""(iii) Sensitivity of the lifts to a 1 mm change of the Panda's d5 (0.384 -> 0.385 m), path fixed in the world."""
import numpy as np
import kin
from kinb import *
from harmonic import Fibres, load_rest, solve, energy
from midpoint import exit_arc, flow_vec

ss, Qrest = load_rest(500); F = Fibres(ss)
QA = np.load('harm_euclid_500.npy'); mp = np.load('midpoint_500.npz'); Qm0 = mp['Qm']
QB = np.load('harm_euclid_500_start-3.npy')

def rest_lift(Q):
    Q = Q.copy()
    for it in range(200):
        e, J = pose_err_b(Q, F.Rg, F.pg)
        JJt = J @ np.transpose(J, (0, 2, 1))
        d = QR - Q
        Jd = np.einsum('mij,mj->mi', J, d)
        dq = pinv_apply(J, e) + d - pinv_apply(J, Jd)
        Q = Q + dq
        if np.abs(dq).max() < 1e-13: break
    return Q

kin.D_[4] = 0.385
Qr1 = rest_lift(Qrest)
QA1, EA1, el, w = solve(F, QA.copy())
QB1, EB1, elB, wB = solve(F, QB.copy())
b, jb = exit_arc(F, Qr1, +1); a, ja = exit_arc(F, Qr1, -1); a = -a
Qm1 = flow_vec(F, Qr1.copy(), 0.5 * (a + b), 0.01)
for name, Q0, Q1 in [('rest-posture lift', Qrest, Qr1), ('harmonic A', QA, QA1), ('harmonic B', QB, QB1), ('midpoint lift', Qm0, Qm1)]:
    print(f"{name:18s}: max_k |q1_k - q0_k| = {np.linalg.norm(Q1 - Q0, axis=1).max():.3e} rad")
print('energies after: A %.6f  B %.6f' % (EA1, EB1), ' EL resid', el, elB)
