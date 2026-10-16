"""E2: holonomy diffeomorphism h: F_0 -> F_0 of the Euclidean (pseudoinverse) connection over one lap, on the
full self-motion circle (joints on T^7, limits ignored)."""
import numpy as np
from kinb import *
from selfmotion import wrap
d = np.load('sm_0.npz'); TH, QQ, ell = d['TH'], d['QQ'], float(d['ell'])
M = 140
idx = np.searchsorted(TH, np.linspace(0, TH[-1], M, endpoint=False))
Q = QQ[idx].copy(); th0 = TH[idx]
nst = 2000; h = L / nst
ss = np.linspace(0, L, 2 * nst + 1)
Rg, pg, Z = path_b(ss)
Ke = 50.0
def rhs(Q, i):
    e, J = pose_err_b(Q, np.broadcast_to(Rg[i], (len(Q), 3, 3)), np.broadcast_to(pg[i], (len(Q), 3)))
    return pinv_apply(J, Z[i][None, :] + Ke * e)
for n in range(nst):
    i = 2 * n
    k1 = rhs(Q, i); k2 = rhs(Q + 0.5 * h * k1, i + 1); k3 = rhs(Q + 0.5 * h * k2, i + 1); k4 = rhs(Q + h * k3, i + 2)
    Q = Q + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
# locate images on the traced circle F_0: nearest traced point (wrapped) then local linear refinement
def locate(q):
    dd = np.linalg.norm(wrap(QQ - q), axis=1); j = int(np.argmin(dd))
    j2 = min(j + 1, len(TH) - 1); t = QQ[j2] - QQ[j]
    a = np.clip(wrap(q - QQ[j]) @ t / (t @ t), -1, 1)
    return TH[j] + a * (TH[j2] - TH[j]), dd[j]
img = np.array([locate(q) for q in Q])
th1, dist = img[:, 0], img[:, 1]
disp = (th1 - th0 + ell / 2) % ell - ell / 2      # displacement h(theta) - theta in (-ell/2, ell/2]
print('fibre length', ell, ' max off-fibre distance of images', dist.max())
print('displacement h(theta)-theta over F_0: min %.4f max %.4f mean %.4f (rad of Euclidean arc length)' % (disp.min(), disp.max(), disp.mean()))
print('sign changes of displacement (fixed points of h):', np.sum(np.sign(disp) != np.sign(np.roll(disp, 1))))
for t, dsp in zip(th0[::7], disp[::7]): print('  theta %.3f  h-theta %+.4f' % (t, dsp))
np.savez('holonomy_euclid.npz', th0=th0, th1=th1, disp=disp, ell=ell)
print('dist quantiles', np.quantile(dist, [0, 0.5, 0.9, 0.95, 1.0]))
bad = np.where(dist > 1e-3)[0]; print('bad idx', bad, 'theta0', th0[bad])
e, J = pose_err_b(Q, np.broadcast_to(Rg[-1], (len(Q), 3, 3)), np.broadcast_to(pg[-1], (len(Q), 3)))
print('end pose err max', np.abs(e).max())
from kinb import fkb
print('sigma_min along circle samples min', np.linalg.svd(fkb(QQ[idx])[2], compute_uv=False)[:, -1].min())
