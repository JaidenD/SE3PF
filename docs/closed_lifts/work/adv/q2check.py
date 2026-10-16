import numpy as np
from kin import QMIN, QMAX
from selfmotion import wrap
for k in [0, 50, 100, 150, 200, 238, 250, 300, 350, 400, 450]:
    d = np.load('sm_%d.npz' % k); TH, QQ, a, b = d['TH'], d['QQ'], float(d['a']), float(d['b'])
    ell = float(d['ell'])
    # admissible arc: theta in [a, b] (a < 0 means wrap from the end)
    th = np.where(TH > ell / 2, TH - ell, TH)
    sel = (th >= a) & (th <= b)
    o = np.argsort(th[sel]); q2 = QQ[sel][o][:, 1]; t = th[sel][o]
    dq2 = np.diff(q2)
    print(f"k={k:3d}: q2 along admissible arc from {q2[0]:+.3f} to {q2[-1]:+.3f}; monotone: {bool(np.all(dq2 > 0) or np.all(dq2 < 0))}, "
          f"min |dq2/dtheta| {np.min(np.abs(dq2 / np.diff(t))):.3f}; q2 at rest lift {QQ[0,1]:+.3f}")
