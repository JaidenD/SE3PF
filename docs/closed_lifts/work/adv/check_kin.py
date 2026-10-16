import numpy as np
from kin import *
D = '../num/data/'
fkc = np.loadtxt(D + 'fk_check.csv', delimiter=',')
e1 = e2 = 0
for row in fkc:
    q = row[:7]; G = row[7:19].reshape(3, 4); Jc = row[19:].reshape(6, 7)
    R, p, J = fk(q, True)
    e1 = max(e1, np.abs(R - G[:, :3]).max(), np.abs(p - G[:, 3]).max()); e2 = max(e2, np.abs(J - Jc).max())
print('FK err', e1, 'Jb err', e2)
print('L', L, 'sigma', SIG)
rest = np.genfromtxt(D + 'rest_500.csv', delimiter=',', names=True)
qs = np.stack([rest['q%d' % i] for i in range(7)], 1)
xs = np.stack([rest['x%d' % i] for i in range(7)], 1)
ss = rest['s']
em = ex = ez = 0
for k in range(0, 500, 7):
    e, J = pose_err(qs[k], ss[k]); em = max(em, np.abs(e).max())
    x, _ = kernel(J); ex = max(ex, np.abs(x - xs[k]).max())
    z = path(ss[k])[2]; ez = max(ez, np.abs(z - np.array([rest['z%d' % i][k] for i in range(6)])).max())
print('rest lift pose err', em, 'kernel err', ex, 'zeta err', ez)
print('q*(0)', qs[0])
