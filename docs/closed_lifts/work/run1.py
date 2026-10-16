import numpy as np, time
from pandageo import *
from lifts import *
path = ConingCircle()
L = path.L
t0=time.time()
m = 1000
s, Qr, gap = rest_lift(path, m)
print('rest lift: m', m, 'closure gap', gap, 'q*(0)', np.round(Qr[0],5), 'time', time.time()-t0)
marg = np.minimum(Qr-QMIN, QMAX-Qr)
print('min joint margin', marg.min(), 'joint', np.unravel_index(marg.argmin(), marg.shape))
np.save('rest_lift_1000.npy', Qr)
# one lap of pseudoinverse transport from q*(0)
for n in (2000, 4000):
    t0=time.time()
    q1 = transport(path, Qr[0], 0.0, L, n=n)
    x0 = kernel_vector(fk(Qr[0])[2])
    print('pinv transport n', n, ': |q(L)-q(0)|', np.linalg.norm(q1-Qr[0]), ' x0^T(q(L)-q(0))', x0@(q1-Qr[0]), 'time', time.time()-t0)
np.save('pinv_one_lap.npy', q1)
