import numpy as np, time
from pandageo import *
from lifts import *
path = ConingCircle(); L = path.L
Qr = np.load('rest_lift_1000.npy')
t0=time.time()
t, C, closed = self_motion_curve(path, 0.0, Qr[0], ds=2e-3)
print('self-motion curve over gamma(0): closed in box?', closed, ' arclength range', t[0], t[-1], 'npts', len(t), 'time', time.time()-t0)
_,_,Jc = fk(C); sv = np.linalg.svd(Jc, compute_uv=False)[:, -1]
print('  sigma_min(J) along the arc: min', sv.min(), 'at t =', t[sv.argmin()])
marg = np.minimum(C-QMIN, QMAX-C).min(axis=1)
print('  limiting joints at the two ends:', np.argmin(np.minimum(C[0]-QMIN, QMAX-C[0])), np.argmin(np.minimum(C[-1]-QMIN, QMAX-C[-1])))
np.save('fibre0_t.npy', t); np.save('fibre0_C.npy', C)
# batch transport
def transport_batch(q0, n=2000, kfb=50.0):
    h = L/n; q=q0.copy(); s=0.0
    def f(q, s):
        e, J = pose_error(q, path, np.full(q.shape[0], s))
        return np.einsum('bij,bj->bi', pinv(J), path.velocity(s)[None,:] + kfb*e)
    inside = np.ones(q.shape[0], bool)
    for i in range(n):
        k1=f(q,s); k2=f(q+0.5*h*k1,s+0.5*h); k3=f(q+0.5*h*k2,s+0.5*h); k4=f(q+h*k3,s+h)
        q=q+(h/6)*(k1+2*k2+2*k3+k4); s+=h
        inside &= np.all((q>QMIN)&(q<QMAX),axis=1)
    return q, inside
idx = np.linspace(5, len(t)-6, 41).astype(int)
t0=time.time()
qL, inside = transport_batch(C[idx])
print('batch transport time', time.time()-t0)
# locate qL on the traced curve (nearest point, then local refinement by projection)
for k,i in enumerate(idx):
    d = np.linalg.norm(C - qL[k], axis=1); j = d.argmin()
    print(f'  t = {t[i]:+.3f}  ->  T(t) = {t[j]:+.3f}   T(t)-t = {t[j]-t[i]:+.4f}   dist to curve {d[j]:.1e}   stayed in box: {inside[k]}')
