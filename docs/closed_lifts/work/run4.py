import numpy as np
from pandageo import *
from lifts import *
path = ConingCircle(); L = path.L; m = 250
s, Qr, gap = rest_lift(path, m); h = L/m
Qv,(a,b),hist = harmonic_section(path, s, Qr.copy(), h, kind='vertical', maxit=60, tol=1e-8)
print('vertical: Newton steps', len(hist)-1, ' E_v', energy_vertical(Qv,h), ' E_full', energy_full(Qv,h))
print('  max_s |q_v - q_rest| =', np.linalg.norm(Qv-Qr,axis=1).max(), ' min sigma_min', np.linalg.svd(fk(Qv)[2],compute_uv=False)[:,-1].min())
np.save('vertical_250.npy', Qv)
# monotonicity of self-motion speed along the vertical lift (Cor. 4.4): v = x.q'
x = kernel_vector(fk(Qv)[2]); dq = (np.roll(Qv,-1,0)-np.roll(Qv,1,0))/(2*h)
v = np.einsum('ma,ma->m', x, dq)
print('  self-motion speed v along vertical lift: min', v.min(), 'max', v.max(), '(Cor. 4.4: never vanishes)')
Qh = np.load('min_from_0.0.npy'); xh = kernel_vector(fk(Qh)[2]); vh = np.einsum('ma,ma->m', xh, (np.roll(Qh,-1,0)-np.roll(Qh,1,0))/(2*h))
print('  for comparison, full-energy lift: v range', vh.min(), vh.max())
