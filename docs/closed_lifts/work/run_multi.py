import numpy as np, time
from pandageo import *
from lifts import *
path = ConingCircle(); L = path.L
m = 250
s, Qr, gap = rest_lift(path, m); h = L/m
xr = kernel_vector(fk(Qr)[2])
# Multi-start: shift the rest lift along its fibres by a constant offset u0, then Newton (full energy).
res=[]
for u0 in (-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0):
    try:
        Q0 = retract(path, s, Qr, xr, np.full(m, u0))
        e,_ = pose_error(Q0, path, s)
        if np.abs(e).max() > 1e-8 or np.any(Q0<QMIN) or np.any(Q0>QMAX):
            print(f'u0 = {u0:+.1f}: start not admissible (max|e| {np.abs(e).max():.1e}, in box {np.all((Q0>QMIN)&(Q0<QMAX))})'); continue
        Qh,(a,b),hist = harmonic_section(path, s, Q0, h, kind='full', maxit=80)
        lam = np.linalg.eigvalsh(cyclic_tridiag_matrix(a,b))
        off = np.einsum('ma,ma->m', xr, Qh-Qr)
        inbox = np.all((Qh>QMIN)&(Qh<QMAX))
        print(f'u0 = {u0:+.1f}: E0 = {energy_full(Q0,h):.4f} -> E = {energy_full(Qh,h):.6f}, |grad| {hist[-1][2]:.1e}, iters {len(hist)}, offset range [{off.min():+.3f},{off.max():+.3f}], lowest Hill eig {lam[0]/h:.3f}, in box {inbox}')
    except Exception as ex:
        print('u0', u0, 'failed', ex)
