import numpy as np, sys
from pandageo import *
from lifts import *
path = ConingCircle(); L = path.L
m = 250
s, Qr, gap = rest_lift(path, m); h = L/m
xr = kernel_vector(fk(Qr)[2])
def offs(Q): return np.einsum('ma,ma->m', xr, Q-Qr)
which = sys.argv[1]
if which == 'a':
    Q0 = retract(path, s, Qr, xr, np.full(m, 2.0))
    Q,(a,b),hist = harmonic_section(path, s, Q0, h, kind='full', maxit=80)
    tr,_ = monodromy_trace(a,b); lam=np.linalg.eigvalsh(cyclic_tridiag_matrix(a,b))/h; o=offs(Q)
    marg=np.minimum(Q-QMIN,QMAX-Q)
    print(f'second local min: E = {energy_full(Q,h):.6f}, offset [{o.min():+.3f},{o.max():+.3f}], Hill eig {lam[0]:.3f}, trace {tr:.4f}, multipliers {np.roots([1,-tr,1])}, min margin {marg.min():.4f} (joint {np.unravel_index(marg.argmin(),marg.shape)[1]+1}), min sigma_min {np.linalg.svd(fk(Q)[2],compute_uv=False)[:,-1].min():.4f}')
    # vertical-energy lift: Floquet trace of its Hessian
    Qv,(av,bv),hv = harmonic_section(path, s, Qr.copy(), h, kind='vertical', maxit=60, tol=1e-8)
    trv,_ = monodromy_trace(av,bv)
    print(f'vertical-energy lift: E_v = {energy_vertical(Qv,h):.6f}, trace {trv:.5f}, multipliers {np.roots([1,-trv,1])}')
    # Prop E corrector constant: delta = (1/k) (x . q_f'') / (1 + (q_f - q_r) . kappa)
    x,kap,J = fibre_jets(path, s, Qr)
    qdd = (np.roll(Qr,-1,0) - 2*Qr + np.roll(Qr,1,0))/h**2
    num = np.einsum('ma,ma->m', x, qdd); den = 1 + np.einsum('ma,ma->m', Qr-QR, kap)
    print('Prop E corrector: max |x.q_f\'\'/Hess| =', np.abs(num/den).max(), ' min Hess^V f =', den.min(), ' max|q_f\'\'^V| =', np.abs(num).max())
    # twist condition on the lifts
    print('twist: min x_j.x_{j+1} on rest lift', np.einsum('ma,ma->m', x, np.roll(x,-1,0)).min())
    # stretch: d/ds log rho = -kappa . Y along the rest lift, Y = J^+ zeta
    Y = np.einsum('mab,mb->ma', pinv(J), path.velocity(s))
    g = -np.einsum('ma,ma->m', kap, Y)
    print('stretch rate -kappa.Y along rest lift: range', g.min(), g.max(), ' integral over a lap', g.sum()*h, ' |kappa| range', np.linalg.norm(kap,axis=1).min(), np.linalg.norm(kap,axis=1).max())
elif which == 'd':
    for sj in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7):
        j = int(round(sj/L*m)) % m
        t, C, closed = self_motion_curve(path, s[j], Qr[j], ds=4e-3)
        lim = [int(np.argmin(np.minimum(C[i]-QMIN, QMAX-C[i]))) + 1 for i in (0, -1)]
        sv = np.linalg.svd(fk(C)[2], compute_uv=False)[:,-1]
        print(f's = {s[j]:.3f}: admissible self-motion component through q_rest: closed circle {closed}, arclength [{t[0]:.2f}, {t[-1]:.2f}], ends at limits of joints {lim}, min sigma_min on it {sv.min():.3f}')
