import numpy as np, time, sys
from pandageo import *
from lifts import *
path = ConingCircle(); L = path.L
m = 250
s, Qr, gap = rest_lift(path, m); h = L/m
xr = kernel_vector(fk(Qr)[2])
def offs(Q): return np.einsum('ma,ma->m', xr, Q-Qr)
def pure_newton(Q, kind='full', k=0.0, maxit=40, tol=1e-10):
    for it in range(maxit):
        x,kap,J = fibre_jets(path, s, Q)
        g,a,b = grad_hess_full(Q,x,kap,h,k)
        if np.abs(g).max()<tol: break
        du = -np.linalg.solve(cyclic_tridiag_matrix(a,b), g)
        nd = np.abs(du).max()
        if nd>0.3: du*=0.3/nd
        Q = retract(path, s, Q, x, du)
    x,kap,J = fibre_jets(path, s, Q); g,a,b = grad_hess_full(Q,x,kap,h,k)
    return Q, g, a, b
which = sys.argv[1]
if which == 'mp':
    for u0 in (1.2, 1.4, 1.6):
        Q0 = retract(path, s, Qr, xr, np.full(m, u0))
        Q, g, a, b = pure_newton(Q0)
        H = cyclic_tridiag_matrix(a,b); lam = np.linalg.eigvalsh(H)/h
        tr,_ = monodromy_trace(a,b)
        o = offs(Q)
        print(f'start offset {u0}: |grad| {np.abs(g).max():.1e}, E = {energy_full(Q,h):.6f}, offset [{o.min():+.3f},{o.max():+.3f}], Hill eigs {lam[:3]}, Morse index {np.sum(lam<0)}, monodromy trace {tr:.4f}, in box {np.all((Q>QMIN)&(Q<QMAX))}')
    # the two minimizers: traces
    for u0 in (0.0, 2.0):
        Q0 = retract(path, s, Qr, xr, np.full(m, u0)); Q,g,a,b = pure_newton(Q0)
        tr,_ = monodromy_trace(a,b); lam=np.linalg.eigvalsh(cyclic_tridiag_matrix(a,b))/h; o=offs(Q)
        marg = np.minimum(Q-QMIN, QMAX-Q)
        print(f'minimizer from {u0}: E = {energy_full(Q,h):.6f}, offset [{o.min():+.3f},{o.max():+.3f}], Hill eig {lam[0]:.3f}, trace {tr:.4f}, min margin {marg.min():.4f} (joint {np.unravel_index(marg.argmin(),marg.shape)[1]+1}), min sigma_min {np.linalg.svd(fk(Q)[2],compute_uv=False)[:,-1].min():.4f}')
        np.save(f'min_from_{u0}.npy', Q)
elif which == 'k':
    Q = Qr.copy()
    for k in (1e4, 1e3, 300, 100, 30, 10, 3, 1, 0.3, 0.0):
        Q, g, a, b = pure_newton(Q, k=k)
        lam = np.linalg.eigvalsh(cyclic_tridiag_matrix(a,b))/h
        d = np.linalg.norm(Q-Qr,axis=1).max()
        print(f'k = {k:8.1f}: |grad| {np.abs(g).max():.1e}, max_s |q_k - q_rest| = {d:.3e}, k*dist = {k*d:.3f}, E (no spring) = {energy_full(Q,h):.5f}, lowest Hill eig {lam[0]:.3f}')
elif which == 'v':
    for start in ('rest','harm'):
        Q0 = Qr.copy() if start=='rest' else np.load('min_from_0.0.npy')
        try:
            Q,(a,b),hist = harmonic_section(path, s, Q0, h, kind='vertical', maxit=60, tol=1e-8)
            o = offs(Q); lam = np.linalg.eigvalsh(cyclic_tridiag_matrix(a,b))/h
            marg = np.minimum(Q-QMIN, QMAX-Q)
            print(f'vertical, start {start}: iters {len(hist)}, |grad| {hist[-1][2]:.1e}, E_v {energy_vertical(Q,h):.6f} (start {energy_vertical(Q0,h):.6f}), E_full {energy_full(Q,h):.5f}, offset [{o.min():+.3f},{o.max():+.3f}], lowest eig {lam[0]:.3f}, min margin {marg.min():.4f} joint {np.unravel_index(marg.argmin(),marg.shape)[1]+1}')
        except Exception as ex:
            print('vertical failed', start, ex)
