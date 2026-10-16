import numpy as np, time, sys
from pandageo import *
from lifts import *
path = ConingCircle(); L = path.L
m = int(sys.argv[1]) if len(sys.argv)>1 else 500
s, Qr, gap = rest_lift(path, m)
h = L/m
print('m', m, 'rest lift energy E =', energy_full(Qr, h))
t0=time.time()
Qh, (a,b), hist = harmonic_section(path, s, Qr.copy(), h, kind='full', verbose=True)
print('time', time.time()-t0)
x, kap, J = fibre_jets(path, s, Qh)
e,_ = pose_error(Qh, path, s)
print('on fibres: max |e| =', np.abs(e).max())
print('harmonic lift energy E =', energy_full(Qh, h), ' rest-lift energy', energy_full(Qr, h))
d = Qh - Qr
print('max_s |q_h - q_rest| (2-norm) =', np.linalg.norm(d,axis=1).max(), ' at s =', s[np.linalg.norm(d,axis=1).argmax()])
xr = kernel_vector(fk(Qr)[2])
print('beta-offset of harmonic lift from rest lift, x_rest^T(q_h - q_rest): range', (np.einsum('ma,ma->m',xr,d)).min(), (np.einsum('ma,ma->m',xr,d)).max())
marg = np.minimum(Qh-QMIN, QMAX-Qh)
print('harmonic lift: min joint margin', marg.min(), 'joint', np.unravel_index(marg.argmin(), marg.shape)[1]+1)
sv = np.linalg.svd(J, compute_uv=False)[:,-1]
svr = np.linalg.svd(fk(Qr)[2], compute_uv=False)[:,-1]
print('min sigma_min(J): harmonic', sv.min(), ' rest', svr.min())
# U = zeta^T (J J^T)^{-1} zeta : horizontal (task-compatibility) term; vertical term
zeta = path.velocity(s)
def split(Q):
    _,_,J = fk(Q); x = kernel_vector(J)
    U = np.einsum('ma,mab,mb->m', zeta, np.linalg.inv(J@np.swapaxes(J,1,2)), zeta)
    dQ = (np.roll(Q,-1,0)-np.roll(Q,1,0))/(2*h)
    v = np.einsum('ma,ma->m', x, dQ)
    return 0.5*h*U.sum(), 0.5*h*(v**2).sum(), U, v
for name,Q in (('rest',Qr),('harmonic',Qh)):
    Eh, Ev, U, v = split(Q)
    print(f'{name:9s}: horizontal part 1/2 int U = {Eh:.5f}, vertical part 1/2 int v^2 = {Ev:.5f}, sum {Eh+Ev:.5f};  max|q\'| = {np.linalg.norm((np.roll(Q,-1,0)-np.roll(Q,1,0))/(2*h),axis=1).max():.3f}')
H = cyclic_tridiag_matrix(a,b)
lam = np.linalg.eigvalsh(H)
print('discrete Hessian (periodic Jacobi/Hill operator) lowest eigenvalues *h:', lam[:4]*h)
print('twist condition: max b_j*h =', (b*h).max(), '(must be < 0)')
tr, det = monodromy_trace(a, b)
print('discrete monodromy: trace', tr, 'det', det, ' -> multipliers', np.roots([1,-tr,det]))
np.save(f'harmonic_full_{m}.npy', Qh); np.save(f'rest_{m}.npy', Qr)
