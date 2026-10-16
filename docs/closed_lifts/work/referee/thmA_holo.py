"""Referee: (1) Theorem A on the paper's path, independently: transport, holonomy angles, frame energy vs sum theta^2/L,
closed perturbations (including one in the other homotopy class and the non-principal logarithm);
(2) pseudoinverse holonomy from the rest-lift node q*(0) (one lap of q' = J^+ (zeta + 50 e))."""
import numpy as np, rk
from numpy.linalg import norm
W = np.diag([rk.WR]*3+[rk.WP]*3); L = rk.L
def dzeta(s, e=1e-5): return (rk.zeta(s+e)-rk.zeta(s-e))/(2*e)
def rhs(s, N): z = rk.zeta(s); return -np.outer(z, dzeta(s)@W@N)
# initial frame: W-Gram-Schmidt on zeta(0), e1..e6
vs = [rk.zeta(0.0)] + list(np.eye(6)); B = []
for v in vs:
    w = v - sum((b@W@v)*b for b in B)
    if np.sqrt(w@W@w) > 1e-8: B.append(w/np.sqrt(w@W@w))
N0 = np.array(B[1:6]).T
n = 2000; hs = L/n; Ns = [N0]; N = N0.copy()
for k in range(n):
    s = k*hs
    k1 = rhs(s, N); k2 = rhs(s+hs/2, N+hs/2*k1); k3 = rhs(s+hs/2, N+hs/2*k2); k4 = rhs(s+hs, N+hs*k3)
    N = N + hs/6*(k1+2*k2+2*k3+k4); Ns.append(N)
Ns = np.array(Ns)
H = N0.T@W@Ns[-1]; U, S, Vt = np.linalg.svd(H); H = U@Vt
ev = np.linalg.eigvals(H); th = np.sort(np.abs(np.angle(ev)))[::-1]
print("holonomy angles:", th[[0, 2]], " det H =", np.linalg.det(H))
from scipy_free_logm import logm_so
from scipy_free_logm import expm
Phi = logm_so(H); print("|expm(Phi)-H| =", norm(expm(Phi)-H), " |Phi|_F^2/(2L) =", norm(Phi)**2/(2*L), " sum th^2/L =", (th[0]**2+th[2]**2)/L)
sg = np.arange(n+1)*hs
def frame_energy(Qfun, nq=4000):
    """E = 1/2 oint |N^T W N'|^2 for N = N0 Q(s), via |Q'|_F (Q from a callable), midpoint quadrature"""
    ss = (np.arange(nq)+0.5)*L/nq; e = 1e-6; tot = 0
    for s in ss:
        dQ = (Qfun(s+e)-Qfun(s-e))/(2*e); tot += 0.5*norm(dQ)**2
    return tot*L/nq
Qp = lambda s: expm(-(s/L)*Phi)
print("paper frame energy (quadrature):", frame_energy(Qp))
rng = np.random.default_rng(3)
for k in (1, 2, 3):
    Xi = rng.standard_normal((5, 5)); Xi = Xi-Xi.T
    for eps in (0.01, 0.1):
        Qe = lambda s: expm(-(s/L)*Phi)@expm(eps*np.sin(2*np.pi*k*s/L)*Xi)
        print(f"  perturbation k={k} eps={eps}: E = {frame_energy(Qe, 2000):.6f}  (pred +eps^2 pi^2 k^2 |Xi|^2/L = {eps**2*np.pi**2*k**2*norm(Xi)**2/L:.6f})")
# a closed frame from a non-principal logarithm: rotate the small plane by theta2 - 2 pi
w, V = np.linalg.eig(Phi)
i2 = np.argmin(np.abs(w.imag-th[2])) ; i2 = [i for i in range(5) if abs(abs(w[i].imag)-th[2]) < 1e-9]
Pl = np.real(sum(V[:, i][:, None]@np.linalg.inv(V)[i][None, :]*w[i] for i in i2))  # component of Phi in plane 2
Phi2 = Phi - Pl + Pl*(th[2]-2*np.pi)/th[2]
print("  non-principal log: |expm(Phi2)-H| =", norm(expm(Phi2)-H), " E =", frame_energy(lambda s: expm(-(s/L)*Phi2), 2000),
      " pred", norm(Phi2)**2/(2*L))
# (2) pseudoinverse holonomy from the rest lift
rest = np.load("../rest_lift_1000.npy"); q = rest[0].copy()
def f(s, q):
    J = rk.jb(q); e = -rk.pose_err(q, s); return np.linalg.pinv(J)@(rk.zeta(s)+50*e)
nn = 4000; hq = L/nn
for k in range(nn):
    s = k*hq; k1 = f(s, q); k2 = f(s+hq/2, q+hq/2*k1); k3 = f(s+hq/2, q+hq/2*k2); k4 = f(s+hq, q+hq*k3)
    q = q + hq/6*(k1+2*k2+2*k3+k4)
x = rk.kernel(rk.jb(rest[0]))
print("pseudoinverse holonomy from q*_rest(0): |q(L)-q(0)| =", norm(q-rest[0]), " x.(q(L)-q(0)) =", x@(q-rest[0]),
      " pose residual", norm(rk.pose_err(q, 0.0)))
