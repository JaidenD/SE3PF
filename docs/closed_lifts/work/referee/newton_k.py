"""Referee: own Newton solver for the discrete E_k = sum |q_{j+1}-q_j|^2/(2h) + (k/2) h sum |q_j - q_r|^2 over closed
sections (K = I), in the re-centred unit-speed fibre chart, with exact gradient and the cyclic-tridiagonal Hessian
(fibre curvature kappa = -J^+ dJ[x] x). Continuation in k from the rest-posture lift (k = infinity) down to k = 0.
Also: the predicted O(1/k) corrector constant max_s |x.q_f''/(1+(q_f-q_r).kappa)| of Prop. E."""
import sys, time
import numpy as np, rk
from numpy.linalg import norm
m = int(sys.argv[1]) if len(sys.argv) > 1 else 500
rest1000 = np.load("../rest_lift_1000.npy"); Qf = rest1000[::1000//m].copy(); h = rk.L/m; S = np.arange(m)*h

def geom(Q):
    X = np.array([rk.kernel(rk.jb(q)) for q in Q])
    for j in range(1, m):
        if X[j]@X[j-1] < 0: X[j] = -X[j]
    if X[0]@X[-1] < 0: raise RuntimeError("kernel orientation flips around the loop")
    K = np.zeros_like(Q); e = 1e-6
    for j in range(m):
        J = rk.jb(Q[j]); dJ = (rk.jb(Q[j]+e*X[j])-rk.jb(Q[j]-e*X[j]))/(2*e); K[j] = -np.linalg.pinv(J)@(dJ@X[j])
    return X, K

def energy(Q, k): return np.sum((np.roll(Q, -1, 0)-Q)**2)/(2*h) + 0.5*k*h*np.sum((Q-rk.QR)**2)

def grad_hess(Q, k, X, Kp):
    lap = 2*Q-np.roll(Q, 1, 0)-np.roll(Q, -1, 0)
    g = np.sum(lap*X, 1)/h + k*h*np.sum((Q-rk.QR)*X, 1)
    a = (2+np.sum(lap*Kp, 1))/h + k*h*(1+np.sum((Q-rk.QR)*Kp, 1))
    b = -np.sum(X*np.roll(X, -1, 0), 1)/h
    H = np.diag(a)+np.diag(b[:-1], 1)+np.diag(b[:-1], -1); H[0, -1] += b[-1]; H[-1, 0] += b[-1]
    return g, H

def retract(Q, du, X, Kp):
    return np.array([rk.project(Q[j]+du[j]*X[j]+0.5*du[j]**2*Kp[j], S[j], x=X[j]) for j in range(m)])

def solve(Q, k, tol=1e-9, maxit=40):
    for it in range(maxit):
        X, Kp = geom(Q); g, H = grad_hess(Q, k, X, Kp)
        if np.abs(g).max()/h < tol: break
        lam = np.linalg.eigvalsh(H)[0]
        Hs = H if lam > 0 else H+(abs(lam)+1e-3)*np.eye(m)
        du = -np.linalg.solve(Hs, g); E0 = energy(Q, k); t = 1.0
        while True:
            Qn = retract(Q, t*du, X, Kp)
            if energy(Qn, k) <= E0 + 1e-14 or t < 1e-3: break
            t *= 0.5
        Q = Qn
    X, Kp = geom(Q); g, H = grad_hess(Q, k, X, Kp)
    return Q, it, np.abs(g).max()/h, np.linalg.eigvalsh(H)[0]/h

# predicted corrector constant from the rest lift (m = 1000, 4th-order periodic differences)
R = rest1000; hh = rk.L/len(R)
d2 = (-np.roll(R, 2, 0)+16*np.roll(R, 1, 0)-30*R+16*np.roll(R, -1, 0)-np.roll(R, -2, 0))/(12*hh**2)
Xr = np.array([rk.kernel(rk.jb(q)) for q in R]); Kr = []
for j in range(len(R)):
    J = rk.jb(R[j]); e = 1e-6; dJ = (rk.jb(R[j]+e*Xr[j])-rk.jb(R[j]-e*Xr[j]))/(2*e); Kr.append(-np.linalg.pinv(J)@(dJ@Xr[j]))
Kr = np.array(Kr); hf = 1+np.sum((R-rk.QR)*Kr, 1)
corr = np.abs(np.sum(Xr*d2, 1)/hf)
print(f"predicted corrector constant max_s |x.q_f''/h_f| = {corr.max():.3f}; min h_f = {hf.min():.3f}", flush=True)

Q = Qf.copy(); t0 = time.time()
for k in [1e5, 1e4, 1e3, 300, 100, 30, 10, 3, 1, 0.3, 0.1, 0.0]:
    Q, it, gn, l0 = solve(Q, k)
    dist = norm(Q-Qf, axis=1).max()
    print(f"k={k:8.1f}: iters {it:2d}, |g|/h {gn:.1e}, lam0/h {l0:10.4f}, max|q_k-q_f| {dist:.3e}, k*dist {k*dist:8.3f}, "
          f"E_I {np.sum((np.roll(Q,-1,0)-Q)**2)/(2*h):.5f}, margin {np.minimum(Q-rk.QMIN, rk.QMAX-Q).min():.4f}  [{time.time()-t0:.0f}s]", flush=True)
np.save(f"newton_k0_m{m}.npy", Q)
A = np.load("../adv/harm_euclid_500.npy") if m == 500 else None
if A is not None: print("k=0 endpoint vs alternatives' minimizer A (500):", norm(Q-A, axis=1).max())
# pure Newton (k = 0) directly from the rest lift
Q0, it, gn, l0 = solve(Qf.copy(), 0.0)
print(f"Newton at k=0 started from the rest lift: iters {it}, |g|/h {gn:.1e}, E_I {np.sum((np.roll(Q0,-1,0)-Q0)**2)/(2*h):.6f}, "
      f"distance to the continuation endpoint {norm(Q0-Q, axis=1).max():.2e}")
