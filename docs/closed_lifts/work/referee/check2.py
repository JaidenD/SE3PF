"""Referee check of stored discrete lifts with the referee's own kinematics (rk.py):
FK residual, discrete K=I energy, discrete E-L residual x.(q_{j+1}-2q_j+q_{j-1})/h^2, joint margin, sigma_min,
discrete Jacobi operator (re-centred unit chart): lowest eigenvalue / h, # negative, monodromy trace and multipliers."""
import sys, os
import numpy as np
import rk
HL = ".."

def analyse(name, Q):
    m = len(Q); h = rk.L/m; s = np.arange(m)*h
    fkres = max(np.linalg.norm(rk.pose_err(Q[j], s[j])) for j in range(m))
    dQ = np.roll(Q, -1, axis=0) - Q
    E = np.sum(dQ**2)/(2*h)
    X = np.array([rk.kernel(rk.jb(q)) for q in Q])
    for j in range(1, m):
        if X[j] @ X[j-1] < 0: X[j] = -X[j]
    lap = (np.roll(Q, -1, 0) - 2*Q + np.roll(Q, 1, 0))
    EL = np.max(np.abs(np.sum(X*lap, 1)))/h**2
    marg = np.minimum(Q - rk.QMIN, rk.QMAX - Q)
    jm = np.unravel_index(np.argmin(marg), marg.shape)
    sig = min(np.linalg.svd(rk.jb(q), compute_uv=False)[-1] for q in Q)
    kap = np.zeros_like(Q); eps = 1e-6
    for j in range(m):
        J = rk.jb(Q[j]); x = X[j]
        dJ = (rk.jb(Q[j]+eps*x) - rk.jb(Q[j]-eps*x))/(2*eps)
        kap[j] = -np.linalg.pinv(J) @ (dJ @ x)
    a = (2 + np.sum((2*Q - np.roll(Q, 1, 0) - np.roll(Q, -1, 0))*kap, 1))/h
    b = -np.sum(X*np.roll(X, -1, 0), 1)/h
    H = np.diag(a) + np.diag(b[:-1], 1) + np.diag(b[:-1], -1)
    H[0, -1] += b[-1]; H[-1, 0] += b[-1]
    ev = np.linalg.eigvalsh(H)
    Mon = np.eye(2)
    for j in range(m):
        T = np.array([[-a[j]/b[j], -b[j-1]/b[j]], [1, 0]])
        Mon = T @ Mon
    mu = np.linalg.eigvals(Mon)
    # gradient of discrete energy in the unit chart (should be ~0 at a critical point): g_j = (2q_j-q_{j-1}-q_{j+1}).x_j/h
    g = np.sum((2*Q - np.roll(Q, 1, 0) - np.roll(Q, -1, 0))*X, 1)/h
    print(f"{name:30s} m={m:4d} FKres={fkres:.1e} E_I={E:.6f} |g|/h={np.abs(g).max()/h:.1e} EL={EL:.2e} "
          f"margin={marg[jm]:.4f}(j{jm[1]+1}) sigmin={sig:.4f} lam0/h={ev[0]/h:.4f} #neg={np.sum(ev<0)} "
          f"trM={np.trace(Mon):.4f} det={np.linalg.det(Mon):.6f} mu={np.round(np.sort(mu.real),4)}", flush=True)
    return Q

items = [("rest lift (lead, 1000)", np.load(f"{HL}/rest_lift_1000.npy")),
         ("A harmonic (lead, 1000)", np.load(f"{HL}/harmonic_full_1000.npy")),
         ("A harmonic (num, 1000)", np.load(f"{HL}/num/data/harm_I_1000_u0.npz")["q"]),
         ("B far (num, start-2, 1000)", np.load(f"{HL}/num/data/far_I_start-2.npz")["q"]),
         ("C far (num, start+2, 1000)", np.load(f"{HL}/num/data/far_I_start+2.npz")["q"]),
         ("B adv start-2 (500)", np.load(f"{HL}/adv/harm_euclid_500_start-2.npy")),
         ("C adv start+2 (500)", np.load(f"{HL}/adv/harm_euclid_500_start+2.npy")),
         ("lead min_from_2.0 (250)", np.load(f"{HL}/min_from_2.0.npy")),
         ("lead vertical (250)", np.load(f"{HL}/vertical_250.npy")),
         ("midsection adv (500)", np.load(f"{HL}/adv/midpoint_500.npz")["Qm"]),
         ("K=M near-rest (num, 1000)", np.load(f"{HL}/num/data/harm_M_1000_u0.npz")["q"]),
         ("K=M far start-2 (num, 1000)", np.load(f"{HL}/num/data/far_M_start-2.npz")["q"])]
for nm, Q in items:
    analyse(nm, Q)
