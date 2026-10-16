"""Referee check of the stored discrete lifts with independent kinematics:
FK residual, discrete energy (K=I), E-L residual, joint margins, sigma_min, discrete Jacobi (Hill) operator,
its lowest eigenvalue (K-normalised) and the discrete monodromy (Floquet multipliers)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import rk

HL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print("L =", rk.L)

def analyse(name, Q):
    m = len(Q); h = rk.L/m; s = np.arange(m)*h
    fkres = max(np.linalg.norm(rk.pose_err(Q[j], s[j])) for j in range(m))
    dQ = np.roll(Q, -1, axis=0) - Q
    E = np.sum(dQ**2)/(2*h)
    X = np.array([rk.kernel(rk.jb(q)) for q in Q])
    # orient consistently with first
    for j in range(1, m):
        if X[j] @ X[j-1] < 0: X[j] = -X[j]
    lap = (np.roll(Q, -1, 0) - 2*Q + np.roll(Q, 1, 0))
    EL = np.max(np.abs(np.sum(X*lap, 1)))/h**2
    marg = np.minimum(Q - rk.QMIN, rk.QMAX - Q)
    jm = np.unravel_index(np.argmin(marg), marg.shape)
    sig = min(np.linalg.svd(rk.jb(q), compute_uv=False)[-1] for q in Q)
    # curvature of self-motion curve: kappa = -J^+ (dJ[x] x)
    kap = np.zeros_like(Q); eps = 1e-6
    for j in range(m):
        J = rk.jb(Q[j]); x = X[j]
        dJ = (rk.jb(Q[j]+eps*x) - rk.jb(Q[j]-eps*x))/(2*eps)
        kap[j] = -np.linalg.pinv(J) @ (dJ @ x)
    a = (2 + np.sum((2*Q - np.roll(Q, 1, 0) - np.roll(Q, -1, 0))*kap, 1))/h
    b = -np.sum(X*np.roll(X, -1, 0), 1)/h    # coupling j, j+1
    H = np.diag(a) + np.diag(b[:-1], 1) + np.diag(b[:-1], -1)
    H[0, -1] = H[-1, 0] = b[-1]
    ev = np.linalg.eigvalsh(H)
    # monodromy of b_{j-1} phi_{j-1} + a_j phi_j + b_j phi_{j+1} = 0
    Mon = np.eye(2)
    for j in range(m):
        T = np.array([[-a[j]/b[j], -b[j-1]/b[j]], [1, 0]])
        Mon = T @ Mon
    mu = np.linalg.eigvals(Mon)
    print(f"{name:28s} m={m} FKres={fkres:.1e} E={E:.6f} EL={EL:.2e} margin={marg[jm]:.4f} (joint {jm[1]+1}) "
          f"sigmin={sig:.4f} lamH/h={ev[0]/h:.4f} (#neg {np.sum(ev<0)}) tr={np.trace(Mon):.4f} det={np.linalg.det(Mon):.6f} mu={np.sort(mu.real)}")
    return dict(E=E, X=X, h=h)

lifts = [("rest_lift_1000 (lead)", "rest_lift_1000.npy"), ("A harmonic_full_1000 (lead)", "harmonic_full_1000.npy"),
         ("A adv 500", "adv/harm_euclid_500.npy"), ("B adv 500 start-2", "adv/harm_euclid_500_start-2.npy"),
         ("C adv 500 start+2", "adv/harm_euclid_500_start+2.npy"), ("lead min_from_2.0 (250)", "min_from_2.0.npy"),
         ("lead vertical_250", "vertical_250.npy")]
for nm, f in lifts:
    analyse(nm, np.load(os.path.join(HL, f)))
d = np.load(os.path.join(HL, "adv/midpoint_500.npz"))
analyse("midsection adv 500", d["Qm"])
for nm in ["far_I_start-2", "far_I_start+2"]:
    z = np.load(os.path.join(HL, "num/data", nm+".npz"))
    print(nm, z.files)
