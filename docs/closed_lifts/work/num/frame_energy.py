"""Task 3: the paper's normal frame N(s) minimizes E[N] = 1/2 oint |N^T W N'|_F^2 ds over closed W-orthonormal
normal frames (Theorem A of the working hypothesis); minimum value sum theta_i^2 / L.

Data: N(s_k), N'(s_k) from Tube::frame at s_k = k L / N (N = 2000 and coarser), the tube's holonomy angles.
Perturbed closed frames: N_eps = N exp(eps f(s) Psi) (f(0) = f(L) = 0) and two-factor non-commuting versions;
their derivative is formed exactly from N, N' and the closed-form derivative of the exponential factors.
Energies by the periodic trapezoid rule (spectrally accurate for smooth periodic integrands; the f = (s/L)(1 - s/L)
integrand has a kink at s = 0, so second order there).
"""
import json, os
import numpy as np
from pandakin import load_frame, W6, META, L, path_vel


def expm_skew(A):
    """exp of a (batch of) skew 5x5 matrices via eigendecomposition of the Hermitian iA."""
    w, V = np.linalg.eigh(1j * A)
    return np.real(np.einsum("...ij,...j,...kj->...ik", V, np.exp(-1j * w), V.conj()))


def dexp_factor(f, df, Psi):
    """E = exp(f Psi), E' = f' Psi E (Psi commutes with exp(f Psi))."""
    E = expm_skew(f[:, None, None] * Psi[None])
    return E, df[:, None, None] * (Psi[None] @ E)


def energy(Nf, dN, h):
    A = np.swapaxes(Nf, 1, 2) @ W6 @ dN
    return 0.5 * h * np.sum(A * A), A


th1, th2 = META["theta1"], META["theta2"]
Emin = (th1 ** 2 + th2 ** 2) / L
res = dict(theta1=th1, theta2=th2, L=L, predicted_min=Emin)
for N in (250, 500, 1000, 2000):
    s, Nf, dN = load_frame(N)
    h = L / N
    E0, A = energy(Nf, dN, h)
    Phi = -L * A.mean(0)
    res["N%d" % N] = dict(E=E0, rel_err=(E0 - Emin) / Emin,
                          orth_err=float(np.abs(np.swapaxes(Nf, 1, 2) @ W6 @ Nf - np.eye(5)).max()),
                          normal_err=float(np.abs(np.einsum("kia,ij,kj->ka", Nf, W6, path_vel(s))).max()),
                          NtWdN_const_err=float(np.abs(A - A.mean(0)).max()),
                          Phi_skew_err=float(np.abs(Phi + Phi.T).max()))
    print("N=%d: E[N] = %.12f, predicted sum theta^2/L = %.12f, rel diff %.2e; N^T W N' constant to %.1e" % (
        N, E0, Emin, (E0 - Emin) / Emin, res["N%d" % N]["NtWdN_const_err"]))

# perturbations on the finest grid
s, Nf, dN = load_frame(2000)
h = L / 2000
E0, A = energy(Nf, dN, h)
Phi = -L * A.mean(0)
Phi = 0.5 * (Phi - Phi.T)
wPhi = np.sort(np.abs(np.linalg.eigvals(Phi).imag))[::-1]
res["Phi_eigen_angles_from_frame"] = wPhi[::2].tolist()
rng = np.random.default_rng(11)
pert = []
fam = {
    "sin(2 pi s/L)": (lambda s: np.sin(2 * np.pi * s / L), lambda s: 2 * np.pi / L * np.cos(2 * np.pi * s / L),
                      lambda P: np.pi ** 2 / L * np.sum(P * P)),
    "(s/L)(1-s/L)": (lambda s: (s / L) * (1 - s / L), lambda s: (1 - 2 * s / L) / L,
                     lambda P: 0.5 * np.sum(P * P) / (3 * L)),
}
for trial in range(3):
    B = rng.standard_normal((5, 5))
    Psi = B - B.T
    for name, (f, df, pred) in fam.items():
        rows = []
        for eps in (1e-3, 1e-2, 1e-1, 0.3):
            E, dE = dexp_factor(eps * f(s), eps * df(s), Psi)
            N2 = Nf @ E
            dN2 = dN @ E + Nf @ dE
            En, _ = energy(N2, dN2, h)
            rows.append((eps, En - E0, (En - E0) / eps ** 2))
        pert.append(dict(trial=trial, f=name, predicted_coeff=pred(Psi), rows=rows))
        print("random Psi %d, f = %-14s: dE/eps^2 =" % (trial, name), ["%.6g" % r[2] for r in rows],
              " predicted %.6g; dE > 0 for all eps: %s" % (pred(Psi), all(r[1] > 0 for r in rows)))
# nonuniform rate along Phi itself: N exp(-eps sin(2 pi s / L)/(2 pi) Phi) = N0 exp(-rho(s) Phi)
rows = []
for eps in (1e-3, 1e-2, 1e-1, 0.5):
    f = -eps * np.sin(2 * np.pi * s / L) / (2 * np.pi)
    df = -eps * np.cos(2 * np.pi * s / L) / L
    E, dE = dexp_factor(f, df, Phi)
    En, _ = energy(Nf @ E, dN @ E + Nf @ dE, h)
    rows.append((eps, En - E0, (En - E0) / eps ** 2))
pert.append(dict(f="nonuniform rate rho(s) = s/L + eps sin(2 pi s/L)/(2 pi)", predicted_coeff=np.sum(Phi * Phi) / (4 * L),
                 rows=rows))
print("nonuniform rate along Phi: dE/eps^2 =", ["%.6g" % r[2] for r in rows], " predicted |Phi|^2/(4L) = %.6g" % (
    np.sum(Phi * Phi) / (4 * L)))
# non-commuting two-factor perturbations N exp(eps f1 Psi1) exp(eps f2 Psi2)
for trial in range(3):
    B1, B2 = rng.standard_normal((5, 5)), rng.standard_normal((5, 5))
    P1, P2 = B1 - B1.T, B2 - B2.T
    f1, df1 = np.sin(2 * np.pi * s / L), 2 * np.pi / L * np.cos(2 * np.pi * s / L)
    f2, df2 = 1 - np.cos(4 * np.pi * s / L), 4 * np.pi / L * np.sin(4 * np.pi * s / L)
    rows = []
    for eps in (1e-3, 1e-2, 1e-1, 0.3):
        E1, dE1 = dexp_factor(eps * f1, eps * df1, P1)
        E2, dE2 = dexp_factor(eps * f2, eps * df2, P2)
        N2 = Nf @ E1 @ E2
        dN2 = dN @ E1 @ E2 + Nf @ dE1 @ E2 + Nf @ E1 @ dE2
        En, _ = energy(N2, dN2, h)
        rows.append((eps, En - E0, (En - E0) / eps ** 2))
    pert.append(dict(trial=trial, f="two-factor non-commuting", rows=rows))
    print("two-factor %d: dE/eps^2 =" % trial, ["%.6g" % r[2] for r in rows], " all dE > 0: %s" % all(r[1] > 0 for r in rows))
# the other homotopy class: rotate the first plane by theta1 - 2 pi instead of theta1
res["other_class_energy"] = ((2 * np.pi - th1) ** 2 + th2 ** 2) / L
res["perturbations"] = pert
print("E_min = %.10f; minimizer in the other homotopy class of closed frames: %.6f" % (Emin, res["other_class_energy"]))
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_frame.json"), "w"), indent=1,
          default=float)
