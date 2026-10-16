"""Task 2: holonomy of the self-motion bundle for the two connections, and 'horizontal lift + uniform drift'.

Horizontal lift for the metric K (connection = K-orthogonal complement of ker J_b):
    q' = J^#(q) (zeta_gamma(s) + ke e(q, s)),   J^# = K^-1 J^T (J K^-1 J^T)^-1,
e = Log_I(FK(q)^-1 gamma(s)) a pose-error feedback (ke = 50 per unit s) that keeps FK(q) on gamma against
integration error; it is horizontal too, so it adds no vertical drift. RK4 with n steps per lap.
Uniform drift: q' = J^#(...) + c nK(q), nK the K-unit kernel vector, c fixed by shooting so that q(L) = q(0).
"""
import json, os
import numpy as np
from pandakin import fk, kernel_vector, path_pose, path_vel, pose_error, L, MjDyn, DATA, constrained_ik, load_rest
from hl import spectral_deriv, trig_interp, Problem

dyn = MjDyn()
KE = 50.0


def rhs(q, s, met, c=0.0):
    R, p, J = fk(q)
    Rt, pt = path_pose(s)
    e = pose_error(R, p, Rt, pt)
    z = path_vel(s) + KE * e
    x = kernel_vector(J)
    if met == "I":
        dq = J.T @ np.linalg.solve(J @ J.T, z)
        nK = x
    else:
        M = dyn.M(q)
        MiJt = np.linalg.solve(M, J.T)
        dq = MiJt @ np.linalg.solve(J @ MiJt, z)
        nK = x / np.sqrt(x @ M @ x)
    return dq + c * nK


def integrate(q0, met, c=0.0, n=4000, laps=1, store=False):
    h = L / n
    q = q0.copy()
    traj = [q.copy()]
    E = 0.0
    for lap in range(laps):
        for i in range(n):
            s = i * h
            k1 = rhs(q, s, met, c)
            k2 = rhs(q + 0.5 * h * k1, s + 0.5 * h, met, c)
            k3 = rhs(q + 0.5 * h * k2, s + 0.5 * h, met, c)
            k4 = rhs(q + h * k3, s + h, met, c)
            q = q + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            if store:
                traj.append(q.copy())
    return q, np.array(traj)


def energy_of_traj(traj, met):
    """Discrete energy (midpoint rule) of a stored closed trajectory (n + 1 points, last = first approx.)."""
    q = traj[:-1]
    n = len(q)
    h = L / n
    d = np.roll(q, -1, 0) - q
    if met == "I":
        return 0.5 * np.sum(d * d) / h
    m = 0.5 * (q + np.roll(q, -1, 0))
    M = dyn.Mbatch(m)
    return 0.5 * np.einsum("ki,kij,kj->", d, M, d) / h


s0, qrest0, x0 = [a[0] for a in load_rest(1000)]
Rt0, pt0 = path_pose(0.0)


def u0(q):
    return float(x0 @ (q - qrest0))


def fibre_length(ua, ub, met, m=24):
    """K-length of the fibre over gamma(0) between chart values ua and ub (Gauss-Legendre in u)."""
    xg, wg = np.polynomial.legendre.leggauss(m)
    us = 0.5 * (ub - ua) * xg + 0.5 * (ua + ub)
    tot = 0.0
    q = qrest0.copy()
    for uu, ww in zip(us, wg):
        q, ok = constrained_ik(Rt0, pt0, x0, x0 @ qrest0 + uu, q)
        _, _, J = fk(q)
        t = np.linalg.solve(np.vstack([J, x0]), np.r_[np.zeros(6), 1.0])
        K = np.eye(7) if met == "I" else dyn.M(q)
        tot += ww * np.sqrt(t @ K @ t)
    return 0.5 * (ub - ua) * tot


res = {}
for met in ("I", "M"):
    H = np.load(os.path.join(DATA, "harm_%s_1000_u0.npz" % met))
    qh, uh = H["q"], H["u"]
    qa = qh[0]
    out = {}
    for n in (2000, 4000):
        qL, _ = integrate(qa, met, n=n)
        out["n%d" % n] = dict(gap_joint2norm=float(np.linalg.norm(qL - qa)), gap_u=u0(qL) - u0(qa))
    qL, traj_h = integrate(qa, met, n=4000, store=True)
    R, p, _ = fk(qL)
    out["constraint_residual_end"] = float(np.abs(pose_error(R, p, Rt0, pt0)).max())
    tau_u = u0(qL) - u0(qa)
    out["tau_u"] = tau_u
    out["tau_joint2norm"] = float(np.linalg.norm(qL - qa))
    out["tau_Klength"] = float(fibre_length(u0(qa), u0(qL), met))
    # holonomy as a map of the fibre: displacement for several start points on the fibre over gamma(0)
    disp = []
    for ust in (-0.6, -0.3, 0.0, 0.3, 0.6, 0.9):
        qs, ok = constrained_ik(Rt0, pt0, x0, x0 @ qrest0 + ust, qrest0.copy())
        if not ok:
            continue
        qe, _ = integrate(qs, met, n=2000)
        disp.append((ust, u0(qe) - ust, float(fibre_length(ust, u0(qe), met))))
    out["holonomy_map_u_start_du_dKlength"] = disp
    # multi-lap drift from the rest-posture lift (the 0.24 rad/lap of the plain pseudoinverse)
    q = qrest0.copy()
    laps = []
    for lap in range(3):
        qn, _ = integrate(q, met, n=2000)
        laps.append(float(np.linalg.norm(qn - q)))
        q = qn
    out["rest_start_gap_per_lap_joint2norm"] = laps

    # uniform vertical drift closing from q_harm(0): shoot on c (secant)
    def gap(c):
        qe, _ = integrate(qa, met, c=c, n=2000)
        return u0(qe) - u0(qa)
    c0 = -out["tau_Klength"] / L
    c1 = c0 * 1.05
    g0, g1 = gap(c0), gap(c1)
    for it in range(8):
        c2 = c1 - g1 * (c1 - c0) / (g1 - g0)
        c0, g0 = c1, g1
        c1, g1 = c2, gap(c2)
        if abs(g1) < 1e-12:
            break
    cstar = c1
    qe, traj_ud = integrate(qa, met, c=cstar, n=4000, store=True)
    out["uniform_drift_c"] = cstar
    out["uniform_drift_closure_joint"] = float(np.linalg.norm(qe - qa))
    out["minus_tauK_over_L"] = -out["tau_Klength"] / L
    # compare with the harmonic lift on the same s grid (trig interpolation of the N = 1000 nodes)
    sg = np.arange(4001) * L / 4000
    qh_g = trig_interp(qh, L, sg)
    out["uniform_drift_vs_harmonic_max_joint2norm"] = float(np.linalg.norm(traj_ud - qh_g, axis=1).max())
    out["horizontal_vs_harmonic_max_joint2norm"] = float(np.linalg.norm(traj_h - qh_g, axis=1).max())
    out["E_uniform_drift"] = float(energy_of_traj(traj_ud, met))
    out["E_horizontal_open"] = float(energy_of_traj(traj_h, met))  # not closed: for reference only
    # harmonic lift energy on the same 4000 grid (re-evaluated from the interpolant)
    out["E_harmonic_on_4000"] = float(energy_of_traj(qh_g, met))
    # vertical velocity of the harmonic lift: v = K(nK, q'), and horizontal energy density
    dq = spectral_deriv(qh, L, 1)
    _, _, J = fk(qh)
    x = kernel_vector(J)
    if met == "I":
        v = np.einsum("ki,ki->k", x, dq)
    else:
        M = dyn.Mbatch(qh)
        v = np.einsum("ki,kij,kj->k", x, M, dq) / np.sqrt(np.einsum("ki,kij,kj->k", x, M, x))
    out["v_min"], out["v_max"], out["v_mean"] = float(v.min()), float(v.max()), float(v.mean())
    out["v_integral"] = float(v.mean() * L)
    out["v_nonuniformity_(max-min)/|mean|"] = float((v.max() - v.min()) / abs(v.mean()))
    res[met] = out
    np.savez(os.path.join(DATA, "holonomy_%s.npz" % met), traj_h=traj_h, traj_ud=traj_ud, qh_g=qh_g, sg=sg, v=v,
             s=H["s"])
    print(met, json.dumps(out, indent=1), flush=True)
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_holonomy.json"), "w"), indent=1)
