"""Ideal passive dragged arm (exact holonomic constraint), in arc-length time s (speed nu = 1; the equations are
homogeneous in nu without gravity and damping, so q(nu t) solves them for every nu):
    K q'' + C(q, q') q' = J_b^T lambda,   FK(q(s)) = gamma(s)
(K = M from MuJoCo with C from RNE, or K = I with C = 0). Index-reduced with Baumgarte stabilization:
    J q'' = zeta'(s) - J' q' + 2 a (zeta - J q') + a^2 e,   e = Log_I(FK(q)^-1 gamma(s)),  a = 20 per unit s,
RK4 with n steps per lap. Coordinates on the 2-dim constraint manifold at s = 0 (= L): u = x_0^T (q - q*_rest(0))
(the chart of task 1) and w = x(q)^T q' (x(q) the unit kernel vector). Shooting (Newton on the Poincare map
(u, w) -> (u, w)) finds the periodic orbit; the 2x2 monodromy is the Jacobian of that map (central differences).
"""
import json, os, sys
import numpy as np
from pandakin import fk, kernel_vector, path_pose, path_vel, path_acc, pose_error, L, MjDyn, DATA, constrained_ik, \
    load_rest
from hl import spectral_deriv, trig_interp

dyn = MjDyn()
A_B = 20.0
_, qrest_all, x_all = load_rest(1000)
qrest0, x0 = qrest_all[0], x_all[0]
Rt0, pt0 = path_pose(0.0)


def accel(q, v, s, met):
    R, p, J = fk(q)
    Rt, pt = path_pose(s)
    e = pose_error(R, p, Rt, pt)
    eps = 1e-6
    nv = np.linalg.norm(v)
    if nv > 0:
        _, _, Jp = fk(q + eps * v / nv)
        _, _, Jm = fk(q - eps * v / nv)
        Jdv = (Jp - Jm) @ v / (2 * eps) * nv
    else:
        Jdv = np.zeros(6)
    z = path_vel(s)
    rhs_c = path_acc(s) - Jdv + 2 * A_B * (z - J @ v) + A_B ** 2 * e
    if met == "I":
        K = np.eye(7); f = np.zeros(7)
    else:
        K = dyn.M(q); f = -dyn.Cv(q, v)
    KKT = np.zeros((13, 13))
    KKT[:7, :7] = K; KKT[:7, 7:] = -J.T; KKT[7:, :7] = J
    sol = np.linalg.solve(KKT, np.r_[f, rhs_c])
    return sol[:7]


def lap(q, v, met, n=2000, s0=0.0, store=False, laps=1):
    h = L / n
    tr = [(q.copy(), v.copy())] if store else None
    for j in range(n * laps):
        s = s0 + j * h
        a1 = accel(q, v, s, met)
        q2, v2 = q + 0.5 * h * v, v + 0.5 * h * a1
        a2 = accel(q2, v2, s + 0.5 * h, met)
        q3, v3 = q + 0.5 * h * v2, v + 0.5 * h * a2
        a3 = accel(q3, v3, s + 0.5 * h, met)
        q4, v4 = q + h * v3, v + h * a3
        a4 = accel(q4, v4, s + h, met)
        q = q + (h / 6) * (v + 2 * v2 + 2 * v3 + v4)
        v = v + (h / 6) * (a1 + 2 * a2 + 2 * a3 + a4)
        if store:
            tr.append((q.copy(), v.copy()))
    return q, v, tr


def state_from_uw(u, w, guess):
    q, ok = constrained_ik(Rt0, pt0, x0, x0 @ qrest0 + u, guess.copy(), tol=1e-14)
    assert ok
    _, _, J = fk(q)
    v = J.T @ np.linalg.solve(J @ J.T, path_vel(0.0)) + w * kernel_vector(J)
    return q, v


def uw_from_state(q, v):
    _, _, J = fk(q)
    return np.array([x0 @ (q - qrest0), kernel_vector(J) @ v])


def poincare(uw, met, n, guess):
    q, v = state_from_uw(uw[0], uw[1], guess)
    qe, ve, _ = lap(q, v, met, n)
    return uw_from_state(qe, ve)


def monodromy(uw, met, n, guess, d=(1e-6, 1e-5)):
    Mo = np.zeros((2, 2))
    for j in range(2):
        e = np.zeros(2); e[j] = d[j]
        Mo[:, j] = (poincare(uw + e, met, n, guess) - poincare(uw - e, met, n, guess)) / (2 * d[j])
    return Mo


if __name__ == "__main__":
    res = {}
    for met in ("I", "M"):
        H = np.load(os.path.join(DATA, "harm_%s_1000_u0.npz" % met))
        qh = H["q"]
        dqh = spectral_deriv(qh, L, 1)
        uw_h = uw_from_state(qh[0], dqh[0])
        out = dict(uw_harmonic_discrete=uw_h.tolist())
        for n in (1000, 2000):
            P = poincare(uw_h, met, n, qh[0])
            out["return_map_defect_n%d" % n] = (P - uw_h).tolist()
        # shooting for the exact periodic orbit
        n = 2000
        uw = uw_h.copy()
        for it in range(6):
            F = poincare(uw, met, n, qh[0]) - uw
            Mo = monodromy(uw, met, n, qh[0])
            print(met, "shooting it %d |P(uw) - uw| = %.3e" % (it, np.abs(F).max()), flush=True)
            if np.abs(F).max() < 1e-11:
                break
            uw = uw - np.linalg.solve(Mo - np.eye(2), F)
        out["uw_periodic_orbit"] = uw.tolist()
        out["shooting_residual"] = float(np.abs(F).max())
        out["monodromy"] = Mo.tolist()
        ev = np.linalg.eigvals(Mo)
        out["floquet_multipliers"] = [str(e) for e in ev]
        out["trace"] = float(np.trace(Mo)); out["det"] = float(np.linalg.det(Mo))
        out["greene_residue"] = float((2 - np.trace(Mo)) / 4)
        # the periodic orbit vs the discrete harmonic lift (N = 1000), max joint distance over the lap
        q, v = state_from_uw(uw[0], uw[1], qh[0])
        _, _, tr = lap(q, v, met, 4000, store=True)
        Q = np.array([a for a, b in tr]); V = np.array([b for a, b in tr])
        sg = np.arange(4001) * L / 4000
        qh_g = trig_interp(qh, L, sg)
        out["periodic_orbit_vs_discrete_harmonic_N1000_max_joint2norm"] = float(np.linalg.norm(Q - qh_g, axis=1).max())
        for N in (250, 500):
            Hn = np.load(os.path.join(DATA, "harm_%s_%d_u0.npz" % (met, N)))
            out["periodic_orbit_vs_discrete_harmonic_N%d" % N] = float(
                np.linalg.norm(Q - trig_interp(Hn["q"], L, sg), axis=1).max())
        out["closure_after_one_lap_n4000_q"] = float(np.linalg.norm(Q[-1] - Q[0]))
        out["closure_after_one_lap_n4000_v"] = float(np.linalg.norm(V[-1] - V[0]))
        np.savez(os.path.join(DATA, "dae_orbit_%s.npz" % met), Q=Q, V=V, sg=sg, uw=uw, Mo=Mo)
        res[met] = out
        print(met, json.dumps(out, indent=1), flush=True)
        json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_dae.json"), "w"), indent=1)
