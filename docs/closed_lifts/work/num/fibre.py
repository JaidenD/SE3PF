"""Topology of the self-motion fibre over gamma(s0): trace the connected component of Pi^-1(gamma(s0)) through the
rest-posture lift node in the torus T^7 (joint limits ignored), by RK4 on q' = x(q) (unit kernel vector, Euclidean
arc length tau) with a Newton projection back onto FK = gamma(s0) after every step. Reports the fibre length, where
it leaves the joint limits, sigma_min(J_b) along it, and where the rest / harmonic lifts sit on it (chart u).
"""
import json, os, sys
import numpy as np
from pandakin import fk, kernel_vector, path_pose, pose_error, QMIN, QMAX, L, DATA, load_rest, constrained_ik

def wrap(q):
    return (q + np.pi) % (2 * np.pi) - np.pi

def project(q, Rt, pt):
    for _ in range(5):
        R, p, J = fk(q)
        e = pose_error(R, p, Rt, pt)
        if np.abs(e).max() < 1e-13:
            break
        q = q + J.T @ np.linalg.solve(J @ J.T, e)
    return q

if __name__ == "__main__":
    res = {}
    for k0 in (0, 250, 500, 750):
        s0, qrest, xr = [a[k0] for a in load_rest(1000)]
        Rt, pt = path_pose(s0)
        q = qrest.copy()
        x_prev = kernel_vector(fk(q)[2])
        h = 0.004
        tau = 0.0
        pts, taus, sig, inside = [q.copy()], [0.0], [], []
        closed = False
        for i in range(20000):
            def f(qq):
                x = kernel_vector(fk(qq)[2])
                return x if x @ x_prev >= 0 else -x   # keep a continuous orientation
            k1 = f(q); k2 = f(q + 0.5 * h * k1); k3 = f(q + 0.5 * h * k2); k4 = f(q + h * k3)
            q = project(q + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4), Rt, pt)
            x_prev = f(q)
            tau += h
            pts.append(q.copy()); taus.append(tau)
            if tau > 0.5 and np.linalg.norm(wrap(q - qrest)) < 0.6 * h:
                closed = True
                break
        P = np.array(pts)
        _, _, J = fk(P)
        sv = np.linalg.svd(J, compute_uv=False)[:, -1]
        Pw = wrap(P)
        ins = np.all((Pw >= QMIN) & (Pw <= QMAX), axis=1)
        # the admissible arc through the start: walk both ways from index 0 (closed curve: index wraps)
        n = len(P) - 1
        fwd = 0
        while fwd < n and ins[fwd + 1]:
            fwd += 1
        bwd = 0
        while bwd < n and ins[(-bwd - 1) % n]:
            bwd += 1
        viol = np.where(~ins)[0]
        first_bad = []
        if len(viol):
            jfw = viol[0]
            first_bad.append(dict(direction="+x", tau=float(taus[jfw]),
                                  joints=[int(j) + 1 for j in np.where((Pw[jfw] < QMIN) | (Pw[jfw] > QMAX))[0]]))
            jbw = viol[-1]
            first_bad.append(dict(direction="-x", tau=float(taus[-1] - taus[jbw]),
                                  joints=[int(j) + 1 for j in np.where((Pw[jbw] < QMIN) | (Pw[jbw] > QMAX))[0]]))
        out = dict(s0=float(s0), closed=closed, length=float(tau), sigma_min_along=float(sv.min()),
                   fraction_inside_limits=float(ins.mean()), admissible_arc_tau_forward=float(fwd * h),
                   admissible_arc_tau_backward=float(bwd * h), exits=first_bad)
        res["s%d" % k0] = out
        print(json.dumps(out), flush=True)
        np.savez(os.path.join(DATA, "fibre_k%d.npz" % k0), P=P, tau=np.array(taus), sig=sv, ins=ins)
    json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_fibre.json"), "w"), indent=1)
