"""Task 1: discrete harmonic lifts for K = I and K = M, N = 250, 500, 1000, several starts.

Writes data/harm_<K>_<N>_<start>.npz and results_harmonic.json.
"""
import json, os, sys, time
import numpy as np
from hl import Problem, spectral_deriv, trig_interp
from pandakin import L, MjDyn, fk, kernel_vector, DATA

dyn = MjDyn()
starts = {"u0": 0.0, "up03": 0.3, "um03": -0.3, "up06": 0.6, "um06": -0.6}
res = {}
for met in ("I", "M"):
    for N in (250, 500, 1000):
        P = Problem(N, met, dyn=dyn if met == "M" else None)
        tol = 2e-9 * (N / 250) ** 2
        Erest = P.energy(P.qrest)
        rr, Fn_r, _ = P.el_residual(P.qrest)
        for name, c in starts.items():
            if N != 1000 and name in ("up06", "um06"):
                continue
            # the +-0.6 starts and random smooth starts only on the finest grid
            u0 = c + 0.0 * P.s
            t0 = time.time()
            try:
                u, q, E, hist = P.solve(u0, tol=tol, maxit=60)
            except RuntimeError as e:
                print(met, N, name, "FAILED:", e, flush=True)
                res["%s_%d_%s" % (met, N, name)] = dict(failed=str(e))
                continue
            dt = time.time() - t0
            g, _ = P.grad_u(u)
            H = P.hessian_u(u, q)
            t, _ = P.tangents(q)
            Kt = t if met == "I" else np.einsum("kij,kj->ki", dyn.Mbatch(q), t)
            D = P.h * np.einsum("ki,ki->k", t, Kt)          # |t_k|_K^2 h : the K-length weight of du_k
            Dm = 1 / np.sqrt(D)
            w, V = np.linalg.eigh(Dm[:, None] * H * Dm[None, :])
            wp = np.linalg.eigvalsh(H)
            # independent check of the smallest eigenvalue: energy second difference along its eigenvector
            phi = Dm * V[:, 0]
            phi /= np.sqrt(np.sum(D * phi ** 2))
            eps = 1e-3
            E0 = P.energy(q)
            Ep = P.energy(P.chart(u + eps * phi, q)); Em = P.energy(P.chart(u - eps * phi, q))
            ray = (Ep - 2 * E0 + Em) / eps ** 2
            r, Fn, dq = P.el_residual(q)
            r4 = P.el_residual_fd4(q)
            summ = P.summary(q)
            out = dict(metric=met, N=N, start=name, u0=c, time=dt, iterations=len(hist) - 1,
                       E=E, Erest=Erest, gradmax_over_h=float(np.abs(g).max() / P.h),
                       umin=float(u.min()), umax=float(u.max()), umean=float(u.mean()),
                       EL_spectral_max=float(np.abs(r).max()), EL_spectral_rel=float(np.abs(r).max() / Fn.max()),
                       EL_fd4_max=float(np.abs(r4).max()),
                       EL_rest_max=float(np.abs(rr).max()), EL_rest_rel=float(np.abs(rr).max() / Fn_r.max()),
                       lam_gen=[float(v) for v in w[:4]], lam_plain_min=float(wp[0]), lam_plain_max=float(wp[-1]),
                       lam_gen_check_energy2nddiff=float(ray),
                       chart_quality=float(P.chart_quality(q)), **summ)
            # periodicity of q*': one-sided 2nd-order differences at the seam (forward at s=0, backward at s=L)
            h = P.h
            fwd = (-3 * q[0] + 4 * q[1] - q[2]) / (2 * h)
            bwd = (3 * q[0] - 4 * q[-1] + q[-2]) / (2 * h)
            out["seam_derivative_jump"] = float(np.abs(fwd - bwd).max())
            res["%s_%d_%s" % (met, N, name)] = out
            np.savez(os.path.join(DATA, "harm_%s_%d_%s.npz" % (met, N, name)), s=P.s, u=u, q=q, qrest=P.qrest,
                     x=P.x, r=r, r4=r4, rr=rr, dq=dq, lam=w, eigvec0=phi, hist=np.array(hist))
            print("%s N=%4d %-5s it=%2d E=%.12f (rest %.9f) u in [%.4f, %.4f] |g|/h=%.1e EL=%.2e (fd4 %.2e) "
                  "lam_gen=%.5g (2nd diff %.5g) margin %.4f j%d sig %.4f max|q'| %.3f dist %.4f  %.1fs" % (
                      met, N, name, out["iterations"], E, Erest, u.min(), u.max(), out["gradmax_over_h"],
                      out["EL_spectral_max"], out["EL_fd4_max"], w[0], ray, summ["minMargin"], summ["marginJoint"],
                      summ["minSigma"], summ["maxDq"], summ["maxDistRest"], dt), flush=True)
        json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_harmonic.json"), "w"), indent=1)

# random smooth starts on N = 500 (uniqueness probe)
rng = np.random.default_rng(3)
for met in ("I", "M"):
    P = Problem(500, met, dyn=dyn if met == "M" else None)
    for trial in range(4):
        a = rng.uniform(-0.3, 0.3, 4)
        u0 = a[0] + a[1] * np.sin(2 * np.pi * P.s / L) + a[2] * np.cos(2 * np.pi * P.s / L) + a[3] * np.sin(
            4 * np.pi * P.s / L)
        try:
            u, q, E, hist = P.solve(u0, tol=8e-9, maxit=60)
        except RuntimeError as e:
            print(met, "random", trial, "FAILED", e)
            continue
        ref = np.load(os.path.join(DATA, "harm_%s_500_u0.npz" % met))
        d = float(np.abs(u - ref["u"]).max())
        res["%s_500_rand%d" % (met, trial)] = dict(E=E, coeffs=a.tolist(), max_du_vs_u0start=d,
                                                   iterations=len(hist) - 1)
        print("%s N=500 random start %d (%s): E=%.12f, max|u - u(start 0)| = %.2e" % (met, trial, np.round(a, 3), E,
                                                                                      d), flush=True)
json.dump(res, open(os.path.join(os.path.dirname(DATA), "results_harmonic.json"), "w"), indent=1)
