"""Reduced picture on the self-motion bundle: in the chart (s, u), a closed lift has energy density
    1/2 |q'|_K^2 = 1/2 a (u' - A)^2 + V,   V(s, u) = 1/2 |J^# zeta_gamma|_K^2 = 1/2 zeta^T (J K^-1 J^T)^-1 zeta,
a = |dq/du|_K^2, A = the horizontal velocity in the chart (u' of the horizontal lift). V is the horizontal kinetic
energy density (operational-space inertia along zeta). This script tabulates V, a, the joint margin and the chart
quality on a grid of (s, u) and plots V with the critical points found (near and far harmonic lifts)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pandakin import fk, kernel_vector, path_pose, path_vel, constrained_ik, load_rest, joint_margin, MjDyn, L, DATA

dyn = MjDyn()
s_all, q_all, x_all = load_rest(1000)
idx = np.arange(0, 1000, 10)
s, qr, xr = s_all[idx], q_all[idx], x_all[idx]
Rt, pt = path_pose(s)
z = path_vel(s)
us = np.linspace(-3.0, 3.0, 121)
i0 = int(np.argmin(abs(us)))
nS, nU = len(s), len(us)
V = {m: np.full((nU, nS), np.nan) for m in ("I", "M")}
a = {m: np.full((nU, nS), np.nan) for m in ("I", "M")}
margin = np.full((nU, nS), np.nan)
cq = np.full((nU, nS), np.nan)
b0 = np.einsum("ki,ki->k", xr, qr)
for direction in (1, -1):
    q = qr.copy()
    rng = range(i0, nU) if direction == 1 else range(i0, -1, -1)
    for i in rng:
        qn, ok = constrained_ik(Rt, pt, xr, b0 + us[i], q, tol=1e-12)
        from pandakin import pose_error
        Rq, pq, _ = fk(qn)
        res_n = np.maximum(np.abs(pose_error(Rq, pq, Rt, pt)).max(1), np.abs(np.einsum("ki,ki->k", xr, qn) - b0 - us[i]))
        good = res_n < 1e-9
        if not good.any():
            print("chart solve failed everywhere at u = %.2f" % us[i]); break
        q = np.where(good[:, None], qn, q)
        bad = ~good
        _, _, J = fk(q)
        x = kernel_vector(J)
        cq[i] = np.abs(np.einsum("ki,ki->k", x, xr))
        margin[i] = joint_margin(q)[0]
        t = np.linalg.solve(np.concatenate([J, xr[:, None, :]], 1), np.r_[np.zeros(6), 1.0][None, :, None].repeat(nS, 0))[..., 0]
        JJ = J @ np.swapaxes(J, 1, 2)
        V["I"][i] = 0.5 * np.einsum("ki,ki->k", z, np.linalg.solve(JJ, z[..., None])[..., 0])
        a["I"][i] = np.einsum("ki,ki->k", t, t)
        M = dyn.Mbatch(q)
        Minv = np.linalg.inv(M)
        Lam = J @ Minv @ np.swapaxes(J, 1, 2)
        V["M"][i] = 0.5 * np.einsum("ki,ki->k", z, np.linalg.solve(Lam, z[..., None])[..., 0])
        a["M"][i] = np.einsum("ki,kij,kj->k", t, M, t)
        for arr in (V["I"][i], V["M"][i], a["I"][i], a["M"][i], margin[i], cq[i]):
            arr[bad] = np.nan
np.savez(os.path.join(DATA, "landscape.npz"), s=s, us=us, VI=V["I"], VM=V["M"], aI=a["I"], aM=a["M"], margin=margin, cq=cq)

INK2 = "#52514e"
fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.8), sharey=True)
for ax, m, title in ((axs[0], "I", "K = I:  V = 1/2 |J^+ zeta|^2"), (axs[1], "M", "K = M:  V = 1/2 zeta^T Lambda zeta")):
    Vm = np.ma.masked_invalid(np.where(margin >= 0, V[m], np.nan))
    pc = ax.pcolormesh(s, us, Vm, shading="auto", cmap="Blues")
    fig.colorbar(pc, ax=ax, label="V(s, u)")
    ax.contour(s, us, margin, levels=[0.0, 0.1], colors=[INK2, INK2], linestyles=["-", ":"], linewidths=0.9)
    ax.plot(s_all, 0 * s_all, color="#2a78d6", lw=1.4, label="rest-posture lift")
    H = np.load(os.path.join(DATA, "harm_%s_1000_u0.npz" % m))
    ax.plot(H["s"], H["u"], color="#eb6834", lw=1.6, ls="--", label="harmonic lift near rest (local min)")
    for c, sty in ((-2.0, "-."), (2.0, ":")):
        f = os.path.join(DATA, "far_%s_start%+g.npz" % (m, c))
        if os.path.exists(f):
            F = np.load(f)
            ax.plot(F["s"], F["u"], color="#1baf7a", lw=1.6, ls=sty, label="critical point from u0 = %+g" % c)
    f = os.path.join(DATA, "limits_%s_far-2.npz" % m)
    if os.path.exists(f):
        F = np.load(f)
        ax.plot(F["s"], F["u"], color=INK2, lw=1.0, ls="-", label="min with 0.10 rad margin enforced (from u0 = -2)")
    ax.set_xlabel("arc length s"); ax.set_title(title, fontsize=9)
    ax.set_xlim(0, L); ax.set_ylim(-3.0, 3.0)
axs[0].set_ylabel("u = x^T (q - q*_rest)  (rad)")
axs[1].legend(fontsize=7, loc="upper center", bbox_to_anchor=(-0.1, -0.16), ncol=3)
fig.text(0.5, -0.2, "white: outside the joint limits (or the linear chart folds); solid / dotted contour: joint margin 0 / 0.10 rad",
         ha="center", fontsize=7, color=INK2)
fig.savefig(os.path.join(os.path.dirname(DATA), "fig_landscape.png"), dpi=160, bbox_inches="tight")
print("chart quality min over admissible grid: %.3f" % np.nanmin(np.where(margin >= 0, cq, np.nan)))
for m in ("I", "M"):
    for lab, uu in (("u=0", 0.0), ("u=-2.3", -2.3)):
        i = int(np.argmin(abs(us - uu)))
        print(m, lab, "mean_s V = %.4f" % np.nanmean(V[m][i]))
