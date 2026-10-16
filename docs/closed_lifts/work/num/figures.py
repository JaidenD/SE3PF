"""Figures for the harmonic-lift numerics (PNG, in this directory)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pandakin import L, DATA, QMIN, QMAX, load_rest
from hl import trig_interp

HERE = os.path.dirname(os.path.abspath(__file__))
C = {"rest": "#2a78d6", "I": "#eb6834", "M": "#1baf7a"}   # validated categorical slots 1-3 (all-pairs)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.6,
                     "legend.frameon": False, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                     "savefig.facecolor": "#fcfcfb"})
lab = {"rest": "rest-posture lift", "I": "harmonic lift, K = I", "M": "harmonic lift, K = M"}
ls = {"rest": "-", "I": "--", "M": "-."}

H = {m: np.load(os.path.join(DATA, "harm_%s_1000_u0.npz" % m)) for m in ("I", "M")}
s = H["I"]["s"]
qrest = H["I"]["qrest"]

# ---------------------------------------------------------------- 1. joint trajectories
fig, axs = plt.subplots(4, 2, figsize=(7.2, 8.0), sharex=True)
for j in range(7):
    ax = axs.flat[j]
    ax.plot(s, qrest[:, j], color=C["rest"], ls=ls["rest"], label=lab["rest"])
    for m in ("I", "M"):
        ax.plot(s, H[m]["q"][:, j], color=C[m], ls=ls[m], label=lab[m])
    lo, hi = ax.get_ylim()
    for lim in (QMIN[j], QMAX[j]):
        if lo - 0.3 < lim < hi + 0.3:
            ax.axhline(lim, color=INK2, lw=0.8, ls=":")
            ax.text(s[-1], lim, " limit", va="center", ha="left", fontsize=7, color=INK2)
    ax.set_ylabel("q%d (rad)" % (j + 1))
axs.flat[6].set_xlabel("arc length s"); axs.flat[5].set_xlabel("arc length s")
axs.flat[7].axis("off")
h, l = axs.flat[0].get_legend_handles_labels()
axs.flat[7].legend(h, l, loc="center", fontsize=9)
fig.suptitle("Closed lifts of the Panda path (N = 1000 nodes)", color=INK, fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_joint_trajectories.png"), dpi=160)
plt.close(fig)

# ---------------------------------------------------------------- 2. beta offsets, vertical speed, E_k family
_, qr1000, x1000 = load_rest(1000)
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3))
ax = axs[0]
ax.axhline(0, color=C["rest"], ls=ls["rest"], label=lab["rest"])
for m in ("I", "M"):
    ax.plot(s, H[m]["u"], color=C[m], ls=ls[m], label=lab[m])
    Hh = np.load(os.path.join(DATA, "holonomy_%s.npz" % m))
    sg = Hh["sg"]
    xg = trig_interp(x1000, L, sg); qg = trig_interp(qr1000, L, sg)
    for key, sty, txt in (("traj_h", ":", "horizontal lift"), ("traj_ud", (0, (1, 1)), None)):
        uu = np.einsum("ki,ki->k", xg, Hh[key] - qg)
        if key == "traj_ud":
            ax.plot(sg, uu, color=C[m], ls="-", lw=0.7, alpha=0.8)
        if key == "traj_h":
            ax.plot(sg, uu, color=C[m], ls=":", lw=1.1)
            ax.annotate("horizontal lift (%s)" % ("K=I" if m == "I" else "K=M"), (sg[-1], uu[-1]), fontsize=7,
                        color=INK2, ha="right", va="top" if m == "I" else "bottom")
ax.set_xlabel("arc length s"); ax.set_ylabel("u = x^T (q - q*_rest)  (rad)")
ax.set_title("(a) offset from the rest-posture lift", fontsize=9, color=INK)
ax.text(0.30, 0.71, "thin: horizontal lift + uniform drift (closed)", fontsize=6.5, color=INK2)
ax.legend(fontsize=7, loc="center left")
ax = axs[1]
R = json.load(open(os.path.join(HERE, "results_holonomy.json")))
for m in ("I", "M"):
    Hh = np.load(os.path.join(DATA, "holonomy_%s.npz" % m))
    ax.plot(Hh["s"], Hh["v"], color=C[m], ls=ls[m], label="v(s), " + lab[m])
    ax.axhline(R[m]["uniform_drift_c"], color=C[m], lw=0.9, ls=":")
    ax.text(0.78, R[m]["uniform_drift_c"], "c (K=%s)" % m, fontsize=7, color=INK2, va="center", ha="left")
ax.set_xlabel("arc length s"); ax.set_ylabel("vertical speed K(n_K, q*')")
ax.set_title("(b) vertical velocity along the harmonic lifts", fontsize=9, color=INK)
ax.set_xlim(0, 0.86)
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1)
ax = axs[2]
F = json.load(open(os.path.join(HERE, "results_family.json")))
for m in ("I", "M"):
    ks = [1, 10, 100, 1e3, 1e4]
    dr = [F["%s_N500_k%g_up" % (m, k)]["max_dist_rest"] for k in ks]
    dh = [F["%s_N500_k%g_up" % (m, k)]["max_dist_harm"] for k in ks]
    ax.loglog(ks, dr, color=C[m], ls="-", marker="o", ms=4, label="to rest lift, K = %s" % m)
    ax.loglog(ks, dh, color=C[m], ls="--", marker="s", ms=4, mfc="none", label="to harmonic lift, K = %s" % m)
kk = np.array([100, 1e4])
ax.loglog(kk, 18 / kk, color=INK2, lw=0.8, ls=":")
ax.text(2.2e2, 18 / 1e3 * 0.35, "slope -1", fontsize=7, color=INK2)
ax.set_xlabel("penalty weight k"); ax.set_ylabel("max_s |q_k - q_ref|  (rad)")
ax.set_title("(c) minimizers of E + (k/2)|q - q_r|^2", fontsize=9, color=INK)
ax.legend(fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_offsets_drift_family.png"), dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- 3. convergence
RH = json.load(open(os.path.join(HERE, "results_harmonic.json")))
RD = json.load(open(os.path.join(HERE, "results_dae.json")))
Ns = np.array([250, 500, 1000])
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.2))
for m in ("I", "M"):
    E = np.array([RH["%s_%d_u0" % (m, N)]["E"] for N in Ns])
    Einf = E[2] + (E[2] - E[1]) / 3
    axs[0].loglog(Ns, abs(E - Einf), color=C[m], marker="o", ms=4, ls=ls[m], label=lab[m])
    el = [RH["%s_%d_u0" % (m, N)]["EL_spectral_max"] for N in Ns]
    axs[1].loglog(Ns, el, color=C[m], marker="o", ms=4, ls=ls[m], label=lab[m])
    dd = [RD[m]["periodic_orbit_vs_discrete_harmonic_N%d" % N] if N != 1000 else
          RD[m]["periodic_orbit_vs_discrete_harmonic_N1000_max_joint2norm"] for N in Ns]
    axs[2].loglog(Ns, dd, color=C[m], marker="o", ms=4, ls=ls[m], label=lab[m])
for ax, y0, t in ((axs[0], 3e-4, "(a) |E_N - E_Richardson|"), (axs[1], 4e-3, "(b) max |x^T (K q'' + Gamma(q',q'))|"),
                  (axs[2], 4e-5, "(c) max |q_N - q_periodic orbit|  (rad)")):
    ax.loglog(Ns, y0 * (Ns / 250.0) ** -2, color=INK2, lw=0.8, ls=":")
    ax.text(500, y0 * 4 ** -0.5 * 1.25, "h^2", fontsize=7, color=INK2)
    ax.set_xlabel("nodes N"); ax.set_title(t, fontsize=9, color=INK)
    ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_convergence.png"), dpi=160)
plt.close(fig)

# ---------------------------------------------------------------- 4. MuJoCo dragged arm
Mj = np.load(os.path.join(DATA, "mujoco_tc0.0005.npz"))
Mx = np.load(os.path.join(DATA, "mujoco_extra_dae_orbit_IC.npz"))
T1 = L / 0.08
lam = float(np.real(complex(RD["M"]["floquet_multipliers"][0])))
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3))
ax = axs[0]
k = "harmonic_M_lead2tc"
ax.semilogy(Mj[k + "_t"] / T1, Mj[k + "_dev"] + 1e-12, color=C["M"], label="start q*_harm(0), N = 1000 nodes")
ax.semilogy(Mx["t"] / T1, Mx["dev"] + 1e-12, color=C["rest"], ls="--", label="start on the exact orbit (DAE)\n(deviation from that orbit)")
k = "harmonic_M_perturbed_w+0.01_lead2tc"
ax.semilogy(Mj[k + "_t"] / T1, Mj[k + "_dev"] + 1e-12, color=C["I"], ls="-.", label="self-motion speed +0.01 at t = 0")
tt = np.linspace(1, 3, 10)
ax.semilogy(tt, 2e-6 * lam ** (tt - 1), color=INK2, lw=0.8, ls=":")
ax.text(1.9, 2e-6 * lam ** 0.9 / 12, "x %.2f per lap (DAE Floquet)" % lam, fontsize=7, color=INK2)
ax.set_xlabel("time (laps of 9.63 s)"); ax.set_ylabel("|q(t) - q*_harm(nu t)|  (rad)")
ax.set_title("(a) passive arm along the K = M harmonic lift", fontsize=9, color=INK)
ax.set_ylim(1e-9, 3)
ax.legend(fontsize=6.5, loc="lower right")
ax = axs[1]
k = "rest_lift_lead2tc"
ax.plot(Mj[k + "_t"], Mj[k + "_dev"], color=C["rest"], label="deviation from the rest-posture lift")
ax.set_xlabel("time (s)"); ax.set_ylabel("|q(t) - q*_rest(nu t)|  (rad)")
ax.set_title("(b) passive arm started on the rest-posture lift", fontsize=9, color=INK)
ax.legend(fontsize=7, loc="upper right")
ax = axs[2]
for k, col, st, txt in (("harmonic_M_lead0", C["I"], "-", "mocap at gamma(nu t)"),
                        ("harmonic_M_lead2tc", C["M"], "--", "mocap at gamma(nu (t + 2 tc))")):
    sel = Mj[k + "_t"] <= T1 + 1e-9
    ax.semilogy(Mj[k + "_t"][sel][1:], Mj[k + "_ep"][sel][1:], color=col, ls=st, label=txt)
ax.set_xlabel("time (s)"); ax.set_ylabel("|p_flange - p_gamma(nu t)|  (m)")
ax.set_title("(c) weld: flange error, tc = 0.5 ms", fontsize=9, color=INK)
ax.set_ylim(1e-9, 1e-3)
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_mujoco_drag.png"), dpi=160)
plt.close(fig)
print("figures written")
