"""Figure for the comparison of P (SE(3) tube + beta, the paper), J (Euclidean joint-space tube around q*) and
L (joint-space coordinates around q* linearized from the paper's outputs at the lift; the reviewer's strongest J).
Rows S1..S4 (time responses), then the S2 offset sweep and the S3r generic pushes."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = sys.argv[1] if len(sys.argv) > 1 else "out_rel"
out = sys.argv[2] if len(sys.argv) > 2 else "compare.png"

# categorical slots 1-3 of the reference palette (validated: CVD and normal-vision separation pass on #fcfcfb);
# L also gets a dashed line, so identity is never color alone
C = {"P": "#2a78d6", "J": "#eb6834", "L": "#1baf7a"}
LS = {"P": "-", "J": "-", "L": (0, (5, 2))}
NAME = {"P": "P: SE(3) tube (eta, xi) + beta (paper)",
        "J": "J: Euclidean joint-space tube around q*(s)",
        "L": "L: joint coordinates around q*(s), linearized from (eta, xi, beta)"}
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.titlecolor": INK, "text.color": INK, "axes.facecolor": SURF, "figure.facecolor": SURF,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
})


def load(name):
    return np.genfromtxt(f"{d}/{name}.csv", delimiter=",", names=True)


rows = [
    ("S1", "S1  paper start: tool off the path (|xi0| = 0.031), beta = 0, at rest", None),
    ("S2", "S2  tool on the path, elbow off the lift (beta0 = 0.2 rad), moving at nu_d", None),
    ("S3", "S3  on the lift at nu_d; elbow push 0.3 rad/s along ker J_b at t = 5 s", 5.0),
    ("S4", "S4  all three given a lift that is wrong by 0.02 rad (joints 2, 4)", None),
]
FLOOR = 1e-9

fig, axes = plt.subplots(5, 2, figsize=(11.5, 15.5), constrained_layout=True)
for i, (sc, title, kick) in enumerate(rows):
    R = {c: load(f"{sc}_{c}") for c in "PJL"}
    axL, axR = axes[i]
    for c in "PJL":
        axL.semilogy(R[c]["t"], np.maximum(R[c]["xiSE3"], FLOOR), color=C[c], ls=LS[c], lw=1.6, label=NAME[c])
        axR.semilogy(R[c]["t"], np.maximum(R[c]["jdist"], FLOOR), color=C[c], ls=LS[c], lw=1.6)
    P = R["P"]
    if sc == "S1":
        t = P["t"]
        x0 = P["xiSE3"][0]
        axL.semilogy(t, np.maximum(x0 * (1 + 4 * t) * np.exp(-4 * t), FLOOR), color=INK2, lw=1.0, ls=(0, (1, 2)),
                     label="designed response |xi0|(1+4t)e^{-4t}")
    axL.set_ylim(1e-8, 0.1)
    axL.set_ylabel("tool error |xi_SE3| (W-distance to path)")
    axL.set_title(title, loc="left", fontsize=9.5)
    axR.set_ylim(1e-8, 1.0)
    axR.set_ylabel("joint distance to true lift, rad")
    axR.set_title(f"{sc}  joint-space distance to the true lift q*", loc="left", fontsize=9.5)
    for ax in (axL, axR):
        ax.set_xlim(0, P["t"][-1])
        ax.set_xlabel("t, s")
        if kick is not None:
            ax.axvline(kick, color=MUTED, lw=0.9, ls=":")
            ax.text(kick + 0.1, ax.get_ylim()[1] * 0.4, "push", color=INK2, fontsize=8, va="top")

    if sc in ("S2", "S3"):
        for c, ytxt, dx in (("J", 3e-3, 1.5), ("L", 3e-6 if sc == "S3" else 2e-5, 2.5 if sc == "S3" else 3.0)):
            k = np.argmax(R[c]["xiSE3"])
            axL.annotate(f"{c} peak {R[c]['xiSE3'][k]:.1e} ({R[c]['posErr'][k]*1e3:.2f} mm, {R[c]['attErrDeg'][k]:.2g} deg)",
                         (R[c]["t"][k], R[c]["xiSE3"][k]), xytext=(R[c]["t"][k] + dx, ytxt), color=INK2, fontsize=8,
                         arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
        axL.text(P["t"][-1] * 0.98, 3e-8, "P stays at the O(dt) sampling floor (~1e-6)", color=INK2, fontsize=8, ha="right")
    if sc == "S4":
        m = R["J"]["t"] >= 15
        axL.text(P["t"][-1] * 0.98, 4e-2, f"J and L steady: up to {R['J']['xiSE3'][m].max():.3f} "
                 f"({R['J']['posErr'][m].max()*1e3:.1f} mm, {R['J']['attErrDeg'][m].max():.2f} deg)",
                 color=INK2, fontsize=8, ha="right")
        axL.text(P["t"][-1] * 0.98, 3e-8, "P: tool error unaffected by the lift error", color=INK2, fontsize=8, ha="right")
        axR.text(P["t"][-1] * 0.98, 0.25, "J and L converge to the wrong lift (0.02 rad)", color=INK2, fontsize=8, ha="right")
    if sc == "S1":
        axR.text(P["t"][-1] * 0.98, 0.25, "all three reach the lift; curves overlap", color=INK2, fontsize=8, ha="right")

# Row 5 left: S2 sweep, peak tool error against the elbow offset
sw = np.genfromtxt(f"{d}/sweep.csv", delimiter=",", names=True)
ax = axes[4, 0]
for c in "PJL":
    ax.loglog(sw["beta0"], sw[c], color=C[c], ls=LS[c], lw=1.6, marker="o", ms=4.5,
              markeredgecolor=SURF, markeredgewidth=1.0)
b = np.array([0.01, 0.8])
ax.loglog(b, 3.0e-3 * b, color=MUTED, lw=0.8, ls=(0, (1, 2)))
ax.loglog(b, 7.0e-3 * b ** 2, color=MUTED, lw=0.8, ls=(0, (1, 2)))
ax.text(0.012, 3.0e-3 * 0.012 * 1.6, "slope 1", color=INK2, fontsize=8)
ax.text(0.035, 7.0e-3 * 0.035 ** 2 * 0.4, "slope 2", color=INK2, fontsize=8)
ax.text(0.79, 1.6e-6, "P: floor (~1.2e-6)", color=INK2, fontsize=8, ha="right", va="bottom")
ax.set_ylim(5e-7, 2e-2)
ax.set_xlabel("elbow offset beta0 (tool on the path), rad")
ax.set_ylabel("peak tool error |xi_SE3|")
ax.set_title("S2 sweep  peak tool error against the elbow offset", loc="left", fontsize=9.5)

# Row 5 right: S3r generic pushes, peak tool error relative to P's
r3 = np.genfromtxt(f"{d}/S3r.csv", delimiter=",", names=True)
ax = axes[4, 1]
idx = np.arange(1, len(r3) + 1)
ax.axhline(1.0, color=C["P"], lw=1.6)
ax.plot(idx, r3["J"] / r3["P"], color=C["J"], ls="none", marker="o", ms=6, markeredgecolor=SURF, markeredgewidth=1.0,
        label="J / P")
ax.plot(idx, r3["L"] / r3["P"], color=C["L"], ls="none", marker="D", ms=5.5, markeredgecolor=SURF, markeredgewidth=1.0,
        label="L / P")
ax.legend(loc="lower left", frameon=False, fontsize=8, ncol=2)
ax.set_ylim(0.9, 1.1)
ax.set_xlim(0.3, len(r3) + 0.7)
ax.set_xticks(idx[::2])
ax.set_xlabel("push number (random joint-velocity direction, 0.3 rad/s)")
ax.set_ylabel("peak tool error / P's peak")
ax.set_title("S3r  generic pushes: peak tool error relative to P", loc="left", fontsize=9.5)
ax.text(len(r3) + 0.5, 1.085, f"J/P {np.min(r3['J']/r3['P']):.2f} to {np.max(r3['J']/r3['P']):.2f};  "
        f"L/P {np.min(r3['L']/r3['P']):.3f} to {np.max(r3['L']/r3['P']):.3f}", color=INK2, fontsize=8, ha="right")

h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="outside upper center", ncol=2, frameon=False, fontsize=9)
fig.savefig(out, dpi=130)

# extra numbers for the report
for sc in ("S1", "S2", "S2small", "S3", "S4"):
    for c in "PJL":
        r = load(f"{sc}_{c}")
        t, x = r["t"], r["xiSE3"]

        def first_below(th):
            above = np.nonzero(x > th)[0]
            return t[above[-1] + 1] if len(above) and above[-1] + 1 < len(t) else (0.0 if not len(above) else float("nan"))
        print(f"{sc} {c}: time to |xi_SE3| < 1e-3: {first_below(1e-3):.3f} s, < 1e-4: {first_below(1e-4):.3f} s; "
              f"jdist final {r['jdist'][-1]:.2e}")
