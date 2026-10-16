#!/usr/bin/env python3
"""Plot the MuJoCo runs of the Panda example from the logs of apps/panda_mujoco (docs/MUJOCO.md).

    python3 tools/plot_mujoco.py [directory]        (default: out/ of the project)

Reads <directory>/mujoco_{A,B,C,C0,D}.csv and <directory>/panda_limits.csv (the torque limits, written by
panda_mujoco) and writes <directory>/mujoco.png and <directory>/mujoco.pdf:
  (a) |xi| against the designed response |xi0| (1 + 4t) e^{-4t}, runs A, B, C, D (first 5 s);
  (b) the deviation ||xi| - design| on a log scale, with the joint-acceleration-level baseline of run A;
  (c) |beta| on a log scale, with the 0.05 rad margin;
  (d) joint torques of run C as a fraction of their limits, with the 60 % margin;
  (e) torque rate of the commanded torque (without G), runs C and C0 (no ramp), with the 600 N m/s margin;
  (f) eta_dot against nu_d = 0.08.
xi and eta are in the units of the metric W = diag(0.01 I, I): metres, with 1 rad counted as 0.1 m.
Needs numpy and matplotlib.
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after choosing the backend)

# Colors as in tools/plot_rigid_body.py: categorical slots in a fixed order (validated for colorblind separation
# with the dataviz palette checker), neutral ink for text, gray for references. Each run keeps its color in every
# panel.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK = "#0b0b0b"
INK_2 = "#52514e"
REFERENCE = "#8a8984"   # designed responses, nu_d, the kinematic baseline (dashed)
LIMIT = "#52514e"       # safety margins (dotted)
GRID = "#e4e3df"

RUNS = {
    "A": ("A: exact model", SERIES[0]),
    "B": ("B: zero-order hold", SERIES[1]),
    "C": ("C: hardware-like (hold, ramp, damping)", SERIES[2]),
    "D": ("D: hand not modeled", SERIES[3]),
    "C0": ("C0: as C, no ramp", SERIES[4]),
}

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

# The hardware bounds drawn as limit lines: HardwareBounds in include/se3pf/panda.hpp (keep them in sync). The joint
# limits themselves are not copied here: panda_mujoco writes PandaArm's to panda_limits.csv.
XI_MAX = 0.04            # HardwareBounds::xi
BETA_MAX = 0.05          # HardwareBounds::beta, rad
TAU_FRACTION = 0.6       # HardwareBounds::tauFrac
TAU_RATE_MAX = 600.0     # HardwareBounds::tauRateFrac * PandaArm::tauRateMax, N m/s


def load_tau_max(directory):
    """tau_max of joints 1..7 (N m), from panda_limits.csv."""
    path = os.path.join(directory, "panda_limits.csv")
    if not os.path.exists(path):
        sys.exit(f"{path} not found; run  ./build-rel/panda_mujoco {directory}  first")
    limits = np.genfromtxt(path, delimiter=",", names=True)
    if limits.shape != (7,) or "tau_max" not in limits.dtype.names \
            or not np.array_equal(limits["joint"], np.arange(1, 8)):
        sys.exit(f"{path}: expected joints 1..7 with a tau_max column; rerun panda_mujoco")
    return limits["tau_max"]


def load(directory, run):
    path = os.path.join(directory, f"mujoco_{run}.csv")
    if not os.path.exists(path):
        sys.exit(f"{path} not found; run  ./build-rel/panda_mujoco {directory}  first")
    with open(path) as f:
        names = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return {name: data[:, i] for i, name in enumerate(names)}


def style(ax, ylabel, title, xlabel="t (s)"):
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.set_ylabel(ylabel, color=INK_2, fontsize=9)
    ax.set_xlabel(xlabel, color=INK_2, fontsize=9)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_2, labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)


def limit_line(ax, y, text, x=0.99):
    ax.axhline(y, color=LIMIT, linewidth=1.0, linestyle=(0, (1, 2)))
    ax.annotate(text, (x, y), xycoords=("axes fraction", "data"), xytext=(0, 3), textcoords="offset points",
                ha="right", va="bottom", fontsize=7.5, color=INK_2)


def positive(x, floor):
    """For a log scale: values below the floor (exact zeros at t = 0) are not drawn."""
    y = np.abs(x).astype(float)
    y[y < floor] = np.nan
    return y


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "out")
    runs = {name: load(directory, name) for name in RUNS}
    tau_max = load_tau_max(directory)
    four = ["A", "B", "C", "D"]

    fig, ax = plt.subplots(3, 2, figsize=(11.0, 10.5), constrained_layout=True)

    # (a) |xi| against its designed response, first 5 s
    a = ax[0, 0]
    for name in four:
        r = runs[name]
        w = r["t"] <= 5.0
        a.plot(r["t"][w], r["xi_norm"][w], color=RUNS[name][1], linewidth=1.6, label=RUNS[name][0])
    r = runs["A"]
    w = r["t"] <= 5.0
    a.plot(r["t"][w], r["xi_design"][w], color=REFERENCE, linewidth=1.8, linestyle=(0, (4, 3)),
           label=r"design $|\xi_0|(1 + 4t)\,e^{-4t}$")
    limit_line(a, XI_MAX, rf"margin $|\xi| \leq {XI_MAX:g}$")
    a.set_xlim(0.0, 5.0)
    a.set_ylim(0.0, 0.25)
    style(a, r"$|\xi|$ (metric units)", r"(a) Transverse output $|\xi|$ and its designed response")
    a.text(0.99, 0.55, "A and B lie on the design;\n"
           "C: ramp over 1 s, then a steady offset\nfrom the uncompensated damping;\n"
           "D sags under the hand's weight",
           transform=a.transAxes, ha="right", va="center", fontsize=8, color=INK_2)

    # (b) deviation from the design, log scale
    b = ax[0, 1]
    for name in four:
        r = runs[name]
        b.semilogy(r["t"], positive(r["xi_norm"] - r["xi_design"], 1e-10), color=RUNS[name][1], linewidth=1.0)
    r = runs["A"]
    b.semilogy(r["t"], positive(r["xi_norm_kin"] - r["xi_design"], 1e-10), color=REFERENCE, linewidth=1.0,
               linestyle=(0, (4, 3)), label="joint-acceleration level (design doc)")
    b.set_xlim(0.0, 20.0)
    b.set_ylim(1e-9, 1.0)
    style(b, r"$||\xi| - $design$|$", r"(b) Deviation from the designed response")
    b.text(0.99, 0.03, "A: %.2e   B: %.2e   baseline: %.2e (max)" % (
        np.max(np.abs(runs["A"]["xi_norm"] - runs["A"]["xi_design"])),
        np.max(np.abs(runs["B"]["xi_norm"] - runs["B"]["xi_design"])),
        np.max(np.abs(r["xi_norm_kin"] - r["xi_design"]))),
        transform=b.transAxes, ha="right", va="bottom", fontsize=7.5, color=INK_2)
    b.legend(fontsize=7.5, loc="upper right", frameon=False)

    # (c) |beta|, log scale
    c = ax[1, 0]
    for name in four:
        r = runs[name]
        c.semilogy(r["t"], positive(r["beta"], 1e-10), color=RUNS[name][1], linewidth=1.0)
    limit_line(c, BETA_MAX, rf"margin $|\beta| \leq {BETA_MAX:g}$ rad", x=0.6)
    c.set_xlim(0.0, 20.0)
    c.set_ylim(1e-8, 0.3)
    style(c, r"$|\beta|$ (rad)", r"(c) Redundancy output $|\beta|$ (zero at the start)")

    # (d) joint torques of run C as a fraction of their limits. Joints are ordered, so an ordinal ramp of one hue
    # (light to dark) with direct labels; the run color is not used here.
    d = ax[1, 1]
    r = runs["C"]
    ramp = ["#a9c9ef", "#86b2e8", "#5f97dd", "#3f80d4", "#2a6bbd", "#1f5598", "#163f73"]
    ends = []
    for i in range(7):
        frac = r[f"tau{i + 1}"] / tau_max[i]
        d.plot(r["t"], frac, color=ramp[i], linewidth=1.1)
        ends.append((frac[-1], i))
    # Direct labels at the line ends, spread apart where several lines end close together: push each label up to
    # keep a minimum gap, then shift the stack down by half its overflow, so the labels stay centred on the lines.
    gap = 0.055
    ends.sort()
    placed = []
    for y, _ in ends:
        placed.append(max(y, placed[-1] + gap) if placed else y)
    shift = 0.5 * (placed[-1] - ends[-1][0])
    for (y, i), y_label in zip(ends, placed):
        d.annotate(f"joint {i + 1}", (r["t"][-1], y), xytext=(r["t"][-1] + 0.5, y_label - shift), textcoords="data",
                   va="center", fontsize=7, color=INK_2,
                   arrowprops=dict(arrowstyle="-", color=GRID, linewidth=0.6, shrinkA=0, shrinkB=0))
    for y in (TAU_FRACTION, -TAU_FRACTION):
        limit_line(d, y, rf"margin $|\tau_i| \leq {TAU_FRACTION:g}\,\tau_{{\max,i}}$" if y > 0 else "")
    d.set_xlim(0.0, 23.0)
    d.set_ylim(-0.7, 0.7)
    style(d, r"$\tau_i / \tau_{\max,i}$", "(d) Run C: joint torques (gravity included) against their limits")

    # (e) torque rate of the commanded torque (tau - G, what libfranka is sent), first 3 s
    e = ax[2, 0]
    for name in ("C0", "C"):
        r = runs[name]
        w = r["t"] <= 3.0
        e.semilogy(r["t"][w], positive(r["tau_rate_cmd"][w], 1e-6), color=RUNS[name][1], linewidth=1.2)
    limit_line(e, TAU_RATE_MAX, "margin 60 % of 1000 N m/s")
    r0 = runs["C0"]
    e.annotate("C0: %.0f N m/s at t = 0\n(the first command)" % r0["tau_rate_cmd"][0], (0.0, r0["tau_rate_cmd"][0]),
               xytext=(14, -22), textcoords="offset points", va="top", fontsize=7.5, color=INK_2)
    e.set_xlim(-0.05, 3.0)
    e.set_ylim(1e-3, 1e4)
    style(e, r"max$_i\,|d(\tau_i - G_i)/dt|$ (N m/s)", "(e) Torque rate: the ramp removes the jump at t = 0")

    # (f) eta_dot against nu_d
    f = ax[2, 1]
    for name in four:
        r = runs[name]
        f.plot(r["t"], r["eta_dot"], color=RUNS[name][1], linewidth=1.2)
    f.axhline(0.08, color=REFERENCE, linewidth=1.8, linestyle=(0, (4, 3)))
    f.annotate(r"$\nu_d = 0.08$", (0.99, 0.08), xycoords=("axes fraction", "data"), xytext=(0, 3),
               textcoords="offset points", ha="right", va="bottom", fontsize=7.5, color=INK_2)
    f.set_xlim(0.0, 20.0)
    f.set_ylim(-0.02, 0.45)
    style(f, r"$\dot\eta$ (metric units / s)", r"(f) Tangential velocity $\dot\eta$ (damping and payload slow it)")

    # One legend for the runs, shared by every panel: each run keeps its color throughout.
    handles = [plt.Line2D([], [], color=RUNS[n][1], linewidth=2.0, label=RUNS[n][0]) for n in ("A", "B", "C", "D",
                                                                                                 "C0")]
    handles.append(plt.Line2D([], [], color=REFERENCE, linewidth=1.8, linestyle=(0, (4, 3)),
                              label="design / reference / kinematic baseline"))
    fig.legend(handles=handles, loc="outside upper center", ncol=3, fontsize=8.5, frameon=False)

    for ext in ("png", "pdf"):
        out = os.path.join(directory, f"mujoco.{ext}")
        fig.savefig(out, dpi=150, facecolor="white")
        print("wrote", out)


if __name__ == "__main__":
    main()
