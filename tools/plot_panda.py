#!/usr/bin/env python3
"""Plot the Panda example (joint-acceleration level, nu_d = 0.08, 1 ms period) from the logs of apps/panda_sim.

    python3 tools/plot_panda.py [directory]        (default: out/ of the project)

Reads <directory>/panda_paper.csv, <directory>/panda_self_motion.csv and <directory>/panda_limits.csv (the joint
limits of PandaArm, written by panda_sim) and writes <directory>/panda.png and <directory>/panda.pdf, six panels:
  (a) transverse output of the paper run: |xi| against the designed response |xi0| (1 + 4t) e^{-4t}, and the
      components xi_1..xi_5;
  (b) the sampling effects of the paper run on a log scale: ||xi| - design| and |beta| (both O(dt));
  (c) the redundancy channel in the self-motion run: beta against beta0 (1 + 4t) e^{-4t}, while the tool stays on
      the path (max |xi| in the panel);
  (d) the tangential speed error |e_v| = |eta_dot - nu_d| on a log scale, both runs, with the designed decay
      |e_v(0)| e^{-k_t t} of the paper run (which starts at rest);
  (e) the distance of each joint to its nearer limit along the paper run, with the hardware margin of 0.10 rad;
  (f) |qd_i| / qd_max,i along the paper run, with the hardware limit of 50 %.
xi, eta and their rates are in the units of the metric W = diag(0.01 I, I): metres, with 1 rad counted as 0.1 m.
Needs numpy and matplotlib.
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after choosing the backend)

# Colors as in tools/plot_rigid_body.py: categorical slots in a fixed order (validated for colorblind separation;
# the light slots need visible labels, so the xi components are labelled directly), neutral ink for text, the
# measured outputs and the references. The joints are not seven hues: the joint that sets the extreme is
# highlighted in slot 1, the others are context in gray.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK = "#0b0b0b"         # measured |xi|, beta, eta_dot (paper run)
INK_2 = "#52514e"       # secondary text, axis labels
REFERENCE = "#8a8984"   # designed responses, references and limits (dashed)
CONTEXT = "#bdbcb6"     # the joints that are not highlighted
GRID = "#e4e3df"
DASH = (0, (4, 3))

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

# The joint limits are not copied here: panda_sim writes PandaArm's (src/panda.cpp) to panda_limits.csv. The two
# hardware margins below are HardwareBounds in include/se3pf/panda.hpp (keep them in sync); the joint-acceleration
# runs of panda_sim do not enforce them, so they are drawn for comparison.
MARGIN_HW = 0.10        # HardwareBounds::jointMargin, rad
VELOCITY_HW = 0.50      # HardwareBounds::dqFrac, fraction of qd_max
KT = 5.0                # tangential gain k_t: PandaExample::kt (include/se3pf/panda_example.hpp); the design line only

# Columns written by apps/panda_sim.cpp
COLUMNS = (["t", "eta", "eta_dot", "ref", "e_tan"] + [f"xi{i}" for i in range(1, 6)]
           + ["xi_norm", "xi_design", "beta", "beta_design"]
           + [f"{name}{i}" for name in ("q", "qd", "qdd") for i in range(1, 8)]
           + ["sigma_min", "joint_margin", "newton"])


def load(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    missing = [c for c in COLUMNS if c not in data.dtype.names]
    if missing:
        sys.exit(f"{path}: missing columns {missing}; rerun panda_sim")
    return data


def load_limits(path):
    """q_min, q_max and qd_max of joints 1..7, each of shape (7,)."""
    limits = np.genfromtxt(path, delimiter=",", names=True)
    if limits.shape != (7,) or not np.array_equal(limits["joint"], np.arange(1, 8)):
        sys.exit(f"{path}: expected joints 1..7; rerun panda_sim")
    return limits["q_min"], limits["q_max"], limits["qd_max"]


def joints(data, name):
    """The seven columns name1..name7 as an array of shape (rows, 7)."""
    return np.column_stack([data[f"{name}{i}"] for i in range(1, 8)])


def style(ax, ylabel, title, zero_line=True):
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.set_ylabel(ylabel, color=INK_2, fontsize=9)
    ax.set_xlabel("t (s)", color=INK_2, fontsize=9)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_2, labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    if zero_line:
        ax.axhline(0.0, color=INK_2, linewidth=0.6)


def note(ax, text, x=0.99, y=0.04, ha="right", va="bottom"):
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=8, color=INK_2)


def joint_margins(data, q_min, q_max, name):
    """Distance of each joint to its nearer limit, shape (rows, 7). Its row minimum must be the log's joint_margin
    column (PandaArm::jointMargin): if the limits file and the log came from different builds, stop."""
    q = joints(data, "q")
    margins = np.minimum(q - q_min, q_max - q)
    err = np.max(np.abs(np.min(margins, axis=1) - data["joint_margin"]))
    if err > 1e-8:  # the logs have 10 significant digits
        sys.exit(f"{name}: panda_limits.csv does not reproduce the joint_margin column (off by {err:.1e} rad); "
                 "rerun panda_sim")
    return margins


def highlight_joints(ax, t, values, joint, label):
    """All seven joints, the one that sets the extreme in slot 1 with a direct label, the others in gray."""
    for i in range(7):
        if i != joint:
            ax.plot(t, values[:, i], color=CONTEXT, linewidth=0.9)
    ax.plot(t, values[:, joint], color=SERIES[0], linewidth=1.8)
    ax.plot([], [], color=CONTEXT, linewidth=0.9, label="other joints")
    ax.plot([], [], color=SERIES[0], linewidth=1.8, label=label)


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "out")
    paper = load(os.path.join(directory, "panda_paper.csv"))
    selfm = load(os.path.join(directory, "panda_self_motion.csv"))
    q_min, q_max, qd_max = load_limits(os.path.join(directory, "panda_limits.csv"))
    margins = joint_margins(paper, q_min, q_max, "panda_paper.csv")
    joint_margins(selfm, q_min, q_max, "panda_self_motion.csv")  # checked only: the plot shows the paper run's
    t = paper["t"]

    fig, ax = plt.subplots(3, 2, figsize=(11.0, 9.6), constrained_layout=True)

    # (a) Transverse output of the paper run, over the first 3 s, where it decays
    a = ax[0, 0]
    window = t <= 3.0
    for i in range(5):
        a.plot(t[window], paper[f"xi{i + 1}"][window], color=SERIES[i], linewidth=1.2, label=rf"$\xi_{i + 1}$")
        # direct label at the start value, in text ink (the color is carried by the line)
        a.annotate(rf"$\xi_{i + 1}$", (0.0, paper[f"xi{i + 1}"][0]), xytext=(-12, 0), textcoords="offset points",
                   va="center", ha="right", fontsize=8, color=INK_2)
    a.plot(t[window], paper["xi_norm"][window], color=INK, linewidth=1.8, label=r"$|\xi|$")
    a.plot(t[window], paper["xi_design"][window], color=REFERENCE, linewidth=1.8, linestyle=DASH,
           label=r"design $|\xi_0|(1 + 4t)\,e^{-4t}$")
    a.set_xlim(-0.3, 3.0)
    style(a, "transverse output (metric units)", r"(a) Paper run: $\xi$ from $|\xi_0| = 0.0312$, poles at $-4$")
    a.legend(ncol=2, fontsize=8, loc="upper right", frameon=False)

    # (b) Deviations from the designed responses (log scale): both are sample-and-hold effects of order dt
    b = ax[0, 1]
    dev = np.abs(paper["xi_norm"] - paper["xi_design"])
    beta = np.abs(paper["beta"])
    k_dev, k_beta = np.argmax(dev), np.argmax(beta)
    # the two peaks are at the same instant, so their values go in the legend rather than beside the markers
    b.semilogy(t[dev > 0], dev[dev > 0], color=INK, linewidth=1.4,
               label=rf"$||\xi| - $design$|$, max {dev[k_dev]:.2e} at t = {t[k_dev]:.2f} s")
    b.semilogy(t[beta > 0], beta[beta > 0], color=SERIES[0], linewidth=1.4,
               label=rf"$|\beta|$ (design: 0), max {beta[k_beta]:.2e} at t = {t[k_beta]:.2f} s")
    b.plot([t[k_dev]], [dev[k_dev]], marker="o", markersize=4, color=INK)
    b.plot([t[k_beta]], [beta[k_beta]], marker="o", markersize=4, color=SERIES[0])
    b.set_xlim(0.0, 20.0)
    b.set_ylim(1e-9, 1e-3)
    style(b, "deviation (log scale)", "(b) Paper run: sampling effects, 1 ms zero-order hold", zero_line=False)
    b.legend(fontsize=8, loc="upper right", frameon=False)

    # (c) Redundancy channel: beta decays as designed while the tool stays on the path
    c = ax[1, 0]
    ts = selfm["t"]
    window = ts <= 3.0
    c.plot(ts[window], selfm["beta"][window], color=INK, linewidth=1.8, label=r"$\beta$")
    c.plot(ts[window], selfm["beta_design"][window], color=REFERENCE, linewidth=1.8, linestyle=DASH,
           label=r"design $\beta_0(1 + 4t)\,e^{-4t}$")
    dev_beta = np.max(np.abs(selfm["beta"] - selfm["beta_design"]))
    note(c, f"max $|\\beta - $design$|$ = {dev_beta:.2e}\nmax $|\\xi|$ = {np.max(selfm['xi_norm']):.1e}: "
            f"the tool stays on the path", y=0.3)
    c.set_xlim(0.0, 3.0)
    style(c, r"redundancy output $\beta$ (rad)", r"(c) Self-motion run: $\beta_0 = 0.2$ on the path, moving at $\nu_d$")
    c.legend(fontsize=8, loc="upper right", frameon=False)

    # (d) Tangential speed error e_v = eta_dot - nu_d, on a log scale so that both runs show: the paper run starts at
    # rest (e_v = -nu_d) and decays as designed, e_v(0) e^{-k_t t} (eq:tan-vel); the self-motion run starts at nu_d.
    # From about 2 s on both stay at a few 1e-6, a sampling effect of order dt like those in (b). As in (a) and (c),
    # the dashed design lies over the measured curve, and the thin self-motion line lies on top of both.
    d = ax[1, 1]
    ev_paper, ev_self = np.abs(paper["e_tan"]), np.abs(selfm["e_tan"])
    design = ev_paper[0] * np.exp(-KT * t)
    d.semilogy(t[ev_paper > 0], ev_paper[ev_paper > 0], color=INK, linewidth=2.0, zorder=1,
               label=r"$|e_v|$, paper run (from rest)")
    d.semilogy(t[design > 1e-12], design[design > 1e-12], color=REFERENCE, linewidth=1.8, linestyle=DASH, zorder=2,
               label=rf"design $|e_v(0)|\,e^{{-{KT:g}t}}$, $|e_v(0)| = \nu_d$")
    d.semilogy(ts[ev_self > 0], ev_self[ev_self > 0], color=SERIES[0], linewidth=1.1, zorder=3,
               label=r"$|e_v|$, self-motion run (starts at $\nu_d$)")
    k_sign = np.argmax(paper["e_tan"] > 0)  # first row with eta_dot above nu_d
    note(d, f"paper run: $e_v < 0$ until t = {t[k_sign]:.2f} s, then $> 0$\n"
            f"final $\\dot\\eta$ = {paper['eta_dot'][-1]:.7f}", y=0.04)
    d.set_xlim(0.0, 20.0)
    d.set_ylim(1e-9, 0.3)
    style(d, r"$|\dot\eta - \nu_d|$ (metric units / s, log scale)", rf"(d) Tangential speed error, $k_t = {KT:g}$",
          zero_line=False)
    d.legend(fontsize=8, loc="upper right", frameon=False)

    # (e) Distance of each joint to its nearer limit along the paper run
    e = ax[2, 0]
    k_min = np.unravel_index(np.argmin(margins), margins.shape)
    j = k_min[1]
    highlight_joints(e, t, margins, j,
                     f"joint {j + 1}: smallest margin {margins[k_min]:.4f} rad at t = {t[k_min[0]]:.2f} s")
    e.axhline(MARGIN_HW, color=REFERENCE, linewidth=1.5, linestyle=DASH, label="hardware margin 0.10 rad")
    e.plot([t[k_min[0]]], [margins[k_min]], marker="o", markersize=4, color=SERIES[0])
    e.set_xlim(0.0, 20.0)
    e.set_ylim(0.0, 4.0)  # room for the legend above the joints far from their limits
    style(e, "distance to the nearer limit (rad)", "(e) Paper run: joint margins")
    e.legend(fontsize=8, loc="upper right", frameon=False, ncol=1)

    # (f) Joint velocities as fractions of their limits
    f = ax[2, 1]
    qd_frac = 100.0 * np.abs(joints(paper, "qd")) / qd_max
    k_max = np.unravel_index(np.argmax(qd_frac), qd_frac.shape)
    j = k_max[1]
    highlight_joints(f, t, qd_frac, j, f"joint {j + 1}: largest {qd_frac[k_max]:.1f} % at t = {t[k_max[0]]:.2f} s")
    f.axhline(100.0 * VELOCITY_HW, color=REFERENCE, linewidth=1.5, linestyle=DASH, label="hardware limit 50 %")
    f.plot([t[k_max[0]]], [qd_frac[k_max]], marker="o", markersize=4, color=SERIES[0])
    f.set_xlim(0.0, 20.0)
    f.set_ylim(0.0, 60.0)
    style(f, r"$|\dot q_i|\,/\,\dot q_{\max,i}$ (%)", "(f) Paper run: joint velocities")
    f.legend(fontsize=8, loc="center right", frameon=False)

    for ext in ("png", "pdf"):
        out = os.path.join(directory, f"panda.{ext}")
        fig.savefig(out, dpi=160, facecolor="white")
        print("wrote", out)


if __name__ == "__main__":
    main()
