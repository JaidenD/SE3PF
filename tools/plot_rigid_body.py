#!/usr/bin/env python3
"""Plot the rigid-body example (the paper's three panels) from the logs of apps/rigid_body_sim.

    python3 tools/plot_rigid_body.py [directory]        (default: out/ of the project)

Reads <directory>/rigid_body_pose.csv and <directory>/rigid_body_velocity.csv and writes <directory>/rigid_body.png
and <directory>/rigid_body.pdf:
  (a) transverse output: |xi| against the designed response |xi0| (1 + 4t) e^{-4t}, and the components xi_1..xi_5
      (pose-reference run; the transverse response does not depend on the tangential law);
  (b) pose reference: eta against eta_d, and the error e_eta;
  (c) velocity reference: eta_dot against nu_d, and the error e_v.
xi, eta and their rates are in the units of the metric W = I, in which radians and metres simply add.
Needs numpy and matplotlib.
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after choosing the backend)

# Colors: categorical slots in a fixed order (validated for colorblind separation), neutral ink for text, the
# measured outputs and the references.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK = "#0b0b0b"         # measured |xi|, eta, eta_dot
INK_2 = "#52514e"       # secondary text, axis labels
REFERENCE = "#8a8984"   # designed response and references (dashed)
GRID = "#e4e3df"

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

# Columns written by apps/rigid_body_sim.cpp
COLUMNS = ["t", "eta", "eta_dot", "ref", "e_tan", "xi1", "xi2", "xi3", "xi4", "xi5", "xi_norm", "xi_design",
           "u1", "u2", "u3", "u4", "u5", "u6"]


def load(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    missing = [c for c in COLUMNS if c not in data.dtype.names]
    if missing:
        sys.exit(f"{path}: missing columns {missing}; rerun rigid_body_sim")
    return data


def style(ax, ylabel, title):
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.set_ylabel(ylabel, color=INK_2)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_2, labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.axhline(0.0, color=INK_2, linewidth=0.6)


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "out")
    pose = load(os.path.join(directory, "rigid_body_pose.csv"))
    vel = load(os.path.join(directory, "rigid_body_velocity.csv"))

    fig, ax = plt.subplots(3, 1, figsize=(7.0, 8.4), constrained_layout=True)

    # (a) Transverse output. Shown over the first 5 s, where it decays; afterwards |xi| < 1e-3.
    a = ax[0]
    window = pose["t"] <= 5.0
    t = pose["t"][window]
    for i in range(5):
        a.plot(t, pose[f"xi{i + 1}"][window], color=SERIES[i], linewidth=1.2, label=rf"$\xi_{i + 1}$")
        # direct label at the start value, in text ink (the color is carried by the line)
        a.annotate(rf"$\xi_{i + 1}$", (0.0, pose[f"xi{i + 1}"][0]), xytext=(-14, 0), textcoords="offset points",
                   va="center", ha="right", fontsize=8, color=INK_2)
    a.plot(t, pose["xi_norm"][window], color=INK, linewidth=1.8, label=r"$|\xi|$")
    a.plot(t, pose["xi_design"][window], color=REFERENCE, linewidth=1.8, linestyle=(0, (4, 3)),
           label=r"design $|\xi_0|(1 + 4t)\,e^{-4t}$")
    dev_p = np.max(np.abs(pose["xi_norm"] - pose["xi_design"]))
    dev_v = np.max(np.abs(vel["xi_norm"] - vel["xi_design"]))
    a.text(0.99, 0.04, f"max $||\\xi| - $design$|$: {dev_p:.2e} (pose ref.), {dev_v:.2e} (velocity ref.)",
           transform=a.transAxes, ha="right", va="bottom", fontsize=8, color=INK_2)
    a.set_xlim(-0.45, 5.0)
    a.set_xlabel("t (s)", color=INK_2, fontsize=9)
    style(a, "transverse output (metric units)", r"(a) Transverse output $\xi$: critically damped, poles at $-4$")
    a.legend(ncol=4, fontsize=8, loc="upper right", frameon=False)

    # (b) Pose reference: eta follows eta_d(t) = L/4 + (L/5)(1 - cos(pi t/10)); eta stays in [0, L) here.
    b = ax[1]
    b.plot(pose["t"], pose["eta"], color=INK, linewidth=1.8, label=r"$\eta$")
    b.plot(pose["t"], pose["ref"], color=REFERENCE, linewidth=1.8, linestyle=(0, (4, 3)), label=r"$\eta_d$")
    b.plot(pose["t"], pose["e_tan"], color=SERIES[0], linewidth=1.5, label=r"$e_\eta = \eta - \eta_d$")
    b.set_xlim(0.0, 20.0)
    b.set_xlabel("t (s)", color=INK_2, fontsize=9)
    style(b, "arc length (metric units)", r"(b) Pose reference, $k_p = 9$, $k_d = 6$")
    b.legend(ncol=3, fontsize=8, loc="lower right", frameon=False)

    # (c) Velocity reference nu_d(t) = 0.5 + 0.2 sin(0.5 t)
    c = ax[2]
    c.plot(vel["t"], vel["eta_dot"], color=INK, linewidth=1.8, label=r"$\dot\eta$")
    c.plot(vel["t"], vel["ref"], color=REFERENCE, linewidth=1.8, linestyle=(0, (4, 3)), label=r"$\nu_d$")
    c.plot(vel["t"], vel["e_tan"], color=SERIES[0], linewidth=1.5, label=r"$e_v = \dot\eta - \nu_d$")
    c.set_xlim(0.0, 20.0)
    c.set_xlabel("t (s)", color=INK_2, fontsize=9)
    style(c, "speed (metric units / s)", r"(c) Velocity reference, $k_t = 5$")
    c.legend(ncol=3, fontsize=8, loc="lower right", frameon=False)

    for ext in ("png", "pdf"):
        out = os.path.join(directory, f"rigid_body.{ext}")
        fig.savefig(out, dpi=160, facecolor="white")
        print("wrote", out)


if __name__ == "__main__":
    main()
