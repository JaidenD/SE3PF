#!/usr/bin/env python3
"""Plot a self-motion run: the arm returning to the lift while the end effector stays on the path.

    python3 tools/plot_self_motion.py LOG.csv [--out FILE.png] [--design LAMBDA] [--sim]

LOG.csv is a log of robot_driver for orbit-selfmotion-hold or orbit-selfmotion-moving (or any run that starts off
the lift, beta0 != 0); its prediction LOG.csv.preflight.csv supplies the lift q*(eta): the predicted samples with
beta = 0. Panels:
  (a) beta against the prediction with the robot's model (the controller with its integral, no friction) and, with
      --design LAMBDA, the designed response beta0 (1 + LAMBDA t) e^{-LAMBDA t} of the PD (poles at -LAMBDA);
  (b) joints 1 and 3, which carry the self-motion, and the wrist joints 5 and 7, as q_i - q*_i(eta): the distance
      of each from the lift;
  (c) the end effector's position (mm) and attitude (deg) error against the path: what the self-motion must not
      disturb. With --sim, on a log scale with the elbow's distance from its lift position (mm) for scale: the error
      is at the level of numerical precision, and a linear axis would show an empty panel.
Writes FILE.png and FILE.pdf (default: LOG without .csv, plus _selfmotion). Needs numpy and matplotlib.
"""
import argparse
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from plot_robot_path import GRID, INK_2, SERIES, desired_pose, fk, joints, logSO3, path_length  # noqa: E402


def lift_of(pre, eta):
    """q*(eta) from the prediction's samples on the lift (|beta| < 1e-4), by interpolation in eta (wrapped)."""
    on = (pre["phase"] == 0) & (np.abs(pre["beta"]) < 1e-4)
    if not on.any():
        raise SystemExit("the prediction never reaches the lift (|beta| < 1e-4): no q* to compare with")
    pe, pq = pre["eta"][on], joints(pre, "q")[on]
    Lp = path_length(pre)
    e = eta if Lp is None else eta % Lp
    order = np.argsort(pe)
    if np.ptp(pe) < 1e-9:  # a hold: one point of the lift
        return np.tile(pq[-1], (len(eta), 1))
    return np.stack([np.interp(e, pe[order], pq[order][:, i]) for i in range(7)], axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("log")
    ap.add_argument("--out")
    ap.add_argument("--design", type=float, help="overlay the PD's designed response with poles at -DESIGN")
    ap.add_argument("--sim", action="store_true", help="the log is a simulation: no separate prediction curve")
    args = ap.parse_args()
    d = np.atleast_1d(np.genfromtxt(args.log, delimiter=",", names=True))
    pre_path = args.log + ".preflight.csv"
    if not os.path.exists(pre_path):
        raise SystemExit(f"{pre_path} not found: the lift comes from the prediction")
    pre = np.atleast_1d(np.genfromtxt(pre_path, delimiter=",", names=True))

    f = d["phase"] == 0
    t, eta, beta = d["t"][f], d["eta"][f], d["beta"][f]
    q = joints(d, "q")[f]
    qs = lift_of(pre, eta)
    R, p = fk(q)
    Rd, pd = desired_pose(d, f)
    posErr = 1e3 * np.linalg.norm(p - pd, axis=1)
    attErr = np.degrees(np.linalg.norm(logSO3(np.swapaxes(R, -1, -2) @ Rd), axis=1))

    beta0 = beta[0]
    print(f"{os.path.basename(args.log)}: {len(t)} following rows; beta {beta0:+.3f} -> {beta[-1]:+.4f} rad "
          f"(min {beta.min():+.4f}); end effector off the path up to {posErr.max():.1f} mm and "
          f"{attErr.max():.2f} deg (end: {posErr[-1]:.1f} mm, {attErr[-1]:.2f} deg)")
    off = q - qs
    print("  joints off the lift, start -> end (rad): " + ", ".join(
        f"{i + 1}: {off[0, i]:+.3f} -> {off[-1, i]:+.3f}" for i in range(7)))

    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})
    fig = plt.figure(figsize=(7, 8))

    # (a) beta
    axb = fig.add_subplot(3, 1, 1)
    axb.plot(t, beta, color=SERIES[0], label=r"$\beta$ (simulation)" if args.sim else r"$\beta$ (robot)")
    if not args.sim:
        axb.plot(pre["t"][pre["phase"] == 0], pre["beta"][pre["phase"] == 0], color=INK_2, linestyle=":",
                 label="prediction")
    if args.design:
        lam = args.design
        axb.plot(t, beta0 * (1 + lam * t) * np.exp(-lam * t), color="#0b0b0b", linestyle="--", linewidth=1.0,
                 label=r"design $\beta_0(1+%gt)e^{-%gt}$" % (lam, lam))
    axb.axhline(0.0, color=INK_2, linewidth=0.6)
    axb.set_ylabel("rad")
    axb.set_title(r"(a) Redundancy output $\beta$", loc="left")
    axb.grid(True, color=GRID); axb.legend(fontsize=8, frameon=False); axb.set_xlabel("$t$ (s)")

    # (b) joints off the lift
    axc = fig.add_subplot(3, 1, 2, sharex=axb)
    for i, col in zip((0, 2, 4, 6), (SERIES[0], SERIES[1], SERIES[2], SERIES[3])):
        axc.plot(t, np.degrees(off[:, i]), color=col, linewidth=1.2, label=f"joint {i + 1}")
    axc.axhline(0.0, color=INK_2, linewidth=0.6)
    axc.set_ylabel("deg")
    axc.set_title(r"(b) Joints off the lift, $q_i - q^\star_i(\eta)$", loc="left")
    axc.grid(True, color=GRID); axc.legend(fontsize=8, frameon=False, ncol=2); axc.set_xlabel("$t$ (s)")

    # (c) end-effector error
    axd = fig.add_subplot(3, 1, 3, sharex=axb)
    if args.sim:
        # Log scale: the error is numerically zero (about a micrometre). The elbow's distance from its lift
        # position, in the same unit (mm), shows how large the self-motion is against it.
        floor = 1e-7
        elbowOff = 1e3 * np.linalg.norm(fk(q, upTo=4)[1] - fk(qs, upTo=4)[1], axis=1)
        axd.semilogy(t, np.maximum(elbowOff, floor), color=INK_2, linestyle=":", label="elbow off the lift (mm)")
        axd.semilogy(t, np.maximum(posErr, floor), color=SERIES[0], label="end-effector position (mm)")
        axd.semilogy(t, np.maximum(attErr, floor), color=SERIES[2], label="end-effector attitude (deg)")
        axd.set_ylim(1e-5, 100.0)
        axd.text(0.98, 0.06, r"end effector: max %.1f $\mu$m, %.1f arcsec" % (1e3 * posErr.max(), 3600.0 * attErr.max()),
                 transform=axd.transAxes, ha="right", va="bottom", fontsize=8)
        axd.legend(fontsize=8, frameon=False, loc="upper right")
    else:
        axd.plot(t, posErr, color=SERIES[0], label="position (mm)")
        axd.plot(t, attErr, color=SERIES[2], label="attitude (deg)")
        axd.set_ylim(0.0, 1.1 * max(posErr.max(), attErr.max(), 1e-9))
        axd.legend(fontsize=8, frameon=False)
    axd.set_title("(c) End-effector error against the path", loc="left")
    axd.grid(True, color=GRID, which="major"); axd.set_xlabel("$t$ (s)")

    fig.tight_layout()
    base = args.out[:-4] if args.out and args.out.endswith(".png") else (
        args.out or os.path.splitext(args.log)[0] + "_selfmotion")
    fig.savefig(base + ".png", dpi=200)
    fig.savefig(base + ".pdf")
    print(f"wrote {base}.png and {base}.pdf")


if __name__ == "__main__":
    main()
