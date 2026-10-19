#!/usr/bin/env python3
"""Plot a run of the robot driver (docs/HARDWARE.md), with its predictions.

    python3 tools/plot_robot_log.py LOG.csv [--overlay PRED.csv[:label]]... [--no-preflight] [--out FILE.png]

LOG.csv is a log of robot_driver (or robot_driver_mock); its columns are in apps/robot/robot_log.hpp (a log without
the tauF columns, from before CTRL-1, still plots). The preflight's
prediction LOG.csv.preflight.csv, which robot_driver writes next to every log, is overlaid automatically (dashed;
--no-preflight leaves it out); --overlay adds more, e.g. the MuJoCo closed loop of the same experiment
(out/mujoco_paper_hand.csv from experiments_mujoco). Writes FILE.png and FILE.pdf (default: next to LOG.csv) with
  (a) |xi| and its bound 0.04;            (b) beta and its bound 0.05 rad;
  (c) eta_dot against nu_d (velocity reference) or the tangential error eta - eta_d (position reference);
  (d) the delivered torque |tau_i + G_i| / tau_max,i of every joint, with the 60 % bound;
  (e) the distance of every joint to its nearer limit, with the 0.10 rad bound;
  (f) the largest commanded |qdd*_i| / qdd_max,i and |qd_i| / qd_max,i, with the 30 % and 50 % bounds;
and prints the peaks, the safety stop (if any), the steady errors and the timing. A safety stop (phase 1: the damped
hold at the configuration of the violation) is shaded in every panel. Needs numpy and matplotlib.
"""
import argparse
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after choosing the backend)

# Colors as in tools/plot_mujoco.py
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#7a5cc9", "#4c4b47"]
INK = "#0b0b0b"
INK_2 = "#52514e"
REFERENCE = "#8a8984"
LIMIT = "#52514e"
GRID = "#e4e3df"
STOP = "#f3c9b5"

# PandaArm's limits (src/panda.cpp) and HardwareBounds (include/se3pf/panda.hpp): keep them in sync
Q_MIN = np.array([-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973])
Q_MAX = np.array([2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973])
DQ_MAX = np.array([2.1750, 2.1750, 2.1750, 2.1750, 2.6100, 2.6100, 2.6100])
DDQ_MAX = np.array([15.0, 7.5, 10.0, 12.5, 15.0, 20.0, 20.0])
TAU_MAX = np.array([87.0, 87.0, 87.0, 87.0, 12.0, 12.0, 12.0])
JOINT_MARGIN, DQ_FRAC, DDQ_FRAC, TAU_FRAC = 0.10, 0.5, 0.3, 0.6
TAU_RATE_MAX = 600.0
XI_MAX, BETA_MAX = 0.04, 0.05


def load(path):
    if not os.path.exists(path):
        sys.exit(f"{path} not found")
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.size == 0:
        sys.exit(f"{path} has no rows")
    return np.atleast_1d(data)


def joints(d, prefix):
    return np.stack([d[f"{prefix}{i}"] for i in range(1, 8)], axis=1)


def is_position_reference(d):
    """A position reference holds eta_d = s0: the ref column is constant and eta_dot stays near zero."""
    following = d["phase"] == 0
    ref = d["ref"][following]
    return ref.size > 0 and np.ptp(ref) < 1e-12 and ref[0] > 0.0


def stop_span(d):
    stop = d["phase"] == 1
    return (d["t"][stop][0], d["t"][stop][-1]) if np.any(stop) else None


def summary(name, d):
    t = d["t"]
    f = d["phase"] == 0
    tau_frac = np.abs(joints(d, "tauG")[f]) / TAU_MAX
    margins = np.minimum(joints(d, "q") - Q_MIN, Q_MAX - joints(d, "q"))[f]
    tau = joints(d, "tau")[f]
    rate = np.abs(np.diff(np.vstack([np.zeros(7), tau]), axis=0)).max() / 1e-3 if tau.size else 0.0
    print(f"{name}: {len(t)} rows, t = {t[0]:.3f} .. {t[-1]:.3f} s")
    if not np.any(f):
        print("  no following rows")
        return
    tf = t[f]
    print(f"  max |xi| {np.max(d['xi_norm'][f]):.4f} [{XI_MAX}], max |beta| {np.max(np.abs(d['beta'][f])):.4f} rad "
          f"[{BETA_MAX}], min sigma_min(A) {np.min(d['sigma_min'][f]):.4f}")
    print(f"  joint margin {margins.min():.3f} rad (joint {np.unravel_index(margins.argmin(), margins.shape)[1] + 1})"
          f" [{JOINT_MARGIN}], |qd| {100 * (np.abs(joints(d, 'dq')[f]) / DQ_MAX).max():.1f} % [50], |qdd*| "
          f"{100 * (np.abs(joints(d, 'qdd')[f]) / DDQ_MAX).max():.1f} % [30], torque {100 * tau_frac.max():.1f} % "
          f"(joint {np.unravel_index(tau_frac.argmax(), tau_frac.shape)[1] + 1}) [60], commanded torque rate "
          f"{rate:.1f} N m/s [{TAU_RATE_MAX:.0f}]")
    last5 = f & (t >= tf[-1] - 5.0)
    print(f"  last 5 s of following: max |xi| {np.max(d['xi_norm'][last5]):.2e}, max |beta| "
          f"{np.max(np.abs(d['beta'][last5])):.2e}; final |xi| {d['xi_norm'][f][-1]:.2e}, beta "
          f"{d['beta'][f][-1]:.2e}, e_tan {d['e_tan'][f][-1]:.2e}")
    if not is_position_reference(d):
        cruise = f & (d["ref"] >= 0.999 * np.max(d["ref"])) & (np.max(d["ref"]) > 0)
        if np.any(cruise):
            print(f"  while cruising (nu_d = {np.max(d['ref']):.3f}): mean eta_dot {np.mean(d['eta_dot'][cruise]):.4f}"
                  f" ({100 * np.mean(d['eta_dot'][cruise]) / np.max(d['ref']):.0f} % of nu_d)")
    if "tauF1" in d.dtype.names:
        tau_f = np.abs(joints(d, "tauF")[f])
        print(f"  friction feedforward |tau_f| up to {tau_f.max():.3f} N m "
              f"(joint {np.unravel_index(tau_f.argmax(), tau_f.shape)[1] + 1})")
    span = stop_span(d)
    if span:
        stop = d["phase"] == 1
        q = joints(d, "q")
        drift = np.max(np.abs(q[stop] - q[stop][0])) if np.any(stop) else 0.0
        print(f"  SAFETY STOP from t = {span[0]:.3f} s to {span[1]:.3f} s; largest joint drift from the stop "
              f"configuration {drift:.4f} rad; final |dq| {np.linalg.norm(joints(d, 'dq')[-1]):.4f} rad/s")
    if np.any(d["compute_us"] > 0):
        print(f"  callback: at most {d['compute_us'].max():.0f} us; {int(np.sum(d['period_ms'][1:] > 1))} periods "
              f"longer than 1 ms; lowest success rate {d['success_rate'].min():.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("log")
    ap.add_argument("--overlay", action="append", default=[], help="PRED.csv[:label], dashed (repeatable)")
    ap.add_argument("--no-preflight", action="store_true", help="do not overlay LOG.csv.preflight.csv")
    ap.add_argument("--out", help="output file (default: LOG without .csv, plus .png and .pdf)")
    args = ap.parse_args()

    run = load(args.log)
    overlays = []
    pre = args.log + ".preflight.csv"
    if not args.no_preflight and os.path.exists(pre):
        overlays.append((load(pre), "preflight (robot's model)"))
    for item in args.overlay:
        path, _, label = item.partition(":")
        overlays.append((load(path), label or os.path.basename(path)))

    summary(os.path.basename(args.log), run)
    for d, label in overlays:
        summary(f"prediction: {label}", d)

    position = is_position_reference(run)
    fig, ax = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    ax = ax.ravel()
    t = run["t"]

    def styled(a, title, ylabel):
        a.set_title(title, color=INK, fontsize=10, loc="left")
        a.set_xlabel("t (s)", color=INK_2)
        a.set_ylabel(ylabel, color=INK_2)
        a.grid(True, color=GRID, linewidth=0.6)
        a.tick_params(colors=INK_2, labelsize=8)
        span = stop_span(run)
        if span:
            a.axvspan(span[0], span[1], color=STOP, alpha=0.6, linewidth=0, label="safety stop")

    # (a) |xi|
    a = ax[0]
    styled(a, "(a) transverse error |xi|", "|xi| (metric units)")
    a.plot(t, run["xi_norm"], color=SERIES[0], label="robot")
    for k, (d, label) in enumerate(overlays):
        a.plot(d["t"], d["xi_norm"], "--", color=SERIES[(k + 1) % len(SERIES)], label=label)
    a.axhline(XI_MAX, color=LIMIT, linestyle=":", label="bound 0.04")
    a.legend(fontsize=8)

    # (b) beta
    a = ax[1]
    styled(a, "(b) redundancy output beta", "beta (rad)")
    a.plot(t, run["beta"], color=SERIES[0], label="robot")
    for k, (d, label) in enumerate(overlays):
        a.plot(d["t"], d["beta"], "--", color=SERIES[(k + 1) % len(SERIES)], label=label)
    for s in (1, -1):
        a.axhline(s * BETA_MAX, color=LIMIT, linestyle=":", label="bound 0.05 rad" if s > 0 else None)
    a.legend(fontsize=8)

    # (c) tangential channel
    a = ax[2]
    if position:
        styled(a, "(c) tangential error eta - eta_d", "e_eta (metric units)")
        a.plot(t, run["e_tan"], color=SERIES[0], label="robot")
        for k, (d, label) in enumerate(overlays):
            a.plot(d["t"], d["e_tan"], "--", color=SERIES[(k + 1) % len(SERIES)], label=label)
    else:
        styled(a, "(c) path speed eta_dot against nu_d", "eta_dot (metric units / s)")
        a.plot(t, run["eta_dot"], color=SERIES[0], label="robot")
        for k, (d, label) in enumerate(overlays):
            a.plot(d["t"], d["eta_dot"], "--", color=SERIES[(k + 1) % len(SERIES)], label=label)
        f = run["phase"] == 0
        a.plot(t[f], run["ref"][f], ":", color=REFERENCE, label="nu_d")
    a.legend(fontsize=8)

    # (d) delivered torques
    a = ax[3]
    styled(a, "(d) delivered torque |tau_i + G_i| / tau_max,i", "fraction of the limit")
    frac = np.abs(joints(run, "tauG")) / TAU_MAX
    for i in range(7):
        a.plot(t, frac[:, i], color=SERIES[i], linewidth=1.0, label=f"joint {i + 1}")
    a.axhline(TAU_FRAC, color=LIMIT, linestyle=":", label="bound 60 %")
    a.legend(fontsize=7, ncol=2)

    # (e) joint margins
    a = ax[4]
    styled(a, "(e) distance of each joint to its nearer limit", "rad")
    q = joints(run, "q")
    margins = np.minimum(q - Q_MIN, Q_MAX - q)
    for i in range(7):
        a.plot(t, margins[:, i], color=SERIES[i], linewidth=1.0, label=f"joint {i + 1}")
    a.axhline(JOINT_MARGIN, color=LIMIT, linestyle=":", label="bound 0.10 rad")
    a.legend(fontsize=7, ncol=2)

    # (f) velocity and commanded acceleration
    a = ax[5]
    styled(a, "(f) max |qd_i|/qd_max,i and commanded |qdd*_i|/qdd_max,i", "fraction of the limit")
    a.plot(t, (np.abs(joints(run, "dq")) / DQ_MAX).max(axis=1), color=SERIES[0], label="|qd| (robot)")
    a.plot(t, (np.abs(joints(run, "qdd")) / DDQ_MAX).max(axis=1), color=SERIES[1], label="|qdd*| (robot)")
    for k, (d, label) in enumerate(overlays):
        a.plot(d["t"], (np.abs(joints(d, "qdd")) / DDQ_MAX).max(axis=1), "--", color=SERIES[(k + 2) % len(SERIES)],
               label=f"|qdd*| ({label})")
    a.axhline(DQ_FRAC, color=LIMIT, linestyle=":", label="bounds 50 % and 30 %")
    a.axhline(DDQ_FRAC, color=LIMIT, linestyle=":")
    a.legend(fontsize=7)

    fig.suptitle(os.path.basename(args.log), color=INK, fontsize=11)
    base = args.out[:-4] if args.out and args.out.endswith(".png") else (args.out or os.path.splitext(args.log)[0])
    fig.savefig(base + ".png", dpi=130)
    fig.savefig(base + ".pdf")
    print(f"wrote {base}.png and {base}.pdf")


if __name__ == "__main__":
    main()
