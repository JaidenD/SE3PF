#!/usr/bin/env python3
"""Plot a robot run in 3D: the desired path with its frames, the flange's realized path and frames, and the elbow.

    python3 tools/plot_robot_path.py LOG.csv [--out FILE.png] [--every SECONDS]

LOG.csv is a log of robot_driver (docs/HARDWARE.md). The flange pose is forward kinematics of the logged q (the
Panda's modified DH, as the library's PandaArm); the desired pose at each instant is gamma(eta), the point of the
coning circle the controller projected the flange onto (so the gap between the two is |xi| made visible: a position
offset and an attitude offset). The elbow (the origin of joint 4) is compared with its position on the lift q*(eta),
read from the preflight prediction LOG.csv.preflight.csv: the gap between the two is what beta measures (the
self-motion: the elbow swinging while the flange stays). Panels:
  (a) 3D: the coning circle (grey), its frames every --every seconds of the reference (dashed triads), the realized
      flange path coloured by time with its frames (solid triads), the realized elbow path and the elbow on the lift;
  (b) the flange's position error (mm) and attitude error (deg) against the path, with |xi| for comparison;
  (c) the elbow's distance from its lift position (mm), with beta for comparison.
Needs numpy and matplotlib.
"""
import argparse
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# The Panda's kinematics (modified DH of PandaArm, src/panda.cpp; the same numbers as docs/closed_lifts/work/num/
# pandakin.py, repeated here so that this tool needs no data files)
DH = np.array([[0.0, 0.333, 0.0], [0.0, 0.0, -np.pi / 2], [0.0, 0.316, np.pi / 2], [0.0825, 0.0, np.pi / 2],
               [-0.0825, 0.384, -np.pi / 2], [0.0, 0.0, np.pi / 2], [0.088, 0.0, np.pi / 2]])
FLANGE = 0.107
QR = np.array([0.0, -np.pi / 4, 0.0, -3 * np.pi / 4, 0.0, np.pi / 2, np.pi / 4])


def hat(w):
    return np.array([[0.0, -w[2], w[1]], [w[2], 0.0, -w[0]], [-w[1], w[0], 0.0]])


def expSO3(w):
    th = np.linalg.norm(w)
    if th < 1e-12:
        return np.eye(3)
    K = hat(w / th)
    return np.eye(3) + np.sin(th) * K + (1.0 - np.cos(th)) * (K @ K)


def logSO3(R):
    """Batched log on SO(3) for angles below pi, as a vector of rotation angles."""
    tr = np.trace(R, axis1=-2, axis2=-1)
    c = np.clip((tr - 1.0) / 2.0, -1.0, 1.0)
    th = np.arccos(c)
    v = np.stack([R[..., 2, 1] - R[..., 1, 2], R[..., 0, 2] - R[..., 2, 0], R[..., 1, 0] - R[..., 0, 1]], -1)
    small = th < 1e-6
    s = np.where(small, 1.0, np.sin(th))
    fac = np.where(small, 0.5 + th ** 2 / 12.0, th / (2.0 * s))
    return fac[..., None] * v


def tdh(a, d, al, th):
    ct, st = np.cos(th), np.sin(th)
    ca, sa = np.cos(al), np.sin(al)
    T = np.zeros(th.shape + (4, 4))
    T[..., 0, 0] = ct; T[..., 0, 1] = -st; T[..., 0, 3] = a
    T[..., 1, 0] = st * ca; T[..., 1, 1] = ct * ca; T[..., 1, 2] = -sa; T[..., 1, 3] = -d * sa
    T[..., 2, 0] = st * sa; T[..., 2, 1] = ct * sa; T[..., 2, 2] = ca; T[..., 2, 3] = d * ca
    T[..., 3, 3] = 1.0
    return T


def fk(q, upTo=7):
    """The frame after joint upTo (the flange for 7: joint 7's frame moved by FLANGE along z): R (B, 3, 3), p (B, 3)."""
    q = np.atleast_2d(np.asarray(q, dtype=float))
    T = np.tile(np.eye(4), (q.shape[0], 1, 1))
    for i in range(upTo):
        T = T @ tdh(DH[i, 0], DH[i, 1], DH[i, 2], q[:, i])
    R = T[:, :3, :3]
    p = T[:, :3, 3] + (FLANGE * T[:, :3, 2] if upTo == 7 else 0.0)
    return R, p

# The paper's path (src/panda_example.cpp, src/path.cpp): the coning circle around p0 with radius r and cone angle
# alpha, R0 the flange attitude at the rest posture, s = sigma theta with L = 2 pi sigma
CENTER = np.array([0.40, 0.0, 0.45])
RADIUS, CONE, L = 0.12, 0.25, 0.77008765
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#7a5cc9", "#4c4b47"]
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def rz(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def path_pose(eta):
    """gamma(eta): attitude R0 Rz(th) A Rz(-th), position c + r (cos th, sin th, 0), th = 2 pi eta / L."""
    R0 = fk(QR)[0][0]
    A = expSO3(np.array([CONE, 0.0, 0.0]))
    th = 2.0 * np.pi * np.atleast_1d(np.asarray(eta, dtype=float)) / L
    R = np.array([R0 @ rz(t) @ A @ rz(-t) for t in th])
    p = CENTER + RADIUS * np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], -1)
    return R, p


def elbow(q):
    """The origin of joint 4's frame (the elbow)."""
    return fk(q, upTo=4)[1]


def joints(d, prefix):
    return np.stack([d[f"{prefix}{i}"] for i in range(1, 8)], axis=1)


def triad(ax, R, p, size, ls, alpha=1.0, lw=1.2):
    for k, color in enumerate(("#d62728", "#2ca02c", "#1f77b4")):  # x red, y green, z blue
        ax.plot([p[0], p[0] + size * R[0, k]], [p[1], p[1] + size * R[1, k]], [p[2], p[2] + size * R[2, k]],
                color=color, linestyle=ls, linewidth=lw, alpha=alpha)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("log")
    ap.add_argument("--out")
    ap.add_argument("--every", type=float, default=2.0, help="seconds between frames in the 3D panel")
    args = ap.parse_args()
    d = np.atleast_1d(np.genfromtxt(args.log, delimiter=",", names=True))
    f = d["phase"] == 0
    t, eta = d["t"][f], d["eta"][f]
    q = joints(d, "q")[f]
    R, p = fk(q)
    Rd, pd = path_pose(eta)
    posErr = np.linalg.norm(p - pd, axis=1)
    attErr = np.linalg.norm(logSO3(np.swapaxes(R, -1, -2) @ Rd), axis=1)
    el = elbow(q)
    # The elbow on the lift at the same eta, from the preflight (its q is q*(eta) on the first loop)
    pre = args.log + ".preflight.csv"
    elLift = None
    if os.path.exists(pre):
        pp = np.atleast_1d(np.genfromtxt(pre, delimiter=",", names=True))
        pe, pq = pp["eta"], joints(pp, "q")
        first = np.arange(len(pe)) < len(pe) // 2
        order = np.argsort(pe[first])
        qLift = np.stack([np.interp(eta % L, pe[first][order], pq[first][order][:, i]) for i in range(7)], axis=1)
        elLift = elbow(qLift)
    print(f"{os.path.basename(args.log)}: {len(t)} following rows, t = {t[0]:.2f} .. {t[-1]:.2f} s, "
          f"eta {eta[0]:.3f} -> {eta[-1]:.3f}")
    print(f"  flange position error: mean {1e3 * posErr.mean():.1f} mm, max {1e3 * posErr.max():.1f} mm; "
          f"attitude error: mean {np.degrees(attErr).mean():.2f} deg, max {np.degrees(attErr).max():.2f} deg "
          f"(|xi| max {d['xi_norm'][f].max():.4f}: 0.01 is about 1 cm or 5.7 deg)")
    if elLift is not None:
        gap = np.linalg.norm(el - elLift, axis=1)
        print(f"  elbow off its lift position: mean {1e3 * gap.mean():.1f} mm, max {1e3 * gap.max():.1f} mm "
              f"(|beta| max {np.abs(d['beta'][f]).max():.4f} rad)")

    fig = plt.figure(figsize=(16, 7))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    th = np.linspace(0.0, L, 400)
    Rc, pc = path_pose(th)
    ax.plot(pc[:, 0], pc[:, 1], pc[:, 2], color=INK_2, linewidth=1.0, label="desired path gamma(s)")
    for s in np.linspace(0.0, L, 12, endpoint=False):
        Rs, ps = path_pose(s)
        triad(ax, Rs[0], ps[0], 0.03, "--", alpha=0.6, lw=0.9)
    step = max(1, int(round(args.every / max(np.median(np.diff(t)), 1e-3))))
    sc = ax.scatter(p[::5, 0], p[::5, 1], p[::5, 2], c=t[::5], cmap="viridis", s=3, label="realized flange")
    for k in range(0, len(t), step):
        triad(ax, R[k], p[k], 0.03, "-", lw=1.2)
    ax.plot(el[:, 0], el[:, 1], el[:, 2], color=SERIES[1], linewidth=1.2, label="realized elbow (joint 4)")
    if elLift is not None:
        ax.plot(elLift[:, 0], elLift[:, 1], elLift[:, 2], color=SERIES[1], linewidth=1.0, linestyle="--",
                label="elbow on the lift q*(eta)")
    ax.plot([0], [0], [0], marker="s", color=INK, markersize=5, linestyle="none", label="base")
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
    ax.set_title("(a) desired and realized flange path and frames (x red, y green, z blue; dashed: reference)",
                 fontsize=9, loc="left")
    ax.legend(fontsize=7, loc="upper left")
    allp = np.vstack([pc, p, el])
    c = allp.mean(axis=0); rad = 0.55 * np.max(np.ptp(allp, axis=0))
    ax.set_xlim(c[0] - rad, c[0] + rad); ax.set_ylim(c[1] - rad, c[1] + rad); ax.set_zlim(max(0.0, c[2] - rad), c[2] + rad)
    fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.1, label="t (s)")

    axb = fig.add_subplot(2, 2, 2)
    axb.plot(t, 1e3 * posErr, color=SERIES[0], label="flange position error (mm)")
    axb.plot(t, np.degrees(attErr), color=SERIES[2], label="flange attitude error (deg)")
    axb.plot(t, 1e3 * d["xi_norm"][f], color=INK_2, linestyle=":", label="|xi| x 1000 (1 unit = 1 cm or 5.7 deg)")
    axb.set_title("(b) the flange against the path: what |xi| means", fontsize=9, loc="left")
    axb.grid(True, color=GRID); axb.legend(fontsize=7); axb.set_xlabel("t (s)")

    axc = fig.add_subplot(2, 2, 4)
    if elLift is not None:
        axc.plot(t, 1e3 * np.linalg.norm(el - elLift, axis=1), color=SERIES[1], label="elbow off its lift position (mm)")
    axc.plot(t, 1e3 * np.abs(d["beta"][f]), color=INK_2, linestyle=":", label="|beta| x 1000 (mrad of self-motion)")
    axc.set_title("(c) the elbow against the lift: what beta means", fontsize=9, loc="left")
    axc.grid(True, color=GRID); axc.legend(fontsize=7); axc.set_xlabel("t (s)")
    fig.suptitle(os.path.basename(args.log), color=INK, fontsize=11)
    fig.tight_layout()
    out = args.out or os.path.splitext(args.log)[0] + "_path.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
