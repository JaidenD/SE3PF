#!/usr/bin/env python3
"""Plot a robot run: the desired and actual end-effector path in 3D, the tracking error, and the elbow.

    python3 tools/plot_robot_path.py LOG.csv [--out FILE.png] [--frames N]

LOG.csv is a log of robot_driver (docs/HARDWARE.md). The flange pose is forward kinematics of the logged q (the
Panda's modified DH, as the library's PandaArm); the desired pose at each instant is gamma(eta), the point of the
coning circle the controller projected the flange onto (so the gap between the two is |xi| made visible: a position
offset and an attitude offset). The elbow (the origin of joint 4) is compared with its position on the lift q*(eta),
read from the preflight prediction LOG.csv.preflight.csv: the gap between the two is what beta measures (the
self-motion: the elbow swinging while the flange stays). Panels:
  (a) 3D: the desired path gamma(s) and the actual end-effector (flange) path, with --frames body frames on each
      (at the same instants, evenly over one loop: desired dashed, actual solid; x red, y green, z blue);
  (b) the transverse error: the end effector's position error (mm) and attitude error (deg) against the path, with
      |xi| for comparison;
  (c) the tangential output against its reference: eta_dot and nu_d for a velocity reference, eta and eta_d for a
      position reference (eta is arc length in the metric W = diag(0.01 I, I), so its "mm" count 1 mrad of
      rotation as 0.1 mm);
  (d) the elbow's distance from its lift position (mm), with beta for comparison.
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


def has_gamma(d):
    """Logs since 2026-11-02 carry the desired pose gamma(eta) (gx, gy, gz, grx, gry, grz)."""
    return "gx" in d.dtype.names and np.isfinite(d["gx"]).any()


def desired_pose(d, rows=slice(None)):
    """gamma(eta) of the given rows: from the logged columns if present, else the coning circle (older logs)."""
    if has_gamma(d):
        p = np.stack([d["gx"][rows], d["gy"][rows], d["gz"][rows]], -1)
        w = np.stack([d["grx"][rows], d["gry"][rows], d["grz"][rows]], -1)
        return np.array([expSO3(x) if np.all(np.isfinite(x)) else np.full((3, 3), np.nan) for x in w]), p
    return path_pose(d["eta"][rows])


def path_length(d):
    """L of the run's path: the coning circle's for older logs; else the largest eta of a run that wrapped past L,
    or None if eta never wrapped (then nothing needs unwrapping)."""
    if not has_gamma(d):
        return L
    eta = d["eta"][np.isfinite(d["eta"])]
    return float(eta.max()) if np.any(np.diff(eta) < -0.5 * np.ptp(eta)) else None


def unwrap(eta, length):
    """eta without its jumps from L back to 0."""
    return eta if length is None else np.unwrap(2.0 * np.pi * eta / length) * length / (2.0 * np.pi)


def is_position_reference(d):
    """A velocity reference logs e_tan = eta_dot - nu_d; a position reference logs e_tan = eta - eta_d."""
    f = d["phase"] == 0
    return f.any() and not np.allclose(d["e_tan"][f], (d["eta_dot"] - d["ref"])[f], atol=1e-6)


def elbow(q):
    """The origin of joint 4's frame (the elbow)."""
    return fk(q, upTo=4)[1]


def joints(d, prefix):
    return np.stack([d[f"{prefix}{i}"] for i in range(1, 8)], axis=1)


AXIS_COLORS = ("#d62728", "#2ca02c", "#1f77b4")  # body x, y, z


def triad(ax, R, p, size, ls, lw):
    """The body frame (R, p): its x, y and z axes as red, green and blue segments of length size."""
    for k, color in enumerate(AXIS_COLORS):
        e = p + size * R[:, k]
        ax.plot([p[0], e[0]], [p[1], e[1]], [p[2], e[2]], color=color, linestyle=ls, linewidth=lw)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("log")
    ap.add_argument("--out")
    ap.add_argument("--frames", type=int, default=12, help="body frames drawn on each path in (a)")
    args = ap.parse_args()
    d = np.atleast_1d(np.genfromtxt(args.log, delimiter=",", names=True))
    f = d["phase"] == 0
    t, eta = d["t"][f], d["eta"][f]
    q = joints(d, "q")[f]
    R, p = fk(q)
    Rd, pd = desired_pose(d, f)
    Lp = path_length(d)
    posErr = np.linalg.norm(p - pd, axis=1)
    attErr = np.linalg.norm(logSO3(np.swapaxes(R, -1, -2) @ Rd), axis=1)
    el = elbow(q)
    # The elbow on the lift at the same eta, from the preflight (its q is q*(eta) wherever its beta is zero)
    pre = args.log + ".preflight.csv"
    elLift = None
    if os.path.exists(pre):
        pp = np.atleast_1d(np.genfromtxt(pre, delimiter=",", names=True))
        # Only the prediction's samples on the lift (|beta| < 1e-4): a run that starts off the lift (self-motion)
        # predicts the same transient, which is not q*
        onLift = (pp["phase"] == 0) & (np.abs(pp["beta"]) < 1e-4)
        pe, pq = pp["eta"][onLift], joints(pp, "q")[onLift]
        first = np.ones(len(pe), dtype=bool)
        order = np.argsort(pe[first])
        Lw = Lp if Lp is not None else path_length(pp)
        wrapped = eta if Lw is None else eta % Lw
        qLift = np.stack([np.interp(wrapped, pe[first][order], pq[first][order][:, i]) for i in range(7)], axis=1)
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

    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    # The desired path: the preflight's gamma (the whole planned motion), else this run's, else the coning circle
    if os.path.exists(pre) and has_gamma(pp):
        pc = desired_pose(pp, pp["phase"] == 0)[1]
    elif has_gamma(d):
        pc = pd
    else:
        pc = path_pose(np.linspace(0.0, L, 400))[1]
    ax.plot(pc[:, 0], pc[:, 1], pc[:, 2], color=INK, linestyle="--", linewidth=1.4, label="desired")
    ax.plot(p[:, 0], p[:, 1], p[:, 2], color=SERIES[0], linewidth=1.0, label="actual")
    # Body frames at the same instants: the desired frame gamma(eta(t)) (dashed) and the actual end-effector frame
    # (solid), so the angle between each pair is the attitude error. The instants are evenly spaced in arc length
    # over one loop centred on the middle of the run (where a run covers less than a loop, over what it covers)
    from matplotlib.lines import Line2D
    u = unwrap(eta, Lp)
    span = (u[-1] - u[0]) if Lp is None else min(Lp, u[-1] - u[0])
    mid = 0.5 * (u[0] + u[-1])
    targets = mid - 0.5 * span + span * (np.arange(args.frames) + 0.5) / args.frames
    for k in np.searchsorted(np.maximum.accumulate(u), targets).clip(0, len(t) - 1):
        triad(ax, Rd[k], pd[k], 0.03, "--", 1.0)
        triad(ax, R[k], p[k], 0.03, "-", 1.4)
    ax.set_xlabel("$x$ (m)"); ax.set_ylabel("$y$ (m)"); ax.set_zlabel("$z$ (m)")
    ax.set_title("(a) End-effector path and body frames", loc="left")
    handles = [Line2D([], [], color=INK, linestyle="--", linewidth=1.4, label="desired"),
               Line2D([], [], color=SERIES[0], linewidth=1.0, label="actual")]
    handles += [Line2D([], [], color=c, linewidth=1.4, label=f"body ${a}$ axis") for c, a in zip(AXIS_COLORS, "xyz")]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=8)
    # Equal scale on all three axes, fitted to the two paths
    allp = np.vstack([pc, p])
    c = 0.5 * (allp.min(axis=0) + allp.max(axis=0)); rad = 0.55 * np.max(np.ptp(allp, axis=0))
    ax.set_xlim(c[0] - rad, c[0] + rad); ax.set_ylim(c[1] - rad, c[1] + rad); ax.set_zlim(c[2] - rad, c[2] + rad)
    ax.set_box_aspect((1, 1, 1))

    axb = fig.add_subplot(3, 2, 2)
    axb.plot(t, 1e3 * posErr, color=SERIES[0], label="position (mm)")
    axb.plot(t, np.degrees(attErr), color=SERIES[2], label="attitude (deg)")
    axb.plot(t, 1e3 * d["xi_norm"][f], color=INK_2, linestyle=":", label=r"$|\xi| \times 10^3$")
    axb.set_title("(b) Transverse tracking error", loc="left")
    axb.grid(True, color=GRID); axb.legend(fontsize=8, frameon=False); axb.set_xlabel("$t$ (s)")

    # (c) Tangential output against its reference
    axt = fig.add_subplot(3, 2, 4, sharex=axb)
    ref = d["ref"][f]
    if is_position_reference(d):
        # eta wraps at L, eta_d does not: unwrap eta and start it on eta_d's side
        ue = unwrap(eta, Lp)
        if Lp is not None:
            ue = ue + Lp * np.round((ref[0] - ue[0]) / Lp)
        axt.plot(t, 1e3 * ue, color=SERIES[0], linewidth=1.0, label=r"actual $\eta$")
        axt.plot(t, 1e3 * ref, color=INK, linestyle="--", linewidth=1.2, label=r"desired $\eta_d$")
        axt.set_ylabel("mm")
    else:
        axt.plot(t, 1e3 * d["eta_dot"][f], color=SERIES[0], linewidth=0.6, label=r"actual $\dot\eta$")
        axt.plot(t, 1e3 * ref, color=INK, linestyle="--", linewidth=1.2, label=r"desired $\nu_d$")
        axt.set_ylabel("mm/s")
    axt.set_title("(c) Tangential tracking", loc="left")
    axt.grid(True, color=GRID); axt.legend(fontsize=8, frameon=False, loc="lower center"); axt.set_xlabel("$t$ (s)")

    axc = fig.add_subplot(3, 2, 6, sharex=axb)
    if elLift is not None:
        axc.plot(t, 1e3 * np.linalg.norm(el - elLift, axis=1), color=SERIES[1], label="elbow deviation (mm)")
    axc.plot(t, 1e3 * np.abs(d["beta"][f]), color=INK_2, linestyle=":", label=r"$|\beta|$ (mrad)")
    axc.set_title("(d) Elbow deviation from the lift", loc="left")
    axc.grid(True, color=GRID); axc.legend(fontsize=8, frameon=False); axc.set_xlabel("$t$ (s)")
    fig.tight_layout()
    base = args.out[:-4] if args.out and args.out.endswith(".png") else (args.out or os.path.splitext(args.log)[0] + "_path")
    fig.savefig(base + ".png", dpi=200)
    fig.savefig(base + ".pdf")
    print(f"wrote {base}.png and {base}.pdf")


if __name__ == "__main__":
    main()
