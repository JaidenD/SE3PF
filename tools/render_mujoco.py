#!/usr/bin/env python3
"""Render a video of a MuJoCo run of the Panda example from the logs of apps/panda_mujoco (docs/MUJOCO.md).

    MUJOCO_GL=egl .venv/bin/python tools/render_mujoco.py [--run C] [--dir out] [--out out/mujoco_panda.mp4]
                                                          [--fps 30] [--size 1280x720] [--stills 0,1,5]

Replays the joint angles q(t) of <dir>/mujoco_<run>.csv in the Menagerie scene.xml (the arm with the Franka Hand,
fingers open, floor and lights) and draws, as user geoms of the MuJoCo scene:
  - the path gamma of the flange (<dir>/mujoco_path.csv) as a thin polyline,
  - the target frame gamma(eta(t)) as a small triad (x red, y green, z blue),
  - the flange's own trail over the last 2 s,
  - a text overlay with t, |xi| and beta.
The hand is drawn in every run; only run D's plant carries it (run C is the arm alone, so the picture of C shows the
hand that run D adds). The flange is the hand's mounting plate: the path runs through the base of the hand.
--stills also writes PNG frames at the given times (s) next to the video, to check the camera and the overlays.
Needs mujoco (same version as the C library), numpy, imageio with imageio-ffmpeg, and Pillow (with matplotlib).
view_mujoco.py imports the scene helpers of this file.
"""
import argparse
import os
import sys

import numpy as np
import mujoco

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SCENE = os.path.join(ROOT, "third_party/mujoco_menagerie/franka_emika_panda/scene.xml")
FLANGE_OFFSET = 0.107  # the flange (Franka frame F) is link7 moved 0.107 m along z (PandaArm::kFlangeOffset)

RUN_TITLES = {
    "A": "A: exact model",
    "A2": "A2: self-motion, exact model",
    "B": "B: zero-order hold (1 ms)",
    "C": "C: hardware-like (hold, ramp, damping)",
    "C0": "C0: hardware-like, no ramp",
    "D": "D: Franka Hand not modeled",
    "DG": "D-G: hand only in G",
}

# Colors (RGBA, 0..1). The path in warm white, readable on the scene's dark blue background; the trail in the
# run color of tools/plot_mujoco.py.
PATH_RGBA = np.array([0.98, 0.96, 0.88, 0.95], dtype=np.float32)
TRAIL_RGBA = {"A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a", "D": "#eda100", "C0": "#e87ba4"}
AXIS_RGBA = [np.array(c, dtype=np.float32) for c in ([0.90, 0.20, 0.15, 1.0], [0.20, 0.80, 0.25, 1.0],
                                                      [0.20, 0.45, 0.95, 1.0])]


def hex_rgba(color, alpha=1.0):
    color = color.lstrip("#")
    return np.array([int(color[i:i + 2], 16) / 255.0 for i in (0, 2, 4)] + [alpha], dtype=np.float32)


# ---------------------------------------------------------------------------------------------------------------
# Logs

def load_run(directory, run):
    """Columns of <directory>/mujoco_<run>.csv as a dict of arrays."""
    path = os.path.join(directory, f"mujoco_{run}.csv")
    if not os.path.exists(path):
        sys.exit(f"{path} not found; run  ./build-rel/panda_mujoco {directory}  first")
    with open(path) as f:
        names = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return {name: data[:, i] for i, name in enumerate(names)}


class Path:
    """The path gamma(s), sampled by panda_mujoco at s_k = k L / n (columns s, p, R row by row)."""

    def __init__(self, directory):
        path = os.path.join(directory, "mujoco_path.csv")
        if not os.path.exists(path):
            sys.exit(f"{path} not found; run  ./build-rel/panda_mujoco {directory}  first")
        data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
        self.s = data[:, 0]
        self.p = data[:, 1:4]
        self.R = data[:, 4:13].reshape(-1, 3, 3)
        self.length = self.s[-1] + (self.s[1] - self.s[0])

    def pose(self, eta):
        """gamma(eta) at the nearest sample (2000 samples: 0.4 mm of arc length apart)."""
        k = int(round((eta % self.length) / self.length * len(self.s))) % len(self.s)
        return self.p[k], self.R[k]


# ---------------------------------------------------------------------------------------------------------------
# Scene helpers: user geoms appended to an mjvScene after mjv_updateScene (or a viewer's user_scn)

def add_geom(scene, geom_type, size, pos, mat, rgba):
    """Appends one geom; silently skipped when the scene is full."""
    if scene.ngeom >= scene.maxgeom:
        return None
    g = scene.geoms[scene.ngeom]
    mujoco.mjv_initGeom(g, geom_type, np.asarray(size, dtype=float), np.asarray(pos, dtype=float),
                        np.asarray(mat, dtype=float).reshape(9), rgba)
    g.category = mujoco.mjtCatBit.mjCAT_DECOR  # decoration: not selectable, no shadow casting needed
    scene.ngeom += 1
    return g


def add_segment(scene, a, b, radius, rgba):
    g = add_geom(scene, mujoco.mjtGeom.mjGEOM_CAPSULE, [radius, 0, 0], np.zeros(3), np.eye(3), rgba)
    if g is not None:
        mujoco.mjv_connector(g, mujoco.mjtGeom.mjGEOM_CAPSULE, radius, np.asarray(a, float), np.asarray(b, float))


def add_polyline(scene, points, radius, rgba, closed=False, stride=1):
    pts = points[::stride]
    for i in range(len(pts) - 1):
        add_segment(scene, pts[i], pts[i + 1], radius, rgba)
    if closed:
        add_segment(scene, pts[-1], pts[0], radius, rgba)


def add_triad(scene, p, R, length=0.08, radius=0.0025):
    for axis in range(3):
        add_segment(scene, p, p + length * R[:, axis], radius, AXIS_RGBA[axis])
    add_geom(scene, mujoco.mjtGeom.mjGEOM_SPHERE, [1.6 * radius, 0, 0], p, np.eye(3),
             np.array([0.95, 0.95, 0.95, 1.0], dtype=np.float32))


def flange_position(model, data):
    link7 = model.body("link7").id
    R = data.xmat[link7].reshape(3, 3)
    return data.xpos[link7] + FLANGE_OFFSET * R[:, 2]


def flange_trail(model, data, log, k, seconds=2.0, every=20):
    """Flange positions at the rows k - seconds/dt .. k (every `every` rows), by forward kinematics of the log's q."""
    dt = log["t"][1] - log["t"][0]
    first = max(0, k - int(round(seconds / dt)))
    q_now = data.qpos[:7].copy()
    points = []
    for j in range(first, k + 1, every):
        data.qpos[:7] = [log[f"q{i + 1}"][j] for i in range(7)]
        mujoco.mj_kinematics(model, data)
        points.append(flange_position(model, data))
    data.qpos[:7] = q_now
    mujoco.mj_kinematics(model, data)
    return np.array(points)


def set_arm(model, data, log, k):
    data.qpos[:7] = [log[f"q{i + 1}"][k] for i in range(7)]
    if model.nq > 7:
        data.qpos[7:9] = 0.04  # fingers open
    mujoco.mj_forward(model, data)


def decorate(scene, model, data, log, k, path, run):
    """Path, target triad gamma(eta) and the flange trail, for row k of the log."""
    add_polyline(scene, path.p, 0.0012, PATH_RGBA, closed=True, stride=4)
    p, R = path.pose(log["eta"][k])
    add_triad(scene, p, R)
    trail = flange_trail(model, data, log, k)
    if len(trail) > 1:
        add_polyline(scene, trail, 0.002, hex_rgba(TRAIL_RGBA.get(run, "#eda100")))


def overlay_text(log, k):
    return ("t = %5.2f s   |xi| = %.2e   beta = %+.2e rad   eta_dot = %.4f" %
            (log["t"][k], log["xi_norm"][k], log["beta"][k], log["eta_dot"][k]))


LOOKAT = [0.30, 0.0, 0.50]


def camera(azimuth=125.0, elevation=-32.0, distance=1.1):
    """A free camera on the coning circle (center (0.40, 0, 0.45) m, radius 0.12 m), from the front left and above,
    so that both the circle and the tilt of the approach axis (the cone) are visible."""
    cam = mujoco.MjvCamera()
    cam.type = mujoco.mjtCamera.mjCAMERA_FREE
    cam.lookat[:] = LOOKAT
    cam.distance = distance
    cam.azimuth = azimuth
    cam.elevation = elevation
    return cam


def translucent_hand(model, alpha=0.3):
    """Draw the Franka Hand (bodies hand, left_finger, right_finger) translucent, so the path and the target frame
    at the flange stay visible through it. Only the drawing changes: masses and collisions are untouched."""
    hand_bodies = {model.body(name).id for name in ("hand", "left_finger", "right_finger")
                   if mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, name) >= 0}
    for g in range(model.ngeom):
        if model.geom_bodyid[g] in hand_bodies:
            rgba = model.mat_rgba[model.geom_matid[g]].copy() if model.geom_matid[g] >= 0 else model.geom_rgba[g]
            model.geom_matid[g] = -1  # with a material, MuJoCo would draw the material's color
            model.geom_rgba[g] = [rgba[0], rgba[1], rgba[2], alpha]


# ---------------------------------------------------------------------------------------------------------------

def draw_text(frame, lines):
    """Text box in the top-left corner (Pillow)."""
    from PIL import Image, ImageDraw, ImageFont
    from matplotlib import font_manager

    image = Image.fromarray(frame)
    draw = ImageDraw.Draw(image, "RGBA")
    size = max(14, frame.shape[0] // 40)
    font = ImageFont.truetype(font_manager.findfont("DejaVu Sans Mono"), size)
    pad, step = size // 2, int(size * 1.35)
    width = max(draw.textlength(line, font=font) for line in lines)
    draw.rectangle([0, 0, width + 2 * pad, len(lines) * step + 2 * pad], fill=(10, 12, 20, 170))
    for i, line in enumerate(lines):
        draw.text((pad, pad + i * step), line, font=font, fill=(240, 240, 235, 255))
    return np.asarray(image)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", default="C", help="run name: A, A2, B, C, C0, D or DG (default C)")
    parser.add_argument("--dir", default=os.path.join(ROOT, "out"), help="directory of the logs (default out)")
    parser.add_argument("--out", help="video file (default <dir>/mujoco_panda.mp4)")
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument("--size", default="1280x720", help="width x height (default 1280x720)")
    parser.add_argument("--stills", default="", help="comma-separated times (s) to save as PNG as well")
    parser.add_argument("--scene", default=SCENE)
    args = parser.parse_args()

    import imageio.v2 as imageio

    width, height = (int(v) for v in args.size.lower().split("x"))
    out = args.out or os.path.join(args.dir, "mujoco_panda.mp4")
    log = load_run(args.dir, args.run)
    path = Path(args.dir)

    model = mujoco.MjModel.from_xml_path(args.scene)
    model.vis.global_.offwidth = max(model.vis.global_.offwidth, width)    # the off-screen buffer must be large
    model.vis.global_.offheight = max(model.vis.global_.offheight, height)  # enough for the frame
    translucent_hand(model)
    data = mujoco.MjData(model)
    cam = camera()
    title = "MuJoCo Panda, run " + RUN_TITLES.get(args.run, args.run)

    dt = log["t"][1] - log["t"][0]
    duration = log["t"][-1]
    frames = int(np.floor(duration * args.fps + 1e-9)) + 1
    stills = [float(s) for s in args.stills.split(",") if s.strip()]

    with mujoco.Renderer(model, height=height, width=width, max_geom=4000) as renderer:
        def render(k):
            set_arm(model, data, log, k)
            renderer.update_scene(data, camera=cam)
            decorate(renderer.scene, model, data, log, k, path, args.run)
            return draw_text(renderer.render(), [title, overlay_text(log, k)])

        writer = imageio.get_writer(out, fps=args.fps, codec="libx264", quality=8, macro_block_size=8)
        for f in range(frames):
            k = min(len(log["t"]) - 1, int(round(f / args.fps / dt)))
            writer.append_data(render(k))
        writer.close()
        print("wrote %s (%d frames, %.1f s at %d fps, %dx%d, %.1f MB)" % (
            out, frames, frames / args.fps, args.fps, width, height, os.path.getsize(out) / 1e6))

        for t in stills:
            k = min(len(log["t"]) - 1, int(round(t / dt)))
            png = os.path.splitext(out)[0] + "_t%05.2f.png" % t
            imageio.imwrite(png, render(k))
            print("wrote", png)


if __name__ == "__main__":
    main()
