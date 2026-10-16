#!/usr/bin/env python3
"""Live replay of a MuJoCo run of the Panda example in MuJoCo's interactive viewer (docs/MUJOCO.md).

    .venv/bin/python tools/view_mujoco.py [--run C] [--dir out] [--speed 1.0] [--loop]

Plays the joint angles of <dir>/mujoco_<run>.csv in real time (or --speed times faster) in the Menagerie scene.xml,
with the overlays of tools/render_mujoco.py: the path gamma, the target frame gamma(eta(t)), the flange's trail and
the outputs t, |xi|, beta in the top-left corner. The camera starts where the video's is; drag to move it as usual.
Keys: space pauses and resumes, backspace restarts. Close the window to quit.
--dry-run builds every frame's overlay without opening a window (a check that needs no display).
Needs a display (GLFW); it uses the same Python packages as render_mujoco.py.
"""
import argparse
import os
import sys
import time

import mujoco
import mujoco.viewer  # (at the top: an import inside main would make the name mujoco local to main)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode = True  # importing render_mujoco would otherwise leave tools/__pycache__/ in the tree
from render_mujoco import (ROOT, RUN_TITLES, SCENE, Path, camera, decorate, load_run,  # noqa: E402
                           set_arm, translucent_hand)

KEY_SPACE = 32       # GLFW key codes, as the viewer's key_callback receives them
KEY_BACKSPACE = 259


def texts(log, k, run):
    """Two columns for the viewer's text overlay: names and values."""
    names = "run\nt\n|xi|\nbeta\neta_dot"
    values = "%s\n%.2f s\n%.2e\n%+.2e rad\n%.4f" % (RUN_TITLES.get(run, run), log["t"][k], log["xi_norm"][k],
                                                   log["beta"][k], log["eta_dot"][k])
    return names, values


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", default="C", help="run name: A, A2, B, C, C0, D or DG (default C)")
    parser.add_argument("--dir", default=os.path.join(ROOT, "out"), help="directory of the logs (default out)")
    parser.add_argument("--speed", type=float, default=1.0, help="replay speed (default 1: real time)")
    parser.add_argument("--loop", action="store_true", help="start again at the end")
    parser.add_argument("--scene", default=SCENE)
    parser.add_argument("--dry-run", action="store_true", help="build the overlays of every frame, no window")
    args = parser.parse_args()

    log = load_run(args.dir, args.run)
    path = Path(args.dir)
    model = mujoco.MjModel.from_xml_path(args.scene)
    translucent_hand(model)
    data = mujoco.MjData(model)
    dt = log["t"][1] - log["t"][0]
    last_row = len(log["t"]) - 1

    if args.dry_run:
        scene = mujoco.MjvScene(model, maxgeom=4000)
        frames = 0
        for k in range(0, last_row + 1, int(round(1.0 / 30.0 / dt))):
            set_arm(model, data, log, k)
            scene.ngeom = 0
            decorate(scene, model, data, log, k, path, args.run)
            texts(log, k, args.run)
            frames += 1
        print("dry run: %d frames, %d user geoms in the last one" % (frames, scene.ngeom))
        return

    state = {"paused": False, "restart": False}

    def on_key(key):
        if key == KEY_SPACE:
            state["paused"] = not state["paused"]
        elif key == KEY_BACKSPACE:
            state["restart"] = True

    with mujoco.viewer.launch_passive(model, data, key_callback=on_key) as viewer:
        start = camera()
        viewer.cam.lookat[:] = start.lookat
        viewer.cam.distance = start.distance
        viewer.cam.azimuth = start.azimuth
        viewer.cam.elevation = start.elevation

        t_log = 0.0
        last = time.perf_counter()
        while viewer.is_running():
            now = time.perf_counter()
            if not state["paused"]:
                t_log += (now - last) * args.speed
            last = now
            if state["restart"]:
                t_log, state["restart"] = 0.0, False
            if t_log > log["t"][-1]:
                t_log = 0.0 if args.loop else log["t"][-1]
            k = min(last_row, int(round(t_log / dt)))

            with viewer.lock():
                set_arm(model, data, log, k)
                viewer.user_scn.ngeom = 0
                decorate(viewer.user_scn, model, data, log, k, path, args.run)
            names, values = texts(log, k, args.run)
            viewer.set_texts((None, mujoco.mjtGridPos.mjGRID_TOPLEFT, names, values))
            viewer.sync()
            time.sleep(1.0 / 60.0)


if __name__ == "__main__":
    main()
