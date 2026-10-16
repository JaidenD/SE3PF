#!/usr/bin/env python3
"""Check the MuJoCo Menagerie Panda against the design doc's kinematics, and print its dynamics facts.

  .venv/bin/python tools/mujoco_model_check.py [--model third_party/.../panda_nohand.xml] [--png out/rest.png]

Kinematics: the modified-DH forward kinematics and body Jacobian of the design doc appendix, written here in numpy
independently of the C++ code, against MuJoCo at 200 random configurations within the joint limits:
  - the frames of bodies link1..link7 against the DH frames 1..7,
  - the flange (DH frame 7 moved 0.107 m along z) against body link7 moved the same way,
  - the body Jacobian, rows (omega, v) in the flange frame, against mj_jac at the flange point.
Dynamics: mj_fullM includes the joint armature, mj_rne does not; mj_inverse equals
M qdd + qfrc_bias - qfrc_passive - qfrc_constraint.
With --png, renders scene.xml (same directory as the model) at the rest posture q_r; set MUJOCO_GL=egl.
"""
import argparse
import os

import numpy as np
import mujoco

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DEFAULT_MODEL = os.path.join(ROOT, "third_party/mujoco_menagerie/franka_emika_panda/panda_nohand.xml")

# Design doc appendix: modified (Craig) DH parameters (a_{i-1}, d_i, alpha_{i-1}); joint i turns about z_i.
A = [0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088]
D = [0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0]
ALPHA = [0.0, -np.pi / 2, np.pi / 2, np.pi / 2, -np.pi / 2, np.pi / 2, np.pi / 2]
D_FLANGE = 0.107  # PandaArm::kFlangeOffset, typed again on purpose: this script is an independent oracle
Q_MIN = np.array([-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973])
Q_MAX = np.array([2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973])
Q_REST = np.array([0.0, -np.pi / 4, 0.0, -3 * np.pi / 4, 0.0, np.pi / 2, np.pi / 4])


def mdh(a, d, alpha, theta):
    ca, sa, ct, st = np.cos(alpha), np.sin(alpha), np.cos(theta), np.sin(theta)
    return np.array([[ct, -st, 0.0, a],
                     [st * ca, ct * ca, -sa, -d * sa],
                     [st * sa, ct * sa, ca, d * ca],
                     [0.0, 0.0, 0.0, 1.0]])


def forward_kinematics(q):
    """Flange pose and the seven joint frames, all in the base frame."""
    T = np.eye(4)
    frames = []
    for i in range(7):
        T = T @ mdh(A[i], D[i], ALPHA[i], q[i])
        frames.append(T.copy())
    return T @ mdh(0.0, D_FLANGE, 0.0, 0.0), frames


def hat(w):
    return np.array([[0.0, -w[2], w[1]], [w[2], 0.0, -w[0]], [-w[1], w[0], 0.0]])


def body_jacobian(q):
    """J_b = Ad_{g^-1} J_s, rows (omega, v); spatial column i is (z_i, p_i x z_i)."""
    g, frames = forward_kinematics(q)
    Js = np.zeros((6, 7))
    for i, T in enumerate(frames):
        z, p = T[:3, 2], T[:3, 3]
        Js[:3, i] = z
        Js[3:, i] = np.cross(p, z)
    R, p = g[:3, :3], g[:3, 3]
    ad_inv = np.zeros((6, 6))
    ad_inv[:3, :3] = R.T
    ad_inv[3:, 3:] = R.T
    ad_inv[3:, :3] = -R.T @ hat(p)
    return ad_inv @ Js, g


def rotation_angle(R):
    """Angle of a rotation matrix, accurate near zero (unlike arccos of the trace)."""
    w = 0.5 * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    return np.arctan2(np.linalg.norm(w), 0.5 * (np.trace(R) - 1.0))


def mujoco_flange(m, d):
    """Flange pose and body Jacobian (omega, v) from MuJoCo: body link7 moved 0.107 m along its z-axis."""
    link7 = m.body("link7").id
    R = d.xmat[link7].reshape(3, 3)
    p = d.xpos[link7] + D_FLANGE * R[:, 2]
    jacp = np.zeros((3, m.nv))
    jacr = np.zeros((3, m.nv))
    mujoco.mj_jac(m, d, jacp, jacr, p, link7)  # world-frame rows; needs mj_kinematics and mj_comPos
    Jb = np.vstack([R.T @ jacr[:, :7], R.T @ jacp[:, :7]])
    return R, p, Jb


def check_kinematics(m, d, samples=200, seed=1):
    rng = np.random.default_rng(seed)
    err = dict(link_pos=0.0, link_rot=0.0, flange_pos=0.0, flange_rot=0.0, jacobian=0.0)
    for _ in range(samples):
        q = Q_MIN + (Q_MAX - Q_MIN) * rng.random(7)
        d.qpos[:7] = q
        mujoco.mj_kinematics(m, d)
        mujoco.mj_comPos(m, d)
        Jb, g = body_jacobian(q)
        _, frames = forward_kinematics(q)
        for i in range(7):
            body = m.body("link%d" % (i + 1)).id
            err["link_pos"] = max(err["link_pos"], np.abs(d.xpos[body] - frames[i][:3, 3]).max())
            err["link_rot"] = max(err["link_rot"], rotation_angle(frames[i][:3, :3].T @ d.xmat[body].reshape(3, 3)))
        R, p, Jb_mj = mujoco_flange(m, d)
        err["flange_pos"] = max(err["flange_pos"], np.abs(p - g[:3, 3]).max())
        err["flange_rot"] = max(err["flange_rot"], rotation_angle(g[:3, :3].T @ R))
        err["jacobian"] = max(err["jacobian"], np.abs(Jb_mj - Jb).max())
    print("Kinematics, %d random q within the limits (max errors):" % samples)
    print("  link1..7 frames: position %.1e m, rotation %.1e rad" % (err["link_pos"], err["link_rot"]))
    print("  flange:          position %.1e m, rotation %.1e rad" % (err["flange_pos"], err["flange_rot"]))
    print("  body Jacobian:   %.1e" % err["jacobian"])
    link0 = m.body("link0").id
    print("  base link0: position %s, |R - I| = %.1e" % (d.xpos[link0], np.abs(d.xmat[link0].reshape(3, 3) - np.eye(3)).max()))
    for name, q in [("q = 0", np.zeros(7)), ("q_r", Q_REST)]:
        d.qpos[:7] = q
        mujoco.mj_kinematics(m, d)
        mujoco.mj_comPos(m, d)
        R, p, _ = mujoco_flange(m, d)
        print("  flange at %-5s p = %s, R rows %s" % (name, np.round(p, 6), np.round(R, 6).tolist()))
    for s in range(m.nsite):
        R = d.xmat[m.body("link7").id].reshape(3, 3)
        Rs = d.site_xmat[s].reshape(3, 3)
        Rrel = R.T @ Rs
        print("  site %s: rotated %.6f deg about the flange z-axis" % (m.site(s).name, np.degrees(np.arctan2(Rrel[1, 0], Rrel[0, 0]))))


def check_dynamics(m, d, seed=3):
    rng = np.random.default_rng(seed)
    d.qpos[:7] = Q_REST
    d.qvel[:7] = 0.5 * rng.normal(size=7)
    mujoco.mj_forward(m, d)
    M = np.zeros((m.nv, m.nv))
    mujoco.mj_fullM(m, d, M)
    bias = d.qfrc_bias.copy()
    # Mass matrix of the recursive Newton-Euler algorithm: columns from unit accelerations.
    M_rne = np.zeros((m.nv, m.nv))
    result = np.zeros(m.nv)
    for j in range(m.nv):
        d.qacc[:] = 0.0
        d.qacc[j] = 1.0
        mujoco.mj_rne(m, d, 1, result)
        M_rne[:, j] = result - bias
    print("Dynamics at q_r with a random qvel:")
    print("  diag(mj_fullM - mj_rne mass matrix) = %s (dof_armature = %s)" % (np.diag(M - M_rne), m.dof_armature))
    qdd = rng.normal(size=m.nv)
    d.qacc[:] = qdd
    mujoco.mj_inverse(m, d)
    # qfrc_constraint is zero for panda_nohand.xml (no contacts, limits inactive); panda.xml has the finger equality.
    tau = M @ qdd + bias - d.qfrc_passive - d.qfrc_constraint
    print("  |qfrc_inverse - (M qdd + qfrc_bias - qfrc_passive - qfrc_constraint)| = %.1e, |qfrc_constraint| = %.1e"
          % (np.abs(d.qfrc_inverse - tau).max(), np.abs(d.qfrc_constraint).max()))
    print("  qfrc_passive = -damping qvel: %.1e" % np.abs(d.qfrc_passive + m.dof_damping * d.qvel).max())


def render_rest(model_path, png):
    import imageio.v3 as iio
    scene = os.path.join(os.path.dirname(model_path), "scene.xml")
    m = mujoco.MjModel.from_xml_path(scene)
    d = mujoco.MjData(m)
    d.qpos[:7] = Q_REST
    mujoco.mj_forward(m, d)
    camera = mujoco.MjvCamera()
    camera.lookat[:] = [0.2, 0.0, 0.4]
    camera.distance, camera.azimuth, camera.elevation = 1.5, 135.0, -15.0
    os.makedirs(os.path.dirname(os.path.abspath(png)), exist_ok=True)
    with mujoco.Renderer(m, height=480, width=640) as renderer:
        renderer.update_scene(d, camera=camera)
        iio.imwrite(png, renderer.render())
    print("wrote %s (MUJOCO_GL=%s)" % (png, os.environ.get("MUJOCO_GL", "default")))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--png", help="render scene.xml at the rest posture to this file")
    args = parser.parse_args()
    m = mujoco.MjModel.from_xml_path(args.model)
    d = mujoco.MjData(m)
    print("MuJoCo %s, %s: nq = %d, nv = %d, nu = %d" % (mujoco.__version__, os.path.basename(args.model), m.nq, m.nv, m.nu))
    check_kinematics(m, d)
    check_dynamics(m, d)
    if args.png:
        render_rest(args.model, args.png)


if __name__ == "__main__":
    main()
