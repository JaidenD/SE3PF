# A mock of libfranka 0.9.2 (software in the loop)

The subset of the libfranka 0.9.2 API that `apps/robot/robot_driver.cpp` uses, with the same headers, names and
signatures, behind which a simulated Panda runs at 1 kHz: `robot_driver_mock` is the driver compiled against these
headers (CMake target, always built with the tests), `tests/sil/test_robot_driver.cpp` runs it end to end. The real
driver `robot_driver` is the same source compiled against the real libfranka (`-DSE3PF_WITH_LIBFRANKA=ON`).

Keep the two in step: every declaration here mirrors one of `libfranka/include/franka/*.h` at tag 0.9.2, and
`tools/ci.sh` compiles the driver against both (the real headers with `-fsyntax-only` when `LIBFRANKA_INCLUDE` is set,
and links it when `FRANKA_DIR` is set).

What the simulated robot does (`mock_arm.h`):

- Kinematics: `se3pf::PandaArm` (the flange is Franka's frame F), so the driver's kinematics check passes exactly.
- Dynamics: not the Panda's. A toy arm, diagonal joint inertias plus a rigid body at the flange (the end effector of
  `RobotState::m_ee`, a Franka Hand: 0.73 kg), so that M(q), C(q, qd) qd and G(q) are consistent (C qd from the
  Newton-Euler equation of the flange body). The torque predictions that matter come from `franka::Model` on the
  robot, and from MuJoCo (`docs/MUJOCO.md`).
- Torque control as libfranka does it: the command is low-pass filtered (against `tau_J_d`, 100 Hz by default) and
  rate limited (1000 N m/s), the "robot" adds G of its model, and the plant integrates at 0.5 ms (with viscous and
  Coulomb joint friction if asked for, `FRANKA_MOCK_FRICTION` and `FRANKA_MOCK_COULOMB`). Reflexes throw
  `franka::ControlException` like the robot: delivered torque above the limits, a joint past its position or velocity
  limit, an external torque above the collision threshold of `setCollisionBehavior`.
- Joint position motions: the arm follows the command; the first command must equal `q_d`; velocity, acceleration
  and jerk of the commanded signal beyond the Panda's limits, or a motion that finishes while moving, throw.
- `q_d` is the robot's last commanded (or held) position and the measured `q` differs from it by a small deflection
  (`FRANKA_MOCK_HOLD_ERROR`, default 1e-5 rad), so a motion planned from `q` instead of `q_d` fails as on the robot
  (start-pose reflex). In torque control `q_d` equals `q`.
- Robot modes as libfranka: a motion (`Move`) and the behavior setters (`setCollisionBehavior`,
  `setJointImpedance`, `setCartesianImpedance`) are rejected with `CommandException` unless the robot is `Idle`, with
  libfranka's message ("... command rejected: command not possible in the current mode ...").
- With `RealtimeConfig::kEnforce` and `FRANKA_MOCK_NO_RT=1`, `control` throws `RealtimeException` before the motion
  starts, as libfranka 0.9.2 does on a kernel without realtime support.

Still different from the robot: the kinematics are PandaArm's own (the driver's kinematics check passes by
construction; `FRANKA_MOCK_POSE_ERROR` makes it fail), the dynamics are a toy arm's, every period is exactly 1 ms
unless `FRANKA_MOCK_LATE` is set, and nothing of the network is simulated.

Faults and other settings for tests, from environment variables (read when the `franka::Robot` is constructed):

| Variable | Meaning |
| --- | --- |
| `FRANKA_MOCK_PUSH=j,tau,t0,T` | an external torque `tau` (N m) on joint `j` (1-7) from `t0` for `T` seconds after the torque loop starts |
| `FRANKA_MOCK_FRICTION=d` | viscous joint friction `d` (N m s/rad) in the plant, not in the model |
| `FRANKA_MOCK_COULOMB=c` or `c1,...,c7` | Coulomb friction (N m) in the plant, every joint or per joint: `-c sign(qd)`, regularized below 1e-4 rad/s and integrated implicitly, so a joint under a net torque below `c` stays put (the stall of 2026-10-19, CTRL-1) |
| `FRANKA_MOCK_PLANT_LOAD=m` | an extra mass `m` (kg) at the flange in the plant only (an end effector missing in Desk) |
| `FRANKA_MOCK_Q0=q1,...,q7` | the initial configuration (default: the rest posture) |
| `FRANKA_MOCK_MODE=idle\|reflex\|user-stopped` | the robot mode reported before any motion |
| `FRANKA_MOCK_POSE_ERROR=e` | an offset `e` (m) of the model's flange position (to test the kinematics check) |
| `FRANKA_MOCK_MAX_SECONDS=s` | the longest control loop (default 600 s) |
| `FRANKA_MOCK_CONFIGURED_LOAD=m[,cx,cy,cz]` | the end effector configured in Desk (and carried by the plant): mass `m` (kg), centre of mass (m, flange frame; default the hand's), the hand's inertia; default a Franka Hand (0.73 kg) |
| `FRANKA_MOCK_HOLD_ERROR=e` | measured `q` minus `q_d` while the robot holds the arm, every joint (rad, default 1e-5) |
| `FRANKA_MOCK_APPROACH_ERROR=e` | a joint position motion settles `e` rad short of its goal on every joint (`q_d` stays the command) |
| `FRANKA_MOCK_LATE=n,ms` | every `n`-th torque period lasts `ms` ms (a lost packet: the robot keeps the last command; the next callback gets the longer period; `control_command_success_rate` drops accordingly) |
| `FRANKA_MOCK_NO_RT=1` | no realtime kernel: `control` with `RealtimeConfig::kEnforce` throws `RealtimeException` |
| `FRANKA_MOCK_USER_STOP_AFTER_READS=n` | the user stop is pressed right after the `n`-th `readOnce` (the robot is `UserStopped` from then on): a state that changed while the operator waited at a prompt |
| `FRANKA_MOCK_NAN=j,t0` | the measured velocity of joint `j` (1-7) is NaN in one torque period, `t0` s after the torque loop starts |
