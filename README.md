# SE(3) path following: a free rigid body and the Franka Panda

A C++17 implementation of the controller of the paper *Pose Path Following using Normal Tubular Coordinates*
(J. Dahlke and C. Nielsen): path following on SE(3) for fully actuated mechanical systems whose output is a pose,
including redundant arms. The library `se3pf` (Eigen 3.4) builds the normal frame and the tubular coordinates
(eta, xi) of a closed path, a lift q* and a redundancy output beta for a seven-joint arm, and the feedback-linearizing
controller. It reproduces the paper's two simulations, a free rigid body and the redundant Franka Panda, and checks
the Panda at the torque level in MuJoCo, the last step before the robot.

Documents:

- [`docs/DESIGN.md`](docs/DESIGN.md): the design document (interfaces, algorithms, tests and the expected numbers).
  It is a snapshot; the living version is the artifact
  <https://claude.ai/code/artifact/d6354d3d-b74d-42e8-bae0-b8a6e5c707e4>. Code comments cite its sections
  ("design doc, Online: projection") and the paper's equation labels (`eq:transport`, `Lemma 1`, ...).
- [`TICKETS.md`](TICKETS.md): the work, ticket by ticket: what is done, and what comes next.
- [`docs/MUJOCO.md`](docs/MUJOCO.md): the MuJoCo verification in detail.
- [`docs/HARDWARE.md`](docs/HARDWARE.md): the lab runbook for the experiments on the real Panda (libfranka 0.9.2).

## Layout

The library headers, in reading order, with the part of the paper each one implements:

| Header (`include/se3pf/`) | What it holds | Paper |
| --- | --- | --- |
| `types.hpp` | fixed-size aliases `Vec5`, `Vec6`, `Vec7`, `Mat65`, `Pose`, `Mat<R, C>` | Notation |
| `se3.hpp` | hat and vee, exp and log on SO(3), J_r, Ad, ad, poses | Notation, Example 1 |
| `metric.hpp` | `Metric`: W, inner product, Exp_I, Log_I, Gamma, distance (product type W = diag(wR I, wp I) only) | eq:metric, Example 1, eq:Gamma, eq:dist |
| `path.hpp` | `Path`; the coning circle (Panda) and the tumbling circle (rigid body), by arc length | eq:path-tangent |
| `spline.hpp` | `PeriodicHermite<T>`: periodic cubic Hermite interpolation | Implementation |
| `tube.hpp` | `Tube`: normal frame N, Psi, outputs (eta, xi), J_Psi; tube-radius estimate | Section IV: Lemma 1, eq:transport, eq:holonomy, eq:frame-ode, eq:Psi, Proposition 1 |
| `redundancy.hpp` | `RedundancyOutput`: beta, d beta and its level sets | eq:beta-flat, eq:dbeta |
| `controller.hpp` | `Controller<System>`, `Gains`, `TangentialReference`, the mu ramp, `stepCount` | Section V: eq:A, eq:transverse, eq:tan-vel, eq:tan-pos, eq:beta-law, eq:fb, Algorithm 1 |
| `rigid_body.hpp` | `RigidBody`: Euler's equations, RK4 on SE(3); a `System` for the controller | Section III, eq:EL |
| `rigid_body_example.hpp` | the paper's rigid-body example and its sampled closed loop | Simulations |
| `panda.hpp` | `PandaArm`: modified-DH kinematics, body Jacobian, limits; `HardwareBounds` | Problem formulation (robot arms) |
| `lift.hpp` | `Lift` q* (rest posture), `KernelFrame` X, `LiftRedundancy`, `CouplingRedundancy`, start configurations | Assumption 1, Lemma 3, eq:beta-flat, eq:dbeta |
| `panda_example.hpp` | the paper's Panda example and its closed loop at the joint-acceleration level | Simulations |
| `safety.hpp` | `SafetyMonitor` (the hardware margins, tighter values only), `Violation`, `ApproachMotion` | Hardware (design doc) |
| `experiments.hpp` | the hardware experiments as data, `ExperimentSetup`, `ExperimentLoop` (one control cycle of a run), logs | Hardware |
| `preflight.hpp` | `DynamicsModel`, `preflight()`: the run simulated with every check before anything moves | Hardware |

The rest: `src/` (one `.cpp` per header that is not a template), `apps/` (the simulations), `apps/robot/` (the robot
driver, the only code that includes libfranka), `mock/franka/` (a simulated Panda behind the libfranka 0.9.2 API),
`sim/` (the MuJoCo layer `se3pf_mujoco`, never part of the library), `tests/` (Catch2 v3; `tests/sil/` runs the
driver against the mock), `tools/` (CI script, plots, MuJoCo fetch and video, libfranka build, the oracle program of
the Panda numbers), `cmake/`.

Conventions: twists are (omega, v), angular first; Jacobians are body Jacobians in the flange frame; indices are
zero-based in code (joint 1 is `q[0]`); eta lies in [0, L). Units: s, L, eta and xi are in the units of the metric W.
With W = diag(l^2 I, I) a rotation by theta counts as l theta metres: for the Panda (l = 0.1 m) a xi of 0.01 is 1 cm
or 0.1 rad; for the rigid body (W = I) radians and metres simply add.

## Build and test

Needs CMake 3.18 or later, a C++17 compiler (GCC 11 is used here) and Eigen 3.4. The first configure downloads
Catch2 v3.5.2 (or pass `-DFETCHCONTENT_SOURCE_DIR_CATCH2=<a local Catch2 source tree>`).

```bash
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release && cmake --build build-rel -j
./build-rel/tests                    # 153 test cases, about 7 s
./build-rel/tests_sil                # 10 test cases, about 45 s: the robot driver against the mock libfranka
cmake -S . -B build-dbg -DCMAKE_BUILD_TYPE=Debug && cmake --build build-dbg -j
./build-dbg/tests                    # 145 test cases, about 1 min: Eigen's checks and the no-allocation tests
bash tools/ci.sh                     # what every change must pass (below)
```

`tools/ci.sh` builds Debug and Release with warnings as errors (`-DSE3PF_WERROR=ON`), runs all tests in both, runs
`rigid_body_sim` and `panda_sim` (their exit codes check the runs) and, in Release, `tests_sil` (the robot driver
against the mock libfranka, every experiment end to end, every kind of violation through the safety stop, faults
and refusals). If MuJoCo was fetched, it builds the MuJoCo
stage and runs `mujoco_smoke`, `tests_mujoco`, `panda_mujoco` (every check of its summary must pass) and
`experiments_mujoco` (every hardware experiment at the torque level). With `LIBFRANKA_INCLUDE=<libfranka>/include` it
compiles the driver against the real libfranka headers, and with `FRANKA_DIR=<libfranka build dir>` it builds and
links `robot_driver` (and checks that the mock-only options are not in it). Its last line, `CI PASSED for source tree
<hash>`, names the validated tree (`tools/tree_hash.sh`). About 6 minutes. Without network access: `CATCH2_SRC=/path/to/Catch2 bash tools/ci.sh`.

Test tags:

- `[slow]`: the 20 s closed loops and the paper's tube radii. They run by default in Release and are hidden in
  Debug, where each closed loop has a short stand-in (`./build-dbg/tests "[slow]"` runs them anyway, for minutes).
- `[.timing]`: timings, meaningful in Release only: `./build-rel/tests "[.timing]"` (the Panda control step takes
  about 12 us; the paper reports 27 us).
- The random inputs are reproducible: the generator is reseeded before every test, so a test gives the same result
  in the full run, alone (`./build-rel/tests "Tube: J_Psi*"`) or in random order (`--order rand`).
- To see a test's numbers: `-s`, e.g. `./build-rel/tests "estimateTubeRadius: the paper*" -s`. `ctest` works too.

## The three examples

Each app writes CSV logs to `out/` (or the directory given as its first argument) and prints a summary with the
expected values in brackets. The plot scripts need numpy and matplotlib.

**Rigid body** (tumbling circle, W = I, 2 ms sample-and-hold, a pose and a velocity reference, 20 s each):

```bash
./build-rel/rigid_body_sim               # out/rigid_body_pose.csv, out/rigid_body_velocity.csv
python3 tools/plot_rigid_body.py         # out/rigid_body.png and .pdf: the paper's three panels
```

**Panda at the joint-acceleration level** (the paper's results: with an exact model, eq:fb makes qdd = qdd*):

```bash
./build-rel/panda_sim                    # out/panda_paper.csv, out/panda_self_motion.csv, out/panda_limits.csv
python3 tools/plot_panda.py              # out/panda.png and .pdf
```

**Panda in MuJoCo** (joint torques on the Menagerie model: masses, inertias, armature, damping):

```bash
bash tools/fetch_mujoco.sh               # MuJoCo 3.15.0 and the Menagerie Panda into third_party/ (once)
python3 -m venv .venv && .venv/bin/pip install -r tools/requirements-mujoco.txt   # once
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release -DSE3PF_WITH_MUJOCO=ON && cmake --build build-rel -j
./build-rel/tests_mujoco                 # about 2 s
./build-rel/panda_mujoco                 # runs A to D (about 15 s), every check in the summary; exit 0 if all pass
.venv/bin/python tools/plot_mujoco.py    # out/mujoco.png and .pdf
MUJOCO_GL=egl .venv/bin/python tools/render_mujoco.py --run C --stills 0,5   # out/mujoco_panda.mp4
.venv/bin/python tools/view_mujoco.py --run D --loop                         # live replay (needs a display)
```

## Results

The paper's numbers next to ours, from the apps' summaries and the tests, which check them (the radii with 1000
directions come from a separate run of `estimateTubeRadius`):

| Quantity | Paper | Ours |
| --- | --- | --- |
| Panda path length L (l = 0.1 m) | 0.770 | 0.77008765 |
| Tube radius, Panda path, l = 0.1 m | 0.125 | 0.1250 |
| Tube radius, Panda path, l = 1 m | 0.279 | 0.307 (the +-e_i and 6 random directions); 0.283 with 1000 directions |
| Tube radius, rigid-body path | 1.89 | 1.93; 1.87 with 1000 directions |
| Rigid body: start off the path | 25.5 deg, 0.44 m | 25.45 deg, 0.4389 m (\|xi0\| = 0.6245) |
| Rigid body: max \|\|xi\| - designed response\| | within 1e-3 | 9.69e-4 (pose reference), 9.51e-4 (velocity reference) |
| Panda: start off the path | 15.1 deg, 16.8 mm, \|xi0\| = 0.031 | 15.07 deg, 16.82 mm, 0.03122 |
| Panda: max \|\|xi\| - designed response\| (nu_d = 0.08, 20 s) | 2.4e-5 | 2.434e-5 |
| Panda: time of one control step | 27 us (2.1 GHz Xeon) | about 12 us (Release) |

Along the Panda run: max |beta| 1.52e-5, final |xi| 7.1e-7, final eta_dot 0.0800025, smallest joint margin 0.3592 rad
(joint 4, set by the path), smallest sigma_min(A) 0.0474, joint velocities at most 26 % of their limits, one Newton
update per cycle. The self-motion run (on the path, beta0 = 0.2, not in the paper) follows the designed beta to
1.51e-4 and ends back on the lift to 1.3e-5. The deviations are sample-and-hold effects: halving the control period
halves them. The l = 1 m tube radius is 9.9 % above the paper's with the default directions (PAPER-1 asks the paper
to state how many it samples).

## What MuJoCo showed

Details and every number in [`docs/MUJOCO.md`](docs/MUJOCO.md).

- **The exact model and the zero-order hold reproduce the design.** With tau = M qdd* + C qdot + G at every MuJoCo
  step, MuJoCo's acceleration equals qdd* to 4e-15 and the kinematics agree with MuJoCo's to 6e-15; the outputs follow
  their designed responses up to MuJoCo's Euler step, an O(h) effect that halves with h. With the torque computed at
  the start of each 1 ms period and held (what the robot does), the deviations are O(dt): 3.4e-5 for |xi| and 2.6e-4
  for |beta|; |beta| and the speed error halve with dt.
- **The hardware-like run keeps every margin.** Held torque, mu ramp over 1 s and joint damping that the controller
  does not know: joint margin 0.354 rad (bound 0.10), joint velocity 13 % (50 %), commanded acceleration 21 % (30 %),
  torque 25 % (60 %), commanded torque rate 1.3 N m/s (600), flange speed 0.08 m/s (0.5), |xi| 0.031 (0.04), |beta|
  0.015 rad (0.05). Without the ramp, the first command asks for 49 % of an acceleration limit and the torque rate
  jumps to 1164 N m/s.
- **Uncompensated joint damping leaves a steady error.** With 1 N m s/rad of damping the path speed eta_dot settles
  near 0.064 instead of nu_d = 0.08, and |xi| and |beta| stay around 1e-2: the controller has no integral action, so a
  constant drag leaves a constant error. The robot's friction will act the same way (CTRL-1).
- **An unmodeled payload makes the arm sag and stall.** With the Franka Hand (0.76 kg) on the arm but not in the
  model, the arm sags (|xi| up to 0.215) and stalls on the path (eta_dot goes to 0). With the hand's weight in G, as
  libfranka adds it when the end effector is configured, every margin holds and the run is almost the run without
  the hand. So configure the end-effector load in Desk before running on the robot.

## Hardware

The robot driver runs the paper's experiments on a Franka Emika Panda (FER) with libfranka 0.9.2, the last version
that supports the Panda (0.10 and later support only the FR3). The lab runbook, from the lab PC's setup to the logs
to bring back: [`docs/HARDWARE.md`](docs/HARDWARE.md).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSE3PF_WITH_LIBFRANKA=ON -DFranka_DIR=$HOME/libfranka/build
cmake --build build -j && ./build/tests && ./build/tests_sil
./build/robot_driver list                              # hold, loops, recover, paper, self-motion
./build/robot_driver <robot-ip> check                  # kinematics against franka::Model, the configured load
./build/robot_driver <robot-ip> preflight all          # every experiment with the robot's own model; no motion
./build/robot_driver <robot-ip> hold --dry-run --log=logs/hold_1.csv   # then without --dry-run; then loops, ...
python3 tools/plot_robot_log.py logs/hold_1.csv        # the log against its preflight prediction
```

Nothing moves unless the preflight, a simulation of exactly that motion with the robot's model, keeps every margin of
`HardwareBounds` at every millisecond (again from the measured start after the approach); every motion waits for
Enter alone (any other answer aborts); in the torque loop the safety monitor checks every margin at 1 kHz and on any
violation the safety stop, a damped hold at the configuration of the violation, brings the arm to rest. Run the
session from the tree `tools/ci.sh` validated (`tools/tree_hash.sh`). Without a robot: `robot_driver_mock` (always built
with the tests) is the same driver against `mock/franka`, and `tests_sil` runs it end to end, with faults;
`experiments_mujoco` runs every experiment at the torque level in MuJoCo. Never run on a robot yet (ticket HW-3).

## Next steps

- **The first session on the robot (HW-3)**: [`docs/HARDWARE.md`](docs/HARDWARE.md). HW-1 (the driver against a
  mock) and HW-2 (the preflight and the safety monitor) are done.
- **Robustness to friction and payload** (CTRL-1): integral action with clamped integrators, friction compensation,
  or the load in the model; the target is run C of MuJoCo with |xi| < 1e-3 and eta_dot within 2 % of nu_d.
- **Milestone 11 (optional, deferred)**: a general left-invariant W, with Exp_I and Gamma from the geodesic ODEs.
  Until then `Metric` accepts the product type only and refuses any other W.
- **The paper** (PAPER-1): state the Gram-Schmidt order of the normal frame, the details of the rest-posture lift
  and the directions of the tube-radius estimate.

Never tested on the robot yet: the driver is built and linked against libfranka 0.9.2 and runs end to end against
the mock, but has not talked to a robot.
