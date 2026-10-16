# Running the experiments on the Panda: the lab runbook

How to run the paper's experiments on a Franka Emika Panda (FER) with libfranka 0.9.2, from an empty lab PC to the logs
you bring back. Read it once completely before the session. Ticket HW-3 in [`TICKETS.md`](../TICKETS.md) is this
session.

**Status.** The code has been run in simulation only, never on a robot. What was checked, and how: the last section.

**Safety, in one paragraph.** Two people: one at the keyboard, one with the user stop (external activation device)
in hand for every motion, nobody inside the reach of the arm. Nothing moves unless the preflight, a simulation of
exactly that motion with the robot's own dynamic model, keeps every margin at every millisecond. Every motion waits
for Enter alone; any other answer aborts. During the torque loop a safety monitor checks every margin at 1 kHz; on the
first violation the arm is held where it is by a damped joint hold until it is at rest. The margins are
`HardwareBounds` (`include/se3pf/panda.hpp`); an experiment can only tighten them, never loosen them. Do not edit
them, and do not skip the preflight.

| Margin (checked at every control period) | Bound |
| --- | --- |
| distance of every joint to its position limits | >= 0.10 rad |
| joint velocity, each joint against its own limit | <= 50 % |
| commanded joint acceleration qdd* | <= 30 % |
| delivered joint torque, tau + G | <= 60 % (52 N m joints 1-4, 7.2 N m joints 5-7) |
| commanded torque rate | <= 600 N m/s (libfranka's rate limiter is on as well) |
| translational speed of the flange | <= 0.5 m/s (0.2 m/s in the moving experiments, 0.1 m/s in `hold` and `self-motion`) |
| transverse error \|xi\| | <= 0.04 (metric units, l = 0.1 m: 4 cm of position, or 0.4 rad = 23 deg of attitude) |
| redundancy output \|beta\| | <= 0.05 rad |
| sigma_min(A) | >= 1e-3 |

**The safety stop** (`SafetyStop`, `include/se3pf/experiments.hpp`): a damped hold at q_stop, the configuration where
the violation was seen, tau = -K_s (q - q_stop) - D qd (the robot adds G), with K_s = (300, 300, 300, 300, 50, 50, 20)
N m/rad and D = (30, 30, 30, 30, 3, 3, 1.5) N m s/rad, each joint clamped to 60 % of its torque limit, until \|qd\| <
0.01 rad/s has held for 0.1 s (at most 2 s); then the torque loop ends and the robot's own controller takes the arm, as
at the end of any libfranka motion. Every gain is at most 60 % of the same gain of libfranka's `joint_impedance_control`
example. Damping alone, the earlier stop, could not hold a weight the model does not know: with the hand missing in Desk,
the flange fell about 0.3 m during the 2 s stop in MuJoCo. The hold keeps it within 0.022 rad on every joint, and a
stop at 50 % of every velocity limit travels at most 0.063 rad, inside the 0.10 rad joint margin (section 9).

## 1. What to bring to the lab PC

Run the session from **exactly the tree that `bash tools/ci.sh` passed on**, not from another checkout: the last line
of `tools/ci.sh` is `CI PASSED for source tree <hash>`. Copy that tree without its build directories and downloads:

```bash
rsync -a --exclude 'build*' --exclude third_party --exclude .venv --exclude out --exclude logs ./ labpc:~/se3pf_hw/
# or: tar czf se3pf_hw.tgz --exclude='./build*' --exclude=./third_party --exclude=./.venv --exclude=./out --exclude=./logs .
```

On the lab PC, `bash tools/tree_hash.sh` must print the same hash. If the tree is merged into another checkout first
(the main repository), run the whole `tools/ci.sh` there (tests, tests_sil, MuJoCo, and with `FRANKA_DIR`) and use
the hash it prints; in particular the non-finite guards that keep a NaN out of the 1 kHz loop must be there
(`PeriodicHermite::eval` in `include/se3pf/spline.hpp`; `Tube::projectionCondition` and `Tube::outputs` in
`src/tube.cpp`).

The first CMake configure downloads Catch2 v3.5.2 (the tests). Without internet on the lab PC, copy a Catch2 v3.5.2
source tree too and configure with `-DFETCHCONTENT_SOURCE_DIR_CATCH2=<it>`; or configure with
`-DSE3PF_BUILD_TESTS=OFF` (then there is no `robot_driver_mock` and no `tests_sil` on that PC). Python 3 with numpy
and matplotlib plots the logs (`tools/plot_robot_log.py`); it can also run on your laptop afterwards.

## 2. The lab PC

1. **Realtime kernel.** libfranka needs a PREEMPT_RT kernel for the 1 kHz loop. Must pass:

   ```bash
   cat /sys/kernel/realtime    # prints 1 (an RT-patched kernel), or:
   uname -v                    # names PREEMPT_RT (or "PREEMPT RT")
   ```

   The driver checks this itself before any motion and refuses without it (exit 7): libfranka 0.9.2 reads only
   `/sys/kernel/realtime`, and reads an uninitialized value when the file is missing. `--no-rt` is accepted only for
   commands that do not move the robot (`check`, `preflight`, `--dry-run`). Franka's guide: "Setting up the
   real-time kernel" in the FCI documentation (<https://frankaemika.github.io/docs/installation_linux.html>).
2. **Realtime permissions** for your user (once, then log out and in). Must pass: `ulimit -r` prints 99.

   ```bash
   sudo addgroup realtime && sudo usermod -a -G realtime $(whoami)
   printf '@realtime soft rtprio 99\n@realtime soft priority 99\n@realtime soft memlock 102400\n@realtime hard rtprio 99\n@realtime hard priority 99\n@realtime hard memlock 102400\n' | sudo tee -a /etc/security/limits.conf
   ulimit -r     # 99 after logging in again
   ```

   Franka also recommends the `performance` CPU frequency governor (`cpupower frequency-set -g performance`, or in
   the BIOS) and no power saving on the network card.
3. **Packages:** `sudo apt install build-essential cmake git libpoco-dev libeigen3-dev`. The project needs Eigen 3.4
   (Ubuntu 22.04 has it; on 20.04 install Eigen 3.4 from source: it is header only).
4. **libfranka 0.9.2**, the last version for the Panda (0.10 and later support only the FR3; CMake refuses them).
   Either with the packages above (then `Franka_DIR=$HOME/libfranka/build`):

   ```bash
   git clone --recursive --branch 0.9.2 https://github.com/frankaemika/libfranka ~/libfranka
   cmake -S ~/libfranka -B ~/libfranka/build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=OFF
   cmake --build ~/libfranka/build -j
   ```

   or without sudo (no libpoco-dev): `bash tools/build_libfranka.sh ~/opt/franka` builds Poco and libfranka from
   source (this is how the libfranka the driver was linked against was built); it prints the `Franka_DIR` to use
   (`~/opt/franka/libfranka/build`). Remove any ROS `libfranka` package that could be found first
   (`dpkg -l | grep -i franka`).
5. **Robot system version.** libfranka 0.9.2 needs Panda system 4.2.1 or newer (robot server version 5); see Franka's
   compatibility table <https://frankaemika.github.io/docs/compatibility.html>. Desk shows the system version;
   `robot_driver <ip> check` prints the server version. If libfranka refuses to connect with "Incompatible library
   version (server version: X, library version: 5)", the driver prints what X needs: **5** libfranka 0.9.2;
   **4** (system 4.0 to 4.2.0) libfranka 0.8.0, configured with `-DSE3PF_ALLOW_LIBFRANKA_08=ON` (the driver compiles
   against the 0.8.0 headers but was never linked or run with them: say so in the logs); **3** (system 3.x)
   libfranka 0.7.1, which this driver does not support: stop, or update the robot's system.
6. **Network.** Connect the PC to the LAN port of the **Control** (not of the arm), with a static address in the
   robot's subnet, e.g. robot 172.16.0.2, PC 172.16.0.1/24 (your lab's addresses may differ). `ping <robot-ip>`
   must answer; Desk opens at `https://<robot-ip>` in a browser. The robot sends its state by UDP: the firewall must
   let it in (`sudo ufw status` inactive, or a rule allowing the robot's address).

## 3. Desk (in the browser)

1. Unlock the joints (open the brakes) and release the user stop: the robot is then ready (libfranka's robot mode
   `Idle`, which `robot_driver <ip> check` prints; the driver moves only from `Idle`).
2. **End effector** (Settings, End effector): the mass, centre of mass and inertia of what is mounted on the flange
   (the Franka Hand: 0.73 kg, centre of mass (-0.01, 0, 0.03) m, as Desk's default for the hand; no end effector:
   0 kg). libfranka adds its weight to the gravity compensation, and `franka::Model` includes it in M and C. An end
   effector missing here makes the arm sag and stall (MuJoCo run D: \|xi\| 0.215 with the 0.76 kg hand not modeled);
   the driver prints the configured load and asks you to confirm it before every run.
3. **Activate FCI** (the menu at the top right). Then nothing else may control the robot (close Desk's apps).
4. The user stop: pressing it stops the arm at once; libfranka then ends the program with an error (exit code 8). To
   continue: release it, then `robot_driver <ip> error-recovery` (or acknowledge in Desk).

## 4. libfranka's own checks (before our code)

From libfranka's `build/examples` (no project code involved):

```bash
./echo_robot_state <robot-ip>      # no motion: prints the robot state (check robot_mode Idle, q, m_ee)
./communication_test <robot-ip>    # MOVES the arm to the ready pose (user stop in hand!), then 10 s of zero torque
```

`communication_test` checks the 1 kHz link: it must not print "WARNING: THIS SETUP IS PROBABLY NOT SUFFICIENT FOR
FCI!" (the success rate should be above 0.95, ideally 1.0). During its 10 s of zero torque (gravity compensated by
the robot) the arm must stay where it is. **If it drifts or sags, press the user stop at once**: the end effector
configured in Desk is wrong. It leaves the arm in the ready pose, which is near every experiment's start (our
approach moves at most 1 rad per joint).

## 5. Build and test on the lab PC

```bash
cd ~/se3pf_hw
bash tools/tree_hash.sh  # the hash of 'CI PASSED for source tree <hash>' (section 1)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSE3PF_WITH_LIBFRANKA=ON -DFranka_DIR=<Franka_DIR of section 2.4>
cmake --build build -j
./build/tests            # about 10 s: the library, the safety monitor and stop, the preflight
./build/tests_sil        # about 45 s: the driver against the simulated robot, every experiment, faults, refusals
./build/robot_driver list
```

CMake prints `libfranka 0.9.2: <dir>`; it stops with an error for libfranka 0.10 or later, and for 0.8.x unless asked
(section 2.5). Both test programs must pass before the robot moves.

## 6. The session

Keep a record of the terminal: `script -f session.txt` (or pipe each command into `tee -a session.txt`). Use a log
directory: `mkdir -p logs`.

1. **Check** (no motion): connection, the realtime setup of this PC, robot mode, q, the configured load, and the
   kinematics: PandaArm (our modified DH model) against `franka::Model` at the flange (pose and body Jacobian,
   libfranka's (v, omega) rows swapped to (omega, v)), at the robot's q, the rest posture, every start and eight points
   of the lift. Must print `OK` (agreement to 1e-5; we expect about 1e-7). It prints the position and rotation errors
   and the omega and v rows of J_b separately, and J_b without the row swap. A `MISMATCH` means stop: do not run
   anything; bring the output back (an error of order 1, or J_b fine only without the swap, is a convention error;
   1e-5 to 1e-3 is numeric; never loosen the tolerance).

   ```bash
   ./build/robot_driver <robot-ip> check
   ```

2. **Preflight of everything** (no motion), with the robot's own model and load. Every experiment must print `OK`.
   The predictions go to `logs/preflight_<experiment>.csv`.

   ```bash
   ./build/robot_driver <robot-ip> preflight all --log-dir=logs
   ```

3. **Each experiment**, in this order, first as a dry run (connects, checks, preflight; no motion), then for real
   with the same `--log` (a run refuses to overwrite an existing log; `--overwrite` allows it):

   ```bash
   ./build/robot_driver <robot-ip> hold --dry-run --log=logs/hold_1.csv
   ./build/robot_driver <robot-ip> hold --log=logs/hold_1.csv
   python3 tools/plot_robot_log.py logs/hold_1.csv          # look before the next experiment
   ```

   A run first checks that its log and `<log>.preflight.csv` can be written (exit 9 otherwise, before anything
   moves), then asks: (1) does the printed load match the mounted end effector (`y`); (2) Enter to move to the start
   configuration (a slow joint motion, at most 1 rad per joint, at least 3 s; it prints how far and how long);
   (3) Enter to start path following (torque control). **Only Enter alone continues**: `n`, `no`, `q`, any other text
   and end of input abort, and whatever was typed before a prompt (an extra Enter during the approach) is discarded.
   After each Enter the driver reads the robot's state again and checks it (Idle, at rest, the arm not moved since
   the prompt). Between (2) and (3) it checks that the arm reached the start (joints within 0.01 rad, \|xi - xi0\| and
   \|beta - beta0\| within 2e-3), prints the measured start next to the bounds (`At the start: |xi| = ... [bound
   0.040], |beta| = ... [bound 0.050]`), and runs the preflight again from the measured start (it must print `OK`;
   that prediction replaces `<log>.preflight.csv`). If the arm is farther than 1 rad (any joint) from the start, it
   refuses: guide it by hand (Desk's guiding mode, or the ready pose of `communication_test`). At the end it prints a
   summary and writes the log; the exit code says how it ended (section 8).

### The experiments, and what to expect

The predictions: the preflight (the closed loop at the joint-acceleration level, which the paper's numbers are) and
MuJoCo (`experiments_mujoco`: the Menagerie Panda at the torque level, the command held for 1 ms, libfranka's filter,
and 1 N m s/rad of joint damping the controller does not know, standing for the robot's friction). The robot's
friction will act like MuJoCo's damping: the path speed lags nu_d and \|xi\|, \|beta\| keep an error of about 1e-2
while moving (the controller has no integral action; ticket CTRL-1).

| # | Experiment | Duration | What moves | Preflight (kinematic) | MuJoCo, torque level with damping | Pass |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | `hold` | 10 s | nothing (holds eta = s0 = 0.1) | everything 0 | errors 1e-14; torque 21.6 % (gravity) | no stop; \|xi\|, \|beta\| < 2e-2 |
| 2 | `loops` | 38 s | two loops at nu_d = 0.05 (ramps of 4 s) | \|qd\| 7.3 %, \|qdd*\| 1.3 %, flange 0.049 m/s | eta_dot 0.035 while cruising; \|xi\| 8e-3, \|beta\| 0.010 rad; \|qdd*\| 8 % | no stop; \|xi\| < 2e-2; \|beta\| < 3e-2 |
| 3 | `recover` | 20 s | from Psi(0.1, xi0/2) (7.5 deg, 8.4 mm off), nu_d 0 -> 0.08 -> 0 | \|xi\| 0.0156 at the start, \|qdd*\| 4.6 % | \|beta\| 0.015, steady \|xi\| 1.1e-2, eta_dot 0.059; \|qdd*\| 11.5 % | no stop; \|xi\| decays in ~2 s, < 2e-2 after |
| 4 | `paper` | 20 s | the paper's start Psi(0.1, xi0) (15.1 deg, 16.8 mm off) | \|xi\| 0.0312, \|qd\| 18.8 %, \|qdd*\| 9.6 % | \|qdd*\| 17 %, torque 25 % (29.5 % with the hand), \|beta\| 0.015, steady \|xi\| 1.1e-2 | no stop; \|xi\| follows the design ~2 s, < 2e-2 after |
| 5 | `self-motion` | 6 s | the elbow, from beta = 0.04 rad back to the lift; the flange stays | \|qd\| 1.6 % | final beta 1e-4, \|xi\| 1e-4 | no stop; flange still; \|beta\| < 1e-2 after 3 s |

Every run starts at rest, ends at rest (the reference ramps down; the run ends once \|qd\| < 0.02 rad/s, at most 5 s
after the reference), and starts its torque from exactly zero (mu is ramped in over 1 s with the quintic smoothstep).
Joint margin along the path: 0.359 rad at least (joint 4, set by the path; the preflight's number, 0.354 in MuJoCo's
closed loop). `robot_driver list` prints each experiment's purpose and pass criteria.

The headroom at the start is thin in two experiments: `paper` starts at \|xi\| = 0.0312 against the 0.04 bound and
`self-motion` at \|beta\| = 0.04 against 0.05. Compare the measured start the driver prints with these numbers, and
`hold`'s steady \|xi\| and \|beta\| (the gravity model's error at rest) with the remaining headroom before running them.

- **hold** tests the computed-torque loop with the robot's model. When the torque loop starts, the robot's joint
  controller lets go and the feedback ramps in over 1 s (the first command is exactly zero), so an error of the
  gravity model or the load can let the arm drift by millimetres in the first fraction of a second before the
  feedback pulls it back. A steady \|xi\| or \|beta\| above 1e-2 at rest points at the gravity model or the
  load: check the end effector in Desk before moving on.
- **loops** measures the friction while moving. If \|beta\| exceeds 3e-2 here, `paper` will probably reach its
  0.05 rad margin and stop (safely); run `recover` first, and report the numbers rather than loosening anything.
- **recover** then **paper**: the paper's experiment, first with half the start offset. Overlay the prediction:
  `python3 tools/plot_robot_log.py logs/paper_1.csv --overlay out/mujoco_paper_hand.csv:MuJoCo` (the MuJoCo file
  comes from `experiments_mujoco`, run on a machine with MuJoCo, e.g. your laptop afterwards; the preflight's
  prediction, `logs/paper_1.csv.preflight.csv`, is overlaid automatically).
- **self-motion** shows the redundancy channel: only the elbow moves.

If you have time left: repeat `paper` two more times (repeatability), and `hold` once more at the end.

## 7. When something stops

- **Safety stop** (exit 5, `SAFETY STOP: <reason> at t = ...; held at the stop configuration until t = ...: largest
  joint drift ... rad, final |dq| ... (at rest)`): our monitor saw a margin violated and the damped hold brought the
  arm to rest where it was (MuJoCo: within about 0.6 s, a few hundredths of a radian); then the robot's own
  controller took it and stays ready (Idle). Plot the log (the stop is shaded) and find why before the next run: a
  push or contact (Acceleration, FlangeSpeed), the load (Tube or FlangeSpeed soon after the start, the arm sagging),
  friction (Redundancy or Tube while moving), timing (NonFinite or Projection after late periods). Running the same
  experiment again is safe; changing a margin is not an option.
- **The arm was not at rest when the stop ended** (the summary's `WARNING`, after 2 s of the hold): something keeps
  pushing it (a contact, or a load far from what Desk has). Do not run again before you know what; check the end
  effector in Desk first.
- **Preflight failed** (exit 4): do not run that experiment. Nothing moved (or, from the measured start, the arm
  stays at the start). Bring back the printed peaks and the prediction; never loosen a margin to make it pass.
- **Robot reflex** (exit 8, `ROBOT REFLEX / CONTROL ERROR`): the robot stopped itself (a collision threshold, a
  limit, the user stop, `communication_constraints_violation` from late packets). The driver still writes the log,
  and libfranka's last states and commands go to `<log>.reflex.csv`. Recover with `robot_driver <ip> error-recovery`
  (or in Desk), then look at the log. Repeated `communication_constraints_violation`: the realtime setup (section 2).
- **Abort** (exit 6) and **not ready** (exit 7: no realtime kernel, robot mode not Idle, the arm moving, too far
  from the start, moved while waiting for Enter, or not at the start after the approach) never start the torque
  loop. After a 7 at the start, run the same command again: the approach starts from where the arm is.
- **Log** (exit 9): the log exists already or cannot be written (checked before anything moves), or, after a run,
  could not be written (the run's data is lost: check the disk before the next run).
- In doubt: press the user stop. It always wins.

## 8. Reference

**Commands** (`robot_driver` = `./build/robot_driver`):

| Command | Motion | What it does |
| --- | --- | --- |
| `robot_driver list` | no | the experiments with their purpose and pass criteria |
| `robot_driver <ip> check [--no-rt]` | no | realtime setup, state, configured load, kinematics check |
| `robot_driver <ip> preflight [experiment\|all] [--log-dir=DIR] [--no-rt]` | no | kinematics check, preflight with the robot's model |
| `robot_driver <ip> error-recovery` | no | `automaticErrorRecovery()` after a reflex (not the experiment `recover`) |
| `robot_driver <ip> <experiment> --dry-run [--log=FILE.csv] [--no-rt]` | no | everything up to the preflight |
| `robot_driver <ip> <experiment> [--log=FILE.csv] [--overwrite]` | yes | approach and the experiment (default log `robot_<experiment>_<date-time>.csv`) |

**Exit codes:** 0 done (or check, preflight, dry run passed); 1 usage (also `--no-rt` for a run that moves the
robot); 2 unknown experiment; 3 kinematics check failed; 4 preflight failed (before the approach, or from the
measured start after it); 5 safety stop; 6 aborted by the operator; 7 not ready or not at the start; 8 libfranka error
(connection, reflex, user stop, realtime); 9 log file (exists or not writable, checked before any motion; or not
written after the run).

**Logs** (bring all of them back, with `session.txt`):

- `<log>.csv`: one row per control callback: t, phase (0 following, 1 safety stop: the damped hold), eta, eta_dot,
  ref (nu_d or eta_d), e_tan, xi1..xi5, xi_norm, beta, beta_dot, sigma_min, ramp, newton, q1..q7, dq1..dq7,
  tau1..tau7 (the commanded torque, without G), tauG1..tauG7 (tau + G: what the joints deliver), qdd1..qdd7 (qdd*),
  compute_us (time in the callback), period_ms (libfranka's period: 1 normally, more after lost packets),
  success_rate.
- `<log>.preflight.csv`: the preflight's prediction of the same run with the robot's model, from the measured start,
  same columns.
- `<log>.reflex.csv`: after a reflex only, libfranka's record of its last states and commands.
- Plots: `python3 tools/plot_robot_log.py <log>.csv [--overlay <prediction>.csv[:label]]` writes `<log>.png` and
  `.pdf` and prints the peaks.

**The software** (the only code that includes libfranka is `apps/robot/robot_driver.cpp`):

- `include/se3pf/safety.hpp`: `Violation`, `SafetyMonitor` (the checks above; tighter values only), `ApproachMotion`
  (the minimum-jerk approach: 20 % velocity, 10 % acceleration, 5 % jerk, at least 3 s).
- `include/se3pf/experiments.hpp`: the experiments as data (`experiments()`), `ExperimentSetup` (start
  configuration, controller with the gains, reference and mu ramp), `ExperimentLoop` (one control cycle: the
  controller step and the checks, shared by the driver, the preflight and the MuJoCo closed loop, so all three run
  the same logic), `SafetyStop` (the damped hold), the log rows and their peaks.
- `include/se3pf/preflight.hpp`: `DynamicsModel` and `preflight()` (from the setup's start or a given one).
- `apps/robot/robot_driver.cpp`: the driver; `franka_convert.hpp` (array conversions), `robot_log.hpp` (CSV).
- `mock/franka/`: a simulated Panda behind the libfranka 0.9.2 API, for `robot_driver_mock` and `tests_sil`.
- `sim/mujoco_experiments.hpp`, `apps/experiments_mujoco.cpp`: the experiments and the safety stop in MuJoCo.

**Changing an experiment:** `src/experiments.cpp` (plain data: start, reference, gains, ramp, tighter limits). Then
run `bash tools/ci.sh` (with MuJoCo) and `robot_driver <ip> preflight` again.

## 9. What was verified before the session, and what was not

Verified (all in `bash tools/ci.sh`, which passes):

- **The safety monitor**, every check for every joint against its own limit, at and beyond each bound, with NaN and
  infinity, and that no experiment can loosen a margin (`tests/test_safety.cpp`); the 1 kHz cycle, the checks and
  the safety stop do not allocate (`tests/test_controller.cpp`, counting operator new and Eigen's malloc check).
- **The safety stop** (`SafetyStop`): its law, clamp, handling of non-finite values and end conditions
  (`tests/test_preflight.cpp`). In MuJoCo (`tests_mujoco`, `experiments_mujoco`): with the hand missing in Desk it
  stops every experiment within 0.12 s and holds every joint within 0.022 rad, at rest 0.6 s later; from 50 % of
  every velocity limit it travels at most 0.063 rad (< 0.10 rad joint margin) and is at rest within 0.83 s (peak
  torque 40 N m); a 12 N m push on joint 4 deflects the hold by 0.04 rad, and the arm returns and rests.
- **The preflight** of every experiment passes: kinematic only, with a test model, with the robot's model in the mock,
  and with MuJoCo's arm and arm-plus-hand models; bad experiments fail with the right reason (a start outside the
  tube or the redundancy margin, huge gains, no mu ramp: acceleration and torque rate, a heavy model: torque, a
  reference that does not end at rest: timeout); from a given start (the measured one) as from the setup's
  (`tests/test_preflight.cpp`, `tests/mujoco/`).
- **The driver**, against the mock (`tests/sil/test_robot_driver.cpp`, 10 test cases): every experiment end to end
  (approach planned from q_d while the measured q differs, start check, preflight from the measured start, torque
  loop, finished at rest, log and prediction written, no allocation in the callback, the run within 1e-4 of its
  preflight in \|xi\|); every kind of violation goes through the safety stop and the arm is held at rest; a 12 N m
  push on joint 4 during `hold`; an end effector missing in Desk (0.73 kg: stopped within 0.2 s, held within
  0.028 rad; 1 kg: 0.035 rad); a NaN joint velocity; 1 N m s/rad of friction; lost packets (every 50th period
  3 ms): `paper` still finishes within its margins. It refuses to move on: a kinematics mismatch, a robot in reflex
  or user-stop mode (the behavior setters, which libfranka rejects outside Idle, come after the Idle check), an arm
  far from the start, a load that fails the preflight with the robot's model (exit 4), an approach that ends off the
  start (exit 7), a user stop pressed while the operator waits at the torque prompt (exit 7: the state is read
again after Enter), no realtime kernel (exit 8), the load not confirmed, any answer but Enter at either motion prompt
  (`n`, `no`, `q`, `abort`, a space) and end of input at any prompt, a log that exists or cannot be written (exit 9);
  a dry run moves nothing and works in any robot mode.
- **MuJoCo at the torque level**: every experiment ends at rest within every margin, with the arm alone and with the
  hand as a load configured in Desk (`tests/mujoco/test_mujoco_experiments.cpp`, `experiments_mujoco`).
- **libfranka 0.9.2**: `robot_driver.cpp` compiles with `-Wall -Wextra -Werror` against the real 0.9.2 headers (tag
  0.9.2, commit f3b8d775), and `robot_driver` is built and linked against a full libfranka 0.9.2 build
  (`tools/build_libfranka.sh`); `robot_driver list` runs, the test-only options `--yes` and `--inject` are not in
  it, `--no-rt` is refused for a run, and connecting to an address without a robot fails cleanly (exit 8). It also
  compiles against the 0.8.0 headers.

Not verified (only the robot can):

- Nothing has run on a real robot, and `robot_driver` has never talked to a real robot or to libfranka's model
  library (`franka::Model` is downloaded from the robot at run time). The mock reproduces libfranka's API, filter,
  rate limiter, reflex types, mode rules (motions and setters only when Idle), q_d and state fields as documented,
  not the robot; its kinematics are PandaArm's, so its kinematics check passes by construction.
- The kinematics check against Franka's own model, and the torque predictions with it (the preflight on the robot
  will show them).
- The safety stop on the robot: each gain is at most 60 % of libfranka's joint impedance example's, and it
  was tuned on MuJoCo's joint inertias, which include a 0.1 kg m^2 armature per joint.
- What the robot does when the torque loop ends: as at the end of any libfranka motion, its own controller takes
  the arm (in the mock and MuJoCo the arm is at rest by then).
- The robot's friction and the steady errors it causes (expected about 1e-2, as MuJoCo's damping), the joint
  elasticity, and the noise of the measured joint velocities in the derivative terms.
- Timing on the lab PC: the control step takes 12-40 us here; libfranka needs the answer within 1 ms.
