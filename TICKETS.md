# Tickets

Work one ticket at a time, in order. A ticket is done when it meets its acceptance criteria **and** the definition
of done below. Then ask for review ("TUBE-1 ready for review").

Status on 2026-10-15: milestones 1 to 10 of the design doc (`docs/DESIGN.md`) and the MuJoCo verification are done
and closed below. Milestone 11 (a general W from the geodesic ODEs) is deferred: Metric accepts the product type
only and refuses any other W. Milestone 12 (the robot): HW-1 and HW-2 are done (in simulation and against a mock and
the real libfranka 0.9.2, revised after the safety review; not yet on a robot), HW-3 is the first session on the
robot. Open: HW-3, CTRL-1, PAPER-1.

## Definition of done (every ticket)

- `bash tools/ci.sh` passes: Debug and Release build with zero warnings, all tests green.
- Every new public function has a test; every function that returns a derivative has a finite-difference test.
- Names and conventions follow the design doc (twists (omega, v), arc length s, zero-based indices, `se3pf`).
- Comments cite the paper's equation labels where the code implements one (`eq:transport`, `Lemma 1`, ...).
- Online code (anything the 1 kHz loop will call) uses fixed-size types and does not throw.
- If the interface differs from the design doc, the doc says why (or is updated).

## Review

- Comments are tagged **blocking** (must change before merge), **nit** (small, should change) or **question**
  (answer it; change the code if the answer is not obvious from the code).
- Address every comment: change the code, or reply why not. Then ask for re-review.
- Approved means merged: the reviewer takes a new snapshot of `main`, and the ticket is closed.

---

## Closed tickets

### Milestone 5: Tube (Section IV)

#### TUBE-1: Normal frame N(s) (Lemma 1)  ·  status: closed

**Result.** Holonomy angles 0.161182 / 8.01e-5 (Panda, l = 0.1 m), 0.707417 / 3.09e-4 (l = 1 m), 2.166040 / 0.492067
(rigid body); node closure 3.7e-11 (Panda) and 4.6e-12 (rigid body); N^T W N = I to 4e-11.

**Context.** The outputs (eta, xi) need a smooth, L-periodic, W-orthonormal frame of the normal space along the path.
Design doc: Module tube, Offline: the normal frame.

**Scope.**
- `include/se3pf/tube.hpp`, `src/tube.cpp`: `class Tube`, constructed from `const Path&`, `const Metric&` and the
  grid size (default 2000). Store the frame in `PeriodicHermite<Mat65>`.
- Steps 1-5 of the design doc: Gram-Schmidt with the skip rule, RK4 transport (`eq:transport`), H, its projection
  onto a rotation and its logarithm, closed frame (`eq:holonomy`) with N' from `eq:frame-ode`.
- Public: `void frame(double s, Mat65& N, Mat65* dN = nullptr) const`, `double length() const`,
  `const Path& path() const`, `const Metric& metric() const`, plus diagnostics for the tests:
  the rotation angles of H, det H, and the node closure gap |N(s_n) - N(s_0)|.
- Guard against the temporary-argument trap: `Tube(const Path&&, const Metric&) = delete;`.

**Acceptance criteria** (Panda path with l = 0.1 m and l = 1 m, rigid-body path with W = I):
- N^T W N = I and N^T W zeta_gamma = 0 to 1e-9 at 50 random s.
- Node closure gap below 1e-9.
- det H = 1 and H^T H = I to 1e-10; Phi skew to 1e-12; exactly two rotation planes.
- Rotation angles of H to 1e-5: Panda l = 0.1: {0.161182, 8.0e-5}; Panda l = 1: {0.707417, 3.09e-4};
  rigid body: {2.166040, 0.492067}.
- N' from the spline against a central difference of N, to 1e-6.

**Out of scope.** Psi, the projection, J_Psi.

#### TUBE-2: Tubular map Psi and the outputs (Definition 1, Proposition 1)  ·  status: closed

**Result.** Round trip to 1.4e-10 (eta) and 4.0e-12 (xi); |xi| equals the sampled distance to 3.6e-9; Newton converges
from +-0.05 (at most 3 updates from +0.02 in the round trip); globalProject gives 0.1, Newton from L/2 the wrong roots
5.614 and 0.485 (a test).

**Context.** (eta, xi) = Psi^{-1}(g) by Newton's method on the projection condition. Design doc: Online: projection
and outputs.

**Scope.**
- `Pose Psi(double s, const Vec5& xi) const` (`eq:Psi`).
- `bool projectionCondition(const Pose& g, double s, double& phi, double& dphi, Vec6& e) const`.
- `bool outputs(const Pose& g, double guess, double& eta, Vec5& xi, int* iterations = nullptr) const`, with the
  Newton counting rule of the design doc (at least one update; count = updates).
- `bool globalProject(const Pose& g, double& eta, int samples = 720) const`.
- `double wrap(double s) const` into [0, L).

**Acceptance criteria.**
- Round trip: outputs(Psi(s, xi)) returns (s, xi) to 1e-9, for 20 random s and |xi| up to 0.06 (Panda, l = 0.1 m).
- |xi| equals the minimum of dist(gamma(s_k), g) over 20000 samples s_k, to 2e-6 (Proposition 1(c)).
- Newton converges from guesses offset by +-0.05; phi'(eta) = -1 to 1e-9 for a pose on the path.
- The rigid-body start Psi(0.1, (0.4, -0.3, 0.2, -0.1, 0.3)) round-trips (|xi| = 0.6245).
- globalProject returns 0.1 for both examples' start poses; Newton started at L/2 converges to the wrong roots
  5.614 (rigid body) and 0.485 (Panda). Write that as a test: it documents the pitfall.

**Out of scope.** J_Psi, the tube radius.

#### TUBE-3: Body Jacobian J_Psi  ·  status: closed

**Result.** J_Psi against central differences to 1.6e-9 (Panda, l = 0.1 m), 9.4e-11 (l = 1 m), 6.6e-10 (rigid
body); exact at xi = 0.

**Context.** The controller's A(q) needs J_Psi^{-1} (`eq:A`). Design doc: Online: J_Psi.

**Scope.** `Mat6 JPsi(double eta, const Vec5& xi) const`, columns ordered (eta, xi_1..xi_5).

**Acceptance criteria.**
- J_Psi d against (Psi(y)^{-1} (Psi(y + h d) - Psi(y - h d)) / 2h)^vee with h = 1e-6, to 1e-6, for 10 random y on
  each of the Panda (|xi| ~ 0.015) and rigid-body (|xi| ~ 0.3) tubes.
- At xi = 0: J_Psi = [zeta_gamma(eta) | N(eta)] exactly (to 1e-14).

#### TUBE-4: Tube-radius estimate  ·  status: closed

**Result.** 0.125041 (Panda, l = 0.1 m) [paper 0.125], 0.306615 (l = 1 m) [0.279; 0.283 with 1000 directions],
1.92649 (rigid body) [1.89], under 1 s each in Release. The test now runs in the Release stage of CI (SLOW_IN_DEBUG).

**Context.** Theorem 2's hypothesis needs |xi| below the tube radius; the paper estimates it numerically.
Design doc: Offline: tube-radius estimate.

**Scope.** `double estimateTubeRadius(const Tube&, const RadiusOptions& = {})` with the options of the design doc
(32 values of s, +-e_i plus 6 random directions with a fixed seed, 1000 path samples plus points near the foot,
doubling then 14 bisections).

**Acceptance criteria.**
- Within 10 % of 1.93 (rigid body) and 0.125 (Panda, l = 0.1 m).
- Both example start poses lie inside: |xi0| = 0.6245 < estimate and 0.0312 < estimate.
- Runs in under 5 s in Release. Tag the test `[.slow]` if it is slow in Debug.

### Milestones 6 to 10 and the MuJoCo verification

#### M6: Rigid body end to end  ·  status: closed

RigidBody, Controller<RigidBody>, rigid_body_sim, plot. Start 25.45 deg and 0.4389 m (|xi0| = 0.6245) [paper 25.5 deg,
0.44 m]; max ||xi| - design| 9.69e-4 (pose reference) and 9.51e-4 (velocity reference) [paper: within 1e-3], halved
with dt (4.75e-4 at 1 ms, 2.37e-4 at 0.5 ms).

#### M7: Tube radius  ·  status: closed

Done as TUBE-4.

#### M8: Panda kinematics  ·  status: closed

FK(0) to 1e-16, FK(q_r) to 4.3e-7 (the design doc's 6 digits), J_b against central differences to 1e-8 and against
MuJoCo to 5e-15; the limits table entry by entry.

#### M9: Lift, kernel frame, beta  ·  status: closed

Projected Newton continuation: closure 1.5e-13, q*(0) as in the table, smallest joint margin 0.3592 rad (joint 4,
s = 0.367), sigma_min(J_b) 0.2070, minimum check at least 4.10e-3; beta = 0 and dbeta X = +1 on the lift to 2e-13;
start configurations to 1e-10.

#### M10: Panda end to end  ·  status: closed

max ||xi| - design| 2.434e-5 [paper 2.4e-5], max |beta| 1.52e-5, final eta_dot 0.0800025; self-motion run: max
|beta - design| 1.51e-4, back on the lift to 1.3e-5; A against finite differences off the lift to 3e-10; no
allocation in step; 12 us per step in Release [paper 27 us].

#### MJ-1: MuJoCo verification at the torque level  ·  status: closed

Exact model: qacc = qdd* to 4e-15, FK and J_b against MuJoCo to 6e-15, the outputs on the design up to MuJoCo's O(h)
Euler step; zero-order hold: O(dt) (both halve); the hardware-like run (hold, ramp, damping) keeps every margin.
Uncompensated damping and an unmodeled hand leave steady errors (CTRL-1). Details: `docs/MUJOCO.md`.

#### CLEAN-1: Final review and cleanup  ·  status: closed

Metric refuses a W that is not of product type; CI builds with -Werror and runs the apps and mujoco_smoke;
`stepCount` checks the run length of all three simulators; PeriodicHermite::eval returns NaN for a non-finite s;
the test random generator is reseeded for every test (failures reproduce alone and in random order); shared test
fixtures (`tests/panda_fixture.hpp`, `tests/tube_fixture.hpp`, `tests/test_helpers.hpp`); `HardwareBounds` in
`se3pf/panda.hpp`; `README.md` and `docs/DESIGN.md`. Every app summary and log is unchanged to the last digit
(panda_limits.csv gained a tau_max column).

---

## Open tickets

### Milestone 12: Hardware (libfranka 0.9.2)

Runbook: [`docs/HARDWARE.md`](docs/HARDWARE.md).

#### HW-1: libfranka 0.9.2 driver against a mock  ·  status: done (2026-10-15; not yet run on a robot)

**Result.** `apps/robot/robot_driver.cpp` (the only code that includes libfranka; `franka_convert.hpp` for the
arrays, `robot_log.hpp` for the CSV log), built as `robot_driver` against libfranka 0.9.2
(`-DSE3PF_WITH_LIBFRANKA=ON -DFranka_DIR=...`; CMake refuses 0.10 and later, and 0.8.x unless
`-DSE3PF_ALLOW_LIBFRANKA_08=ON`) and as `robot_driver_mock` against `mock/franka/` (always, with the tests).
Commands `list`, `<ip> check`, `<ip> preflight [experiment|all]`, `<ip> error-recovery`,
`<ip> <experiment> [--dry-run] [--log=] [--overwrite]`; `--no-rt` only for commands that do not move the robot (the
real build refuses it for a run, and checks for a PREEMPT_RT kernel itself before any motion); exit codes documented.
The mock is a simulated Panda behind the 0.9.2 API (Robot, Model, RobotState, Duration, Torques, JointPositions,
MotionFinished, the exceptions with ControlException's log, RealtimeConfig, the low-pass filter and rate limiter as
libfranka applies them, reflexes, motions and setters only when Idle, q_d distinct from q), with faults from
environment variables (pushes, friction, a load missing in Desk or a heavy configured load, approach and hold errors,
late periods, no realtime kernel, a NaN velocity).

**Safety review (2026-10-15), fixed.** The stop is a damped hold (HW-2). Only Enter alone continues at a motion
prompt (any other answer and end of input abort; typed-ahead input is discarded). The state is read again after each
Enter, and the checks repeated. The behavior setters come after the Idle check. The log path is checked before
anything moves, and an existing log is never overwritten unasked. The preflight is run again from the measured
start. The kinematics check prints position, rotation, J_b omega and v rows, and J_b without the row swap.
`tools/tree_hash.sh` ties the lab PC's tree to the one `tools/ci.sh` validated.

**Verified.**
- Compiles with `-Wall -Wextra -Werror` against the real 0.9.2 headers (tag 0.9.2, f3b8d775) and against the 0.8.0
  headers; built and linked against a full libfranka 0.9.2 build (`tools/build_libfranka.sh`, Poco 1.11.8);
  `robot_driver list` runs, `--yes` and `--inject` do not exist in it, `--no-rt` is refused for a run, a connection
  to an address without a robot fails cleanly (exit 8, libfranka's message printed once).
- Kinematics check: PandaArm against `model.pose(Frame::kFlange)` and `model.bodyJacobian(Frame::kFlange)` at the
  robot's q, the rest posture, every start and 8 points of the lift, rows (v, omega) swapped, column-major arrays;
  fails above 1e-5 (the mock with a 0.1 mm flange error: exit 3, diagnosed as a position error). The conversions are
  tested against finite differences of the flange pose (`tests/test_franka_convert.cpp`).
- Torque loop: `ExperimentLoop::follow` (Controller::step, the monitor), tau = M qdd* + C qd from `franka::Model`'s
  RobotState overloads (the robot adds G), log preallocated and written after the loop (also after a reflex). The
  mock runs follow their preflight to 1e-4 in |xi| (`paper`: about 3e-5, the hold and filter effects, O(dt) as MuJoCo
  run B's 3.4e-5). The callback allocates nothing (operator new counted inside it, mock build: 0).
- `tests_sil` (10 test cases): every experiment end to end; every kind of violation through the stop path
  (`--inject`, mock only); a push, a missing load (0.73 and 1 kg), a NaN velocity, friction, late periods; refusals:
  kinematics, robot mode, distance, a configured load that fails the preflight (exit 4), an approach that ends off
  the start (exit 7), the approach planned from q_d with q 1 mm off, a user stop pressed while waiting at a prompt
(the state read again after Enter), no realtime kernel, every operator answer that
  is not Enter at both motion prompts and end of input at each, the log path (missing directory, existing log).
  Mutation check: of the review's 30 mutants five survived; the four realistic ones are now killed (the preflight
  result ignored, confirm continuing at end of input, the start check skipped, the approach planned from q; the
  fifth made the rate check stricter). Of 11 new mutants of the fixes, 9 are killed (confirm accepting text, setters
  before the Idle check, the log path or --overwrite unchecked, the stop without stiffness, without its settle time,
  without its clamp or its NaN latch, no re-check after Enter); two survive by design: ignoring the preflight from
  the measured start (it cannot fail within the 2e-3 start tolerance; `preflight()` from a start is unit tested) and
  ignoring a failed prediction write (the log path is checked before).

**Not verified.** Never connected to a robot; `franka::Model` (downloaded from the robot) never used; the robot's
behavior when the torque loop ends; timing on the lab PC unknown.

#### HW-2: Preflight and safety monitor  ·  status: done (2026-10-15; not yet run on a robot)

**Result.** `se3pf/safety.hpp` (`Violation`, `SafetyMonitor` from `HardwareBounds`, tighter per-experiment values
only, allocation-free noexcept checks; `ApproachMotion`), `se3pf/preflight.hpp` (`DynamicsModel`, `preflight()`:
exactly the robot's run, through its end condition, every check at every period, from the setup's start or a given
one), `se3pf/experiments.hpp` (the five experiments `hold`, `loops`, `recover`, `paper`, `self-motion`, all starting
and ending at rest; `ExperimentLoop`, the control cycle shared by the driver, the preflight and the MuJoCo closed
loop; `SafetyStop`). In the driver: preflight with the robot's model and load before every run and again from the
measured start, minimum-jerk approach (20 / 10 / 5 % of the limits, at least 3 s, at most 1 rad per joint), start
check, Enter before every motion, and on any violation the safety stop: a damped hold at the configuration of the
violation, tau = -K_s (q - q_stop) - D qd, K_s = (300, 300, 300, 300, 50, 50, 20) N m/rad, D = (30, 30, 30, 30, 3,
3, 1.5) N m s/rad, clamped to 60 % of each torque limit, until at rest for 0.1 s (at most 2 s).

**Changed by the safety review (needs the author's approval).** The brief asked for joint damping, tau = -D qd with
D = (10, 10, 10, 10, 3, 3, 1.5). Damping cannot hold a weight that G does not know: with the hand missing in Desk,
the flange fell 0.3 m during the 2 s stop in MuJoCo, and joint 2 drifted 0.5 rad in the mock. With the stiffness
added, D = 10 on joints 1-4 left the hold ringing (1.9 s to rest after a stop at 50 % speed). Every gain is at
most 60 % of the same gain of libfranka's joint impedance example (K 600/250/150/50, D 50/30/25/15).

**Verified.**
- Every check for every joint against its own limit, boundaries, NaN, and that a margin cannot be loosened
  (`tests/test_safety.cpp`, 9 test cases); no allocation in the cycle and the stop (`tests/test_controller.cpp`).
- The stop's law, clamp, non-finite handling and end conditions (`tests/test_preflight.cpp`).
- Every experiment passes the preflight kinematic only, with a test model, with the mock's model, with MuJoCo's arm
  and arm + hand models. The paper's motion: peak |qdd*| 9.6 % with the ramp (the nu_d ramp lowers the design doc's
  11.5 %), 40.4 % without it at t = 0 (49.4 % with nu_d = 0.08 from t = 0, as the design doc: refused), and without
  the ramp the torque-rate margin fails; bad experiments and bad starts fail with the right reason
  (`tests/test_preflight.cpp`).
- MuJoCo at the torque level (hold, libfranka's filter, 1 N m s/rad damping): every experiment ends at rest within
  every margin, also with the hand configured; `paper`: |qdd*| 17 %, torque 25 % (29.5 % with the hand), |beta|
  0.015 rad, steady |xi| 1.1e-2, eta_dot 0.059 instead of 0.08 (`experiments_mujoco`, `tests_mujoco`).
- The stop holds the arm: in MuJoCo, an unmodeled hand stops every experiment within 0.12 s and the hold keeps every
  joint within 0.022 rad, at rest 0.6 s later; from 50 % of every velocity limit it travels at most 0.063 rad (inside
  the 0.10 rad joint margin) and rests within 0.83 s; a 12 N m push on joint 4 is stopped within 20 ms. In the mock:
  0.028 rad for a 0.73 kg missing load, 0.035 rad for 1 kg.
- `self-motion` starts at beta = 0.04 rad (inside the 0.05 margin; the start tolerance after the approach is 2e-3,
  and the preflight is run again from the measured start).

**Not verified.** The torque predictions with the robot's own model (the preflight on the robot will print them);
the margins against the robot's friction (the `loops` experiment measures it); the stop's gains on the robot's own
joint inertias (tuned on MuJoCo's, which include a 0.1 kg m^2 armature).

#### HW-3: First session on the robot  ·  status: open (today)

**Context.** HW-1 and HW-2 are done; the runbook is [`docs/HARDWARE.md`](docs/HARDWARE.md). MuJoCo run D showed that
an unmodeled hand makes the arm sag, so the end effector must be configured in Desk (the driver prints the load and
asks for confirmation; the monitor stops a sagging arm and the stop holds it).

**Scope.**
- Lab PC: the tree `tools/ci.sh` validated (`tools/tree_hash.sh` prints its hash), PREEMPT_RT
  (`/sys/kernel/realtime` = 1 or `uname -v` names PREEMPT_RT), `ulimit -r` 99, the firewall open to the robot's UDP,
  libfranka 0.9.2, the project built with `-DSE3PF_WITH_LIBFRANKA=ON`, `tests` and `tests_sil` green. Desk: joints
  unlocked, end effector set, FCI active, user stop in hand.
- libfranka's `echo_robot_state` and `communication_test` (link quality, and the load: no drift under zero torque).
- `robot_driver <ip> check` (kinematics against franka::Model at 15 configurations), `preflight all`, then the five
  experiments, each as a dry run and for real, in this order: `hold`, `loops`, `recover`, `paper`, `self-motion`;
  plot every log against its prediction before the next run.

**Acceptance criteria.**
- FK and J_b agree with franka::Model to 1e-5 at every check configuration.
- Each of the five experiments passes its preflight with the robot's model (also from the measured start) and
  completes without a reflex or a safety stop, meeting its pass criterion (docs/HARDWARE.md, section 6): `hold` and
  `self-motion` without visible flange motion, `self-motion` with |beta| < 1e-2 after 3 s; `loops` with |xi| < 2e-2
  and |beta| < 3e-2; `recover` and `paper` with |xi| following the designed response for about 2 s, then < 2e-2.
- The measured start (|xi|, |beta| against their bounds), the steady |xi|, |beta| and eta_dot reported next to
  MuJoCo's (`experiments_mujoco`): friction is expected to act like its damping.
- The logs (`*.csv`, `.preflight.csv`, any `.reflex.csv`) and the terminal record brought back. A safety stop or a
  failed preflight is a valid, safe outcome to report, not a reason to change a margin.

**Out of scope.** Fixing the steady errors (CTRL-1); hardware figures for the paper.

### Controller

#### CTRL-1: Robustness to joint friction and payload  ·  status: open

**Context.** The controller has no integral action. In MuJoCo run C (joint damping 1 N m s/rad, not compensated)
eta_dot settles near 0.064 instead of nu_d = 0.08, and over the last 5 s |xi| reaches 1.05e-2 and |beta| 1.32e-2;
in run D (the 0.76 kg hand not modeled) the arm sags to |xi| = 0.215 and stalls (eta_dot -> 0), while D-G (the hand
in G) is almost run C. The robot's friction will act like run C's damping.

**Scope.** Options, to be compared in MuJoCo:
- integral action on xi, beta and the tangential error, with clamped integrators, as in the reference controller
  (`~/Desktop/panda_path_following/include/pf/controller.hpp`, off by default there);
- friction compensation (viscous and Coulomb terms; identified on the robot later);
- the end-effector load in the model (M, C qdot and G), and in Desk.
New gains default to off, so the paper's closed loops are unchanged.

**Acceptance criteria.**
- MuJoCo run C: over the last 5 s |xi| < 1e-3 and eta_dot within 2 % of nu_d (0.0784 to 0.0816), with every
  hardware margin kept (ramp included).
- With the new terms off, every existing test passes unchanged (the paper's numbers).
- The integrators are clamped (a test drives them into the clamp) and reset with Controller::reset.
- Run D-G still keeps every margin.

**Out of scope.** Identification of the robot's friction (HW-3 data); a general W.

### Paper

#### PAPER-1: State the Gram-Schmidt order and the rest-posture lift details in the paper  ·  status: open

**Context.** The normal frame is unique only up to a constant rotation, and xi0 is a coordinate vector in it, so the
start poses depend on the Gram-Schmidt order (reversed columns put the rigid body 15.50 deg and 0.563 m off instead
of 25.45 deg and 0.4389 m). The lift "nearest the rest posture" needs q_r, the solver and the continuation to be
reproduced. The tube-radius estimate depends on the sampled directions (Panda l = 1 m: 0.307 with the +-e_i and 6
random directions, 0.283 with 1000, paper 0.279).

**Scope.** The Implementation section of the paper (in the paper's repository, not here): Gram-Schmidt on
zeta_gamma(0), e_1, ..., e_6 in that order with the skip rule; q_r = (0, -pi/4, 0, -3pi/4, 0, pi/2, pi/4), the
Newton-type method (projected step, 1e-12), the continuation from q_r at s = 0 on 2000 intervals and the closure
check; the directions and samples of the tube-radius estimate.

**Acceptance criteria.** Every bracketed number of the paper (0.770; 0.125 / 0.279; 1.89; 25.5 deg, 0.44 m;
within 1e-3; 15.1 deg, 16.8 mm, 0.031; 2.4e-5; 27 us) can be reproduced from the stated procedure with this code,
or is updated to what the code gives.

**Out of scope.** Code changes.
