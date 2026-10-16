# The Panda in MuJoCo: the controller at the torque level

The paper's Panda results, and the design doc's Panda simulation, are at the joint-acceleration level: with an exact
model, the feedback eq:fb makes q̈ = q̈*, so the plant is a double integrator. The robot is not. This simulation
drives the MuJoCo Menagerie model of the Panda (masses, inertias, joint armature, joint damping) with joint torques,
as libfranka will, and checks the hardware margins of the design doc before anything moves. It is the last step
before the robot.

Nothing here is part of the `se3pf` library, which never includes MuJoCo. The MuJoCo layer is the separate library
`se3pf_mujoco` (`sim/mujoco_panda.hpp`, `.cpp`) with the app `panda_mujoco` and the test program `tests_mujoco`,
all built only with `-DSE3PF_WITH_MUJOCO=ON`.

## Fetch, build, run

```bash
bash tools/fetch_mujoco.sh                     # MuJoCo 3.15.0 and the Menagerie Panda into third_party/ (once)
python3 -m venv .venv && .venv/bin/pip install -r tools/requirements-mujoco.txt   # for the Python tools (once)

cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release -DSE3PF_WITH_MUJOCO=ON
cmake --build build-rel -j
./build-rel/tests_mujoco                       # about 2 s in Release (about 23 s in Debug, [slow] tests hidden there)
./build-rel/panda_mujoco                       # all runs (about 15 s): logs and panda_limits.csv in out/, summary;
                                               # exit code 0 if all pass
.venv/bin/python tools/plot_mujoco.py          # out/mujoco.png and out/mujoco.pdf
MUJOCO_GL=egl .venv/bin/python tools/render_mujoco.py --run C --stills 0,5   # out/mujoco_panda.mp4 (+ PNG stills)
.venv/bin/python tools/view_mujoco.py --run D --loop                         # live replay in MuJoCo's viewer
bash tools/ci.sh                               # with MuJoCo there: also builds build-ci-mujoco, runs mujoco_smoke,
                                               # tests_mujoco and panda_mujoco (all its checks must pass; logs in
                                               # build-ci-mujoco/out)
```

CMake finds MuJoCo in `third_party/mujoco-*` (or `-DMUJOCO_DIR=...`) and the models in
`third_party/mujoco_menagerie/franka_emika_panda` (or `-DMENAGERIE_PANDA_DIR=...`); `panda_mujoco [out dir] [model dir]`
takes both directories as arguments too. `mujoco_smoke` is the raw C-API tour of the same facts, without se3pf.

Each log `out/mujoco_<run>.csv` has one row per control period: t, η, η̇, ξ (5), |ξ| and its design, β and its
design, q, q̇, τ (the full joint torque, `ctrl`), G (`grav1`..`grav7`, the gravity part of τ: τ − G is what libfranka
would be sent), q̈*, σ_min(A), the joint margin, the flange speed, the torque rates of τ − G and of τ, ncon, nefc, the
ramp and the baseline's |ξ|, β and |q − q_baseline|.

The viewer (`view_mujoco.py`) needs a display: space pauses, backspace restarts, `--speed 0.25` slows it down. The
video and the viewer draw the path γ of the flange (white), the target frame γ(η(t)) (x red, y green, z blue), the
flange's trail over the last 2 s (in the run's color) and t, |ξ|, β. The hand is drawn translucent in every run so the
flange stays visible; only the plants of runs D and D-G carry it.

## The model, and how the code uses it

| Fact | Where it matters |
| --- | --- |
| `panda_nohand.xml`: nq = nv = nu = 7; joint i is `qpos[i]`, dof i and `ctrl[i]` (checked in the constructor) | every run except D |
| `panda.xml`: adds the Franka Hand (0.73 kg + two 15 g fingers), two finger joints, a finger actuator, an equality constraint between the fingers | run D (the plant), the video |
| Link frames `link1`..`link7` are the modified-DH frames; the flange is `link7` moved 0.107 m along its own z-axis (the `attachment_site` is turned 135° about z: its attitude is not used) | `MujocoPanda::flange`, run (0) |
| The actuators are position servos. The constructor turns the seven arm actuators into torque motors: gain 1, no bias, `ctrlrange` = `forcerange` = ±87 / ±12 N m, so `ctrl` = τ | every run |
| Joint armature 0.1 kg m² (part of M), joint damping 1 N m s/rad (`qfrc_passive` = −D q̇), no friction loss. Both integrators treat the damping explicitly: MuJoCo's Euler would integrate it implicitly, with (M + hD)⁻¹, so the constructor sets `mjDSBL_EULERDAMP` | M from `mj_fullM`; damping on only in C and D |
| `qfrc_bias` = C(q, q̇) q̇ + G(q); MuJoCo 3.15 has no `d->qM`: `mj_fullM` gives the dense M, `mj_mulM` the product M v | `MujocoPanda::dynamics`, `exactStep` |
| Timestep 1 ms (the model's default is 2 ms, implicitfast); automatic reset on a bad state switched off, so a blow-up shows as a warning | all runs report MuJoCo warnings (0) |
| No contacts along the lift (the run reports `ncon` and `nefc` at every step) | runs A–D: ncon = 0; nefc = 0 except in D: 1 to 3 rows, the finger coupling (always) and the fingers' open stops (the finger servo's target, 0.04 m, is exactly the joint limit) |

Two ways to apply the torque, both in `MujocoPanda`:

- **Exact step** (`exactStep`, runs A and A2): `mj_step1` evaluates M, C q̇ + G and the passive forces at the current
  state, τ = M q̈* + C q̇ + G − passive goes to `ctrl`, and `mj_step2` integrates. Then MuJoCo's `qacc` equals q̈*
  to round-off (4e-15), and so does the acceleration it integrates, also with the damping on: with MuJoCo's default
  implicit Euler damping, `qacc` would still read q̈* but the velocity update would miss q̇ + h q̈* by h × 0.030
  rad/s² at |q̇| ≈ 0.3 rad/s (3e-14 with the explicit damping the wrapper sets). `mj_step1`/`mj_step2` always
  integrate with semi-implicit Euler, q̇ ← q̇ + h q̈ and then q ← q + h q̇ (with the new q̇): q moves by h² q̈ per
  step instead of the exact h²/2 q̈ of a held acceleration. That is an O(h) difference to the design doc's update,
  verified by halving h (below). q̈* is held over the 1 ms control period; with `substeps` > 1 the torque is
  recomputed at every MuJoCo step. The exact step is exact for `panda_nohand.xml` only: in `panda.xml` the finger
  dofs couple to the arm through M, so `simulateMujoco` refuses an exact run with the hand.
- **Zero-order hold** (`holdStep`, runs B, C, D): τ = M q̈* + C q̇ + G from the state at the start of the period, held
  for 1 ms while RK4 integrates (ctrl stays constant over the four RK4 stages). This is what the robot does. The
  passive forces are not compensated: the damping plays the role of the robot's unmodeled friction.

## What each run checks, and the results

Measured with `panda_mujoco` (Release); expected values from the design doc's Testing section in brackets. Every
check in the app's summary passes.

**(0) Kinematics.** `PandaArm::output` against MuJoCo's flange pose and J_b = [Rᵀ jacr; Rᵀ jacp] (`mj_jac` at the
flange point on `link7`), at 1000 random q within the joint limits: pose 2.0e-15, J_b 5.2e-15; along run A (20001
configurations): 6.4e-15 [about 1e-15; fail above 1e-10].

**(A) Exact model**, `panda_nohand.xml`, damping off, the paper's start, 20 s. Next to it, the app runs the same
controller at the joint-acceleration level (q ← q + Δt q̇ + ½Δt² q̈*, q̇ ← q̇ + Δt q̈*), the *baseline*, which
reproduces the design doc to its printed precision.

| | MuJoCo, h = 1 ms | baseline | [design doc] |
| --- | --- | --- | --- |
| max \|qacc − q̈*\| | 4.3e-15 | | [round-off] |
| max \|\|ξ\| − design\| | 3.89e-5 (t = 0.29 s) | 2.434e-5 | [2.43e-5; paper 2.4e-5] |
| max \|β\| | 3.58e-5 | 1.52e-5 | [1.5e-5] |
| final \|ξ\| | 1.01e-5 | 7.09e-7 | [7.1e-7] |
| final η̇ | 0.0800049 | 0.0800025 | [0.0800025] |
| smallest joint margin | 0.3592 rad, joint 4 | | [0.3592] |
| smallest σ_min(A) | 0.0474 | | [0.0474] |
| largest \|q̇ᵢ\| / q̇_max | 26.2 % | | [26 %] |
| Newton updates per cycle | 1 | | [1] |
| contacts, constraint rows | 0, 0 | | |

The MuJoCo column differs from the baseline only by MuJoCo's Euler step: max |q_MuJoCo − q_baseline| is 3.90e-4,
1.95e-4 and 9.74e-5 rad for h = 1, 0.5 and 0.25 ms (ratios 2.00 and 2.00), and the deviation from the design falls
to 3.08e-5 and 2.73e-5 on its way to the baseline's 2.43e-5. The dynamics themselves are cancelled exactly.

**(A2) Self-motion**, exact model: start on the path at γ(0.1) with β₀ = 0.2 and A(q₀) q̇₀ = (0.08, 0, …, 0), one loop
(9627 steps).

| | MuJoCo, h = 1 ms | baseline | [design doc] |
| --- | --- | --- | --- |
| max \|ξ\| | 1.91e-5 | 1.24e-6 | [1.2e-6] |
| max \|β − β₀(1 + 4t)e^{−4t}\| | 2.37e-4 | 1.51e-4 | [1.51e-4, accept < 2e-4] |
| final maxᵢ \|qᵢ − q*ᵢ(η)\|, 2-norm | 5.6e-5, 6.5e-5 | 9.9e-6, 1.3e-5 | [9.9e-6, 1.3e-5 (< 2e-5)] |

Again the excess is Euler's: max |q_MuJoCo − q_baseline| is 1.28e-4 at h = 1 ms and 6.40e-5 at 0.5 ms (ratio 2.00),
where the three values above become 9.8e-6, 1.88e-4 and 2.8e-5. The design doc's acceptance bounds are for the
joint-acceleration level; the baseline meets them.

**(B) Zero-order hold**, damping off, 20 s: max ||ξ| − design| 3.44e-5, max |β| 2.61e-4, final |ξ| 3.28e-5, final
η̇ 0.0801011; joint margin 0.3591 rad, σ_min(A) 0.0474, |q̇| 26.2 % as in A. Within a period the torque was right
only at its start (M and C q̇ + G change with the state), so the outputs leave their designed responses by O(Δt):
with Δt = 0.5 ms, max |β| is 1.30e-4 (ratio 2.00) and η̇ − ν_d 5.0e-5 instead of 1.0e-4 (ratio 2.02). These are
small next to every margin; they would grow with a slower loop.

**(C) Hardware-like**: hold, μ ramp over T = 1 s (`Controller::setRampTime`), joint damping 1 N m s/rad on and not
compensated, 20 s. Every margin of the hardware section (design doc, Hardware; the bounds are `HardwareBounds` in
`include/se3pf/panda.hpp`, used by the app and the tests):

| margin | peak | bound |
| --- | --- | --- |
| distance to the joint limits | 0.3537 rad (joint 4) | ≥ 0.10 rad |
| \|q̇ᵢ\| / q̇_max | 13.4 % (joint 5) | ≤ 50 % |
| commanded \|q̈*ᵢ\| / q̈_max | 20.7 % (joint 5, t = 0.67 s) | ≤ 30 % |
| delivered \|τᵢ\| / τ_max (gravity included: MuJoCo's ctrl is the robot's τ_cmd + G) | 25.1 % (joint 2) | ≤ 60 % |
| torque rate of the commanded torque, \|d(τ − G)/dt\| | 1.3 N m/s | ≤ 600 N m/s |
| torque rate of the full torque, \|dτ/dt\| (for comparison) | 5.9 N m/s | |
| flange speed | 0.079 m/s | ≤ 0.5 m/s |
| \|ξ\| | 0.0312 (the start) | ≤ 0.04 |
| \|β\| | 0.0151 rad | ≤ 0.05 rad |
| σ_min(A) | 0.0474 | ≥ 1e-3 |
| contacts, clamped torques, MuJoCo warnings | 0, 0, 0 | 0 |

Every joint is measured against its own limit (`hardwareMargins`): the torque peaks per joint are 0.2, 25.1, 4.6,
22.0, 6.0, 15.7 and 1.9 % (87 N m for joints 1–4, 12 N m for 5–7), so joint 2 sets the torque margin on this path.
A wrong wrist limit would not show here; the unit tests of `hardwareMargins` let every joint in turn set every
margin, so it cannot hide there.

The damping is felt, and the app checks it: η̇ settles near 0.064 instead of ν_d = 0.08 (k_t = 5 cannot hold the
speed against a constant drag without integral action; check: < 0.07), and over the last 5 s |ξ| reaches 1.05e-2 and
|β| 1.32e-2 rad. Without the damping (the app's run "C without damping") η̇ ends at 0.0801 and these are 2.76e-5 and
1.47e-4 (checks: ratios above 10; measured 380 and 90): the damping, not the hold or the ramp, causes them. Were the
damping missing from the plant, or compensated in the torque, run C would be the run without damping and these
checks would fail. The baseline column of C confirms the design doc's preflight numbers for the ramp at the
joint-acceleration level: 11.5 % of the acceleration limit and 22.4 % of the velocity limit [11.5 %, 22 %].

**(C0) No ramp**: the first command asks for 49.4 % of a joint's acceleration limit [49 %] and the commanded torque
jumps from 0 to M q̈*(0) in one period, 1164 N m/s: two margins violated, as the design doc predicts. With the ramp
(C) the first command is exactly zero and the peaks are 20.7 % and 1.3 N m/s.

**(D) Model mismatch**: the plant is `panda.xml` (the hand, fingers held open by their actuator); M and C q̇ + G come
from a second `MujocoPanda` on `panda_nohand.xml`, set to the plant's arm state every period; hold, ramp and damping
as in C. The hand's weight is not compensated, and the feedback (poles at −4, no integral action) is soft: the
closed loop's stiffness at the flange is roughly m_eff λ², some tens of N/m, against 7.5 N of weight. The arm sags
while the ramp is still small (|ξ| up to 0.215 at t = 0.72 s), the commanded accelerations reach 110 % of the limit
and the joint velocities 60 %, and the arm then stalls on the path: η̇ → 0, because on that part of the path the
hand's weight pulls back harder than the tangential law can push (at most k_t ν_d = 0.4, in units of η̈, at η̇ = 0).
Steady state: |ξ| 0.0915, |β| 0.0118 rad. Margins violated: |ξ| (0.215 > 0.04), |q̇| (60 % > 50 %, joint 2), |q̈*|
(110 % > 30 %, joint 2); the rest hold (joint margin 0.121 rad, torque 32 %, flange speed 0.46 m/s, |β| 0.040 rad). No contacts.

**(D-G)** is D with the hand's weight in G, as on the robot when the hand is configured in Desk (libfranka then adds
G with the end-effector load). M and C q̇ still miss the hand. All margins hold, and the run is almost run C: last
5 s |ξ| 1.10e-2 and |β| 1.31e-2 (C: 1.05e-2, 1.32e-2), η̇ 0.0642 (C: 0.0643), torque 29.5 % (C: 25.1 %). The hand's
inertia hardly matters; its weight does.

The figure (`out/mujoco.png`): (a) |ξ| against its design for A–D, (b) the deviation on a log scale with the
baseline, (c) |β|, (d) run C's joint torques against their limits, (e) the torque rate of C and C0, (f) η̇ against
ν_d.

## The tests (`tests_mujoco`)

- the model in torque mode: options (Euler's damping explicit), the hand model, refusal of a model that is not a Panda;
- kinematics: FK and J_b against `PandaArm::output` to 1e-12 at 200 random q (also with the hand model);
- M symmetric positive definite with the armature, the bias at rest equals G, passive = −D q̇, and M q̈* + bias −
  passive equals the torque of `exactStep`;
- `exactStep` gives qacc = q̈* and the Euler velocity update q̇ ← q̇ + h q̈* to 1e-12, damping on or off;
- τ = `qfrc_bias` holds q_r for 2 s to 1e-10 (Euler and RK4), τ = 0 lets the arm fall;
- contacts of a folded arm, the limit row of a joint past its limit;
- `hardwareMargins` with synthetic inputs: every joint in turn sits at 50 % of its own velocity, acceleration and
  torque limit (the others at 20 %) and must set that margin (one limit for all joints puts a wrist joint at 7 %);
  the rates of τ − G and of τ (G of joint 2 rising at 10 N m/s, the commanded torque of joint 6 at 3); the first
  instant after rest; the flange speed of a twist with |v| = 0.5 m/s;
- the argument checks of `simulateMujoco`, among them an exact run with the hand (refused);
- 2 s exact-model loop from the paper's start: ||ξ| − design| < 5e-5 (measured 3.89e-5) and |β| < 5e-5 (3.58e-5),
  the baseline's 2.43e-5 and 1.52e-5, the first command 49.4 %, the log against the summary;
- 2 s zero-order hold: ||ξ| − design| < 5e-5 (3.2e-5), 1e-4 < |β| < 4e-4 (2.6e-4);
- 2 s hardware-like loop (run C): every margin, and each peak recomputed from the log (per joint, and both torque
  rates from the logged τ and G); the commanded torque rate below 2 N m/s (1.27); τ recomputed from the log with a
  model of its own, M q̈* + C q̇ + G and no passive term (to 1e-12); the damping felt: over each period the plant's
  mean acceleration minus q̈* against the damping's M⁻¹(−D q̇), least-squares factor 1 ± 0.01 (0.9995; it would be
  0 with the damping off or compensated), largest residual 3.1e-3 against 2.9 rad/s² of damping acceleration; in
  the outputs |β| > 5e-3 (1.39e-2; 2.6e-4 without damping) and η̇(2 s) < 0.07 (0.053; 0.080); the ramp's zero
  first command; without the ramp the acceleration and torque-rate margins fail;
- 1 s with the hand as payload (runs D and D-G): D's |ξ| reaches 0.215 (> 0.1), D-G keeps every margin (peaks
  recomputed from the log); τ of both recomputed from the log (arm-only M and C q̇ + G at the plant's state, in D-G
  with the plant's G) and the logged G;
- `[slow]` (hidden in Debug): halving h halves |q_MuJoCo − q_baseline|, and the self-motion run.

The bounds of the closed loops sit 25–50 % above the measured values. They are not tight to the 1e-15 level because
the Euler step and the hold add O(h) terms; a wrong torque (no M, M without its off-diagonal terms, no C q̇ + G, a
wrong sign of the passive force) or wrong kinematics fails them by far (each of these was tried).

Upper bounds alone are not enough: a run without the damping, or with the damping compensated, passes every margin
of run C more easily, and a wrong wrist torque limit does not show where joint 2 sets the peak. The review's mutation
tests found such holes (the damping off in the plant or compensated, a torque rate against a stale previous torque,
one torque limit for all joints, the separate model of D and D-G tested only by the app). They are closed by the
checks above that the damping is felt, by τ and G recomputed from the log, by the rates recomputed from the logged
G, by `hardwareMargins` with synthetic inputs and by the D/D-G test; a second round of 24 mutants of the wrapper,
the loop and the margins (those holes and their neighbours: half the damping compensated, the rate of τ for the
rate of τ − G, the flange's angular speed for its linear speed, implicit Euler damping, an exact run with the hand,
...) all fail `tests_mujoco`. CI also runs `panda_mujoco`, whose summary checks the full 20 s runs.

## The hardware experiments (`experiments_mujoco`)

The five experiments of the robot driver (`se3pf/experiments.hpp`, [`HARDWARE.md`](HARDWARE.md)) run on the MuJoCo
Panda exactly as the driver runs them: `ExperimentLoop` (the controller step and the safety monitor's checks), the
safety stop on a violation (`SafetyStop`, a damped hold), the same end conditions; the commanded torque M q̈* + C q̇ goes through libfranka's
low-pass filter (100 Hz) and rate limiter, G is added (as the robot does), and the sum is held for 1 ms on the plant
with its joint damping on and not compensated (as run C). Two plants: the arm (`panda_nohand.xml`), and the arm with
the Franka Hand as a load whose weight is in G but not in M and C q̇ (as run D-G: the hand configured in Desk). Each
experiment's preflight is also run with MuJoCo's dynamics (`MujocoDynamics`, a `DynamicsModel`) for both models.

```bash
./build-rel/experiments_mujoco          # about 25 s; out/mujoco_<experiment>_<arm|hand>.csv in the robot log format
python3 tools/plot_robot_log.py out/mujoco_paper_arm.csv --no-preflight
```

Every experiment ends at rest within every margin, with both plants (2026-10-15):

| experiment | \|q̈*\| | \|q̇\| | torque (arm / with hand) | \|ξ\| max | \|β\| max | steady \|ξ\| (last 5 s) | η̇ while cruising |
| --- | --- | --- | --- | --- | --- | --- | --- |
| hold | 0.0 % | 0.0 % | 21.6 / 25.7 % | 1e-14 | 1e-15 | 1e-14 | |
| loops | 8.1 % | 6.3 % | 25.1 / 29.4 % | 0.0081 | 0.0102 | 3.6e-3 | 0.035 (ν_d 0.05) |
| recover | 11.5 % | 9.2 % | 25.1 / 29.5 % | 0.0156 | 0.0152 | 1.1e-2 | 0.059 (ν_d 0.08) |
| paper | 17.4 % | 11.3 % | 25.1 / 29.5 % | 0.0312 | 0.0152 | 1.1e-2 | 0.059 (ν_d 0.08) |
| self-motion | 1.6 % | 0.9 % | 21.6 / 25.7 % | 1e-4 | 0.0400 (start) | | |

Commanded torque rate at most 0.9 N m/s, flange speed at most 0.079 m/s, joint margin at least 0.354 rad, no
contacts, no clamped torques, no MuJoCo warnings. The damping leaves the errors of run C (the path speed lags ν_d by
about 27 %), which on the robot friction will cause too. `tests_mujoco` checks all of this
(`tests/mujoco/test_mujoco_experiments.cpp`), and the safety stop: on a 12 N m push on joint 4 the monitor stops the
arm within 20 ms and the damped hold brings it back to rest (0.04 rad of deflection while pushed); with an unmodeled
hand the monitor stops every experiment within 0.12 s and the hold keeps every joint within 0.022 rad of where it
stopped, at rest 0.6 s later (damping alone, the earlier stop law, let the flange fall about 0.3 m in 2 s); from 50 %
of every velocity limit the hold travels at most 0.063 rad, inside the 0.10 rad joint margin, and is at rest within
0.83 s. `experiments_mujoco` prints these two stops for every experiment.

## What is still different on the real robot

- **libfranka adds G itself.** The torque sent is τ_cmd = M q̈* + C q̇ (from `franka::Model`); the robot adds its own
  G, which includes the end-effector load configured in Desk. In MuJoCo `ctrl` is the full joint torque, so the
  delivered-torque margin applies to `ctrl` and the torque-rate margin to `ctrl − G`. Run D shows why the hand (or any
  payload) must be configured: unmodeled weight is the largest error this controller can meet.
- **The model is not the robot's.** The Menagerie inertias come from `franka_description`; `franka::Model` is
  Franka's own. The preflight must use `franka::Model` for its torques, not these.
- **Friction.** Real joints have Coulomb and viscous friction, and the Panda's joint torque control (with its torque
  sensors) removes part of it. What remains acts like run C's damping: the path speed lags ν_d and |ξ| and |β| settle
  at a steady offset (about 1e-2 in run C, with the model's 1 N m s/rad). Nothing in the controller removes a steady
  error; friction compensation or integral action on the outputs would be the next step if the robot shows it
  (ticket CTRL-1 in TICKETS.md).
- **Joint elasticity, torque sensor noise, measured q̇.** The robot's q̇ comes from encoders (and is filtered), its
  torque from sensors; MuJoCo gives the exact state. The derivative terms (K_d ξ̇, the central difference for Ȧq̇)
  will see that noise.
- **Timing.** libfranka calls the control callback every 1 ms and expects the answer within it (in Release, run A
  with its baseline, two controller steps and one MuJoCo step per period, takes 0.8 s for 20 000 periods); late
  packets and jitter do not exist here. libfranka can also rate-limit and low-pass filter the commanded torque (the
  `limit_rate` and `cutoff_frequency` arguments of `Robot::control`); the simulation does neither.
- **Reflexes.** The robot stops on its own collision and limit thresholds; MuJoCo only reports contacts.

What the hardware preflight still has to do (design doc, Hardware): simulate exactly the motion the robot will run
(same start, reference, gains, ramp) with `franka::Model` and the configured load, check every margin above including
the torque rate of τ_cmd, refuse to move if one fails, approach the start with a slow joint motion, and keep the
online safety monitor (damp to rest on any violation) in the torque loop.
