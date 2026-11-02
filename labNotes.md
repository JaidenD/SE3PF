# Lab notes

A running record of robot experiments and what we learned from them. The newest day goes at the bottom. Logs are in
`logs/`; each run's `<log>.csv.preflight.csv` is the prediction with the robot's model. Plot with
`tools/plot_robot_path.py` (paper figure) and `tools/plot_robot_log.py` (diagnostics).

Units: eta, xi and the path speed are in the metric W = diag(0.01 I, I), so 1 rad counts as 0.1 m. "mm/s" of
path speed is metric mm/s, close to real mm/s on these paths.

---

## 2026-11-02

### Setup on this day

- **Robot:** Panda at 129.97.71.19, Franka Hand configured in Desk (0.73 kg). FCI must be activated in Desk first.
- **Plant model** (identification of 2026-10-20): (M(q) + J_d) q̈ + C(q, q̇) q̇ + G(q) + τ_f(q̇) = u.
  - J_d = diag(1.0, 1.0, 1.0, 1.0, 0.13, 0.12, 0.10) kg m², the drive inertia missing from franka::Model.
  - Coulomb friction F_c = (0.65, 0.40, 0.20, 0.30, 0.60, 0.30, 0.35) N m. Joints 2 and 4 are strongly direction
    dependent; that is not modelled.
- **Controller for the paper:** no friction model, J_d kept as part of M, integral action on (eta, xi, beta) to
  absorb friction. This is the `loops-pid` configuration: poles at -6, k_t = 15, K_i = K_i,beta = 36,
  k_i,eta = 25.
- **Safety limits — watch this.** `include/se3pf/panda.hpp` has uncommitted loosened limits:
  - velocity and commanded acceleration at 100 % (committed: 50 % and 30 %);
  - |xi| and |beta| at 1 (committed: 0.04 and 0.1).

  The torque (60 %) and torque-rate (600 N m/s) checks are unchanged. The unit tests do not compile with the
  loosened file; they were run on a copy with the committed limits. Several runs below exceed the 30 %
  acceleration margin.

### Runs

| Log | Experiment | Path | Controller | Outcome |
|---|---|---|---|---|
| `loops_pid_1` | loops-pid | coning circle | PID, no friction model | Full 38 s. Speed 100 % of nu_d (7.8 mm/s rms swing), pos/att error rms 6.7 mm / 4.5°, qdd* 59 % (joint 7). |
| `orbit_velocity_1` | orbit-velocity | inspection orbit | PID, no friction model | Full 39 s, no stop. Speed 100 %, \|xi\| ≤ 0.022, \|beta\| ≤ 0.091, torque rate ≤ 474 N m/s, qdd* 51 %. |
| `orbit_position_1` | orbit-position | inspection orbit | PID, no friction model | Full 43 s, no stop. \|xi\| ≤ 0.018, \|beta\| ≤ 0.069, final \|xi\| 3.6e-3, torque rate ≤ 416 N m/s, qdd* 49 %. |
| `dome_velocity_1` | dome-velocity | dome figure-eight | PID, no friction model | **Torque-rate stop at 21.6 s.** Tracking good until then: speed 102 %, \|xi\| ≤ 0.019, \|beta\| ≤ 0.050. |
| `dome_velocity_lo_1` | dome-velocity-lo | dome figure-eight | PID, paper's gains (poles -4, k_t = 5) | **Torque-rate stop at 22.2 s, worse tracking:** \|xi\| 0.041, \|beta\| 0.109, speed 6–173 mm/s. |
| `dome_position_1` | dome-position | dome figure-eight | PID, no friction model | Full 49 s, no stop. \|xi\| ≤ 0.017, \|beta\| ≤ 0.028, final \|xi\| 1.4e-3, torque rate ≤ 566 N m/s (close to 600), qdd* 53 %. |

Earlier runs used for comparison (same coning circle):

| Log | Controller | Speed (% of nu_d) | Speed error rms | Pos / att error rms |
|---|---|---|---|---|
| `loops_inertia_1` | PD, no friction model | 15 % (stalls) | 43 mm/s | 5.7 mm / 1.9° |
| `loops_measured_1` | PD + friction model | 76 % | 12.8 mm/s | 4.7 mm / 3.8° |
| `loops_pid_1` | PID, no friction model | 100 % | 7.8 mm/s | 6.7 mm / 4.5° |
| `loops_14` | PID + friction model | 100 % | 2.8 mm/s | 3.3 mm / 0.9° |

### Observations

**1. PID without a friction model keeps the arm moving, but friction shows as stick-slip at joint reversals**
(`loops_pid_1`).
- The integrator removes the mean speed lag of the PD (15 % → 100 %).
- A ±15–20 mm/s speed swing remains, repeating with path position: the two loops correlate 0.98 at the same
  path point.
- The two largest position spikes (16 mm at 11.6 s and 27 s) are joint 3 reversing.
- Integral action cancels only constant disturbances; Coulomb friction flips sign at every reversal.

**2. Without friction feedforward the wrist needs large commanded accelerations.**
- Friction can only be overcome through (M + J_d) q̈*, and the wrist's inertia is small (about 0.14 kg m²).
- So 0.6 N m of friction needs about 4 rad/s² (≈ 30 % of the limit) before the joint moves.
- Every PID run without friction model reached qdd* 45–59 %. With the committed 30 % limit they would all
  stop (`loops_integral_1` did, at 2.6 s).

**3. Raising the tangential gains** was prepared (`loops-pid-hi`: k_t 30, k_i,eta 100) but not run.

**4. The new paths are much richer in orientation than the coning circle.**

| Path | Rotation over a loop | Joint margin | Tube radius | Peak joint speed at 50 mm/s |
|---|---|---|---|---|
| Coning circle | 29° | 0.36 rad | 0.125 | 7 % |
| Inspection orbit | 52° | 0.41 rad | 0.093 | 12.5 % |
| Dome figure-eight | 92° | 0.34 rad | 0.053 | 14.4 % |

- The dome was meant as a tennis-ball seam. A full seam would need the tool pointing up, so it became a
  figure-eight over a dome. It twists ±46° about the tool axis so it does not cross itself in SE(3).
- The dome has no lift with `hardwareLiftWeights()` (weight 10 on joint 7); the driver uses the unweighted lift.

**5. The torque-rate stops are noise, not motion** (`dome_velocity_1`).
- Velocity noise passes through the derivative gains (K_d = 12, k_t = 15) and joint 2's large command inertia
  (M_22 + 1.0 ≈ 2.7 kg m²).
- Joint 2's command jitters by about 0.09 N m rms per ms, i.e. 100–150 N m/s; the correlation with qdd*₂ noise is
  0.94.
- Orbit and dome have the same noise floor (median about 90 N m/s, 99th percentile about 285–290). Both reach
  450–566 N m/s regularly, so the orbit runs passed on margin, not because the orbit is quieter.
- The dome runs a bit hotter on the half of its loop with more wrist reversals and larger planned joint 7
  accelerations (99th percentile about 330–350 N m/s).
- libfranka smooths the command (100 Hz low-pass, rate limiter) before the robot, but our check looks at the raw
  command.

**6. Softer gains made the dome worse** (`dome_velocity_lo_1`).
- Torque-rate noise did drop (mostly 150–300 N m/s), but friction rejection collapsed:
  - joints stuck longer and the path speed swung from 6 to 173 mm/s;
  - the integrals hit their bounds (int_eta −1.0, int_beta 3.0);
  - at about 21 s the stored effort released: joints 2, 4 and 6 lunged at about 1 rad/s and the torque rate
    tripped.
- Conclusion: keep the gains and deal with the noise elsewhere.

**7. Elbow-deviation spikes in the position runs come from sticky wrist joints, not elbow joints.**
- When a wrist joint sticks, the controller keeps the flange on the path with the shoulder and elbow joints,
  which moves the elbow. When the joint breaks free, the arm snaps back.
- Nearly every spike is 2.5–3.5 s into a station move, when a wrist joint's planned speed passes near zero.
- Share of move time stuck (|q̇| < 0.01 while the plan asks > 0.02 rad/s):

  | Run | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
  |---|---|---|---|---|---|---|---|
  | orbit_position_1 | 8 % | 1 % | 3 % | 4 % | 16 % | 14 % | 34 % |
  | dome_position_1 | 3 % | 2 % | 5 % | 3 % | 19 % | 8 % | 20 % |

- **Sideways swings** (largest): orbit 11.2 s (11.4 mm), dome 27.0 s (10.3 mm).
  - Joints 5 and/or 7 lag by 0.14–0.17 rad; joints 1 and 3 compensate in opposite directions (−0.05 / +0.05 rad).
  - This is self-motion, so |beta| is large (−0.066 on the orbit). A planned reversal of joint 1 or 7 is often
    at the same moment.
- **Up-and-down swings:** orbit 29.4 s and 35.5 s, dome 3.1 s, 12.5 s and 42.6 s.
  - Joint 6 (wrist pitch) sticks; joints 2 and 4 (parallel axes) compensate, with joint 2 about +0.03 rad off plan.
  - Joint 2 itself follows its planned speed.
- Possible fixes: friction compensation on joints 5–7 only, a stiffer beta loop (for the sideways swings), or
  paths that keep wrist reversals out of the moves.

**8. Self-motion demo prepared** (`orbit-selfmotion-hold`, `orbit-selfmotion-moving`; not yet run on the robot).
- Start: flange on the orbit at s0 (xi0 = 0), arm off the lift by beta0 = 0.3 rad. Joints 1 and 3 are
  counter-rotated by about 11° and 13°; the start pose is within 0.03 mm of gamma(s0).
- Variants: hold eta = s0 for 8 s, or follow one loop at 50 mm/s.
- The elbow point moves only about 10 mm for beta0 = 0.3 (17 mm for 0.5): on the orbit the upper arm is nearly
  vertical, so the elbow sits close to the axis it swings around. The self-motion is visible in joints 1 and 3.
- The beta integral overshoots. In the preflight (exact model) with K_i,beta = 36, beta crossed zero and
  overshot by 46–68 % of beta0: the integral winds up while beta is large. With K_i,beta = 0 there is no overshoot
  and beta follows beta0 (1 + 6t) e^{-6t}, so the demos run without the beta integral; the eta and xi integrals
  stay. On the robot, expect a small residual beta from friction on joints 1 and 3.
- Preflight and mock-robot runs pass:
  - hold: |xi| ≤ 4e-4 (0.35 mm), qd 18 %, qdd* 14 %;
  - moving: joint margin 0.41 rad.
- Needs the beta bound above 0.3 rad. With the committed HardwareBounds::beta = 0.1 the monitor stops the run at
  once; the loosened file currently allows it.
- `plot_robot_path.py` fix: the "elbow on the lift" reference now uses only preflight samples with |beta| < 1e-4.
  Before, a run that starts off the lift got a wrong reference.

**9. First self-motion run** (`orbit_selfmotion_hold_1`): full 8 s, no stop, but friction spoiled the demo.
- beta followed the plan until about 0.75 s (0.235 vs 0.229 at 0.5 s), lagged by 1 s (0.088 vs 0.033), then
  stuck at about 0.04 from 2 s and was still 0.029 at 8 s (plan: 0).
- Joints 1 and 3 stopped about 0.024 and 0.027 rad short of the lift: static friction beats the small beta
  feedback near zero, and there is no beta integral to build up.
- Wrist joints 5 and 7 stuck about 0.03 rad off plan during the swing, so the flange moved 4.8 mm and 2.4° (plan:
  0.35 mm). The xi integral pulled it back slowly (int_xi 0 → 0.88).
- At 6.8 s joint 7 broke free (−0.08 rad/s), overshot by 0.06 rad, and |xi| jumped back to 6.3e-3.
- Figure: the path plot does not suit a hold. Panel (a) zooms into a 5 mm box with 30 mm frame axes, and panel (d)
  hides the 10 mm elbow motion under beta in mrad. A dedicated self-motion figure is needed (joints 1 and 3 vs
  the lift, beta vs the design, flange error).
- Next: a beta integral with a small bound (break joints 1 and 3 free with little overshoot), beta0 = 0.5, a
  self-motion figure.

**10. Self-motion demo retuned** (not yet run on the robot).
- beta0 0.3 → 0.5 rad: joints 1 and 3 swing back about 19° and 21°.
- beta poles −6 → −4 (K_p,beta 16, K_d,beta 8). At beta0 = 0.5 with poles at −6, joint 3 would reach 31 % of
  its speed limit, and the preflight refused it.
- Capped beta integral: K_i,beta = 16, bound 1.5.
- Tuned on the mock robot with the identified Coulomb friction and drive inertia
  (`FRANKA_MOCK_COULOMB=0.65,0.40,0.20,0.30,0.60,0.30,0.35 FRANKA_MOCK_ROTOR=1,1,1,1,0.13,0.12,0.10`). The mock
  reproduces `orbit_selfmotion_hold_1` qualitatively: beta stalls at 0.02 (robot 0.04), |xi| ≤ 0.0065
  (robot 0.0063). Real friction looks about twice the mock's.
- Sweep at beta0 = 0.5, poles −4 (final beta; undershoot where there is any):

  | beta integral (gain / bound) | no friction | 1× | 1.5× | 2× |
  |---|---|---|---|---|
  | off | 0.000 | 0.053 | 0.076 | – |
  | 16 / 0.5 | 0.000 (5 % under) | 0.022 | 0.046 | – |
  | 16 / 1.0 | 0.000 (10 %) | −0.007 | 0.016 | – |
  | **16 / 1.5** | 0.000 (15 %) | **0.001 (6 %)** | −0.012 (3 %) | 0.004 |
  | 16 / 2.0 | 0.000 (19 %) | 0.002 (11 %) | – | −0.012 (4 %) |
  | 36 / 1.0 | 0.000 (9 %) | 0.002 | 0.016 | – |

- Chosen: 16 / 1.5. On the mock with 1× friction:
  - hold: beta 0.001 at the end, undershoot 6 %, max |xi| 0.0081 (5.3 mm / 3.8°), end |xi| 0.0003, joint 3 at
    26 % of its speed limit;
  - moving: undershoot 16 %, max |xi| 0.0135.
- The flange disturbance (about 5 mm / 4°) is the wrist joints 5 and 7 sticking while they counter-rotate. The
  beta integral does not change it; only wrist friction compensation would.
- New `tools/plot_self_motion.py`, three stacked panels (the 3D arm panel was dropped at the operator's request):
  - (a) beta against the prediction;
  - (b) joints 1, 3, 5 and 7 minus their lift values;
  - (c) the end-effector error.

**11. Stiffer end-effector gains against the wrist disturbance** (mock, hold demo, beta0 0.5, beta poles −4).

| xi poles (K_p / K_d / K_i) | friction | max end-effector error | final beta | qdd* |
|---|---|---|---|---|
| −6 (36 / 12 / 36) | 1.5× | 6.8 mm / 4.7° | −0.012 | 23 % |
| −8 (64 / 16 / 64) | 1.5× | 5.0 mm / 3.4° | −0.009 | 24 % |
| −10 (100 / 20 / 100) | 1.5× | 3.7 mm / 2.6° | −0.009 | 26 % |
| −10 | 1× | 2.9 mm / 2.0° (−6: 5.3 mm / 3.8°) | +0.003 | 19 % |

- The disturbance shrinks roughly in proportion to the stiffness, but does not vanish: the wrist still has to
  break away.
- beta poles at −2 do not work with K_i,beta = 16: s³ + 4s² + 4s + 16 is marginally stable, so beta never
  settles. That would need its own retune.
- The limit is noise. The mock has no velocity noise, so the torque rate is estimated by scaling
  `orbit_selfmotion_hold_1` (K_d 12: 99.9th percentile 256, max 349 N m/s) with K_d: about 465 max at −8, 523
  at −9, 581 at −10, against the 600 N m/s limit.

- Applied: transverse poles at −8 (K_p 64, K_d 16, K_i 64) in both self-motion demos. On the mock with 1×
  friction:
  - hold: beta 0.5 → −0.0004 (undershoot 6 %), end effector off the path up to 4.2 mm / 2.6° (was 5.3 mm / 3.8°),
    0.6 mm at the end;
  - moving: beta → 0.000 (undershoot 15 %), up to 4.4 mm / 3.6° (was |xi| 0.0135), 0.0 mm at the end.

**12. Self-motion demo in simulation** (`figures/selfmotion_sim/`).
- Exact model, no friction, PD only (no integrals), no controller fade-in; transverse poles −8, beta poles −4.
- beta follows beta0 (1 + 4t) e^{−4t} to within 4e-4 rad; the joints are on the lift by about 2 s.
- The end effector stays on the path to within 1.2 µm and 0.4 arcsec (|xi| ≤ 1.2e-6), hold and moving alike.
- With the robot's 1 s fade-in left on, beta lags the design by about 0.3 s. The fade-in is a hardware measure,
  so it is off in the simulation.
- `plot_self_motion.py` gained `--design LAMBDA` (overlay the PD design) and `--sim`. The error panel's scale is
  at least 1 mm / 1°, so that zero reads as zero.
- Panel (c) with `--sim` is on a log scale and shows the elbow's distance from its lift position (mm) for scale:
  the elbow starts about 17 mm off the lift, the end effector's error peaks at about 1 µm and 0.4 arcsec, roughly
  four orders of magnitude apart. (A linear axis showed an empty panel; |beta| in rad next to µm was misleading.)
- The difference from the robot run (about 4–5 mm, 3–4°) is entirely friction, mostly wrist joints 5 and 7
  sticking.

**13. Rail experiments for pushing by hand prepared** (`rail-damper`, `rail-spring`; not yet run on the robot).
- Path `PathKind::Rail`: a smooth racetrack at constant attitude (tool down).
  - Centre (0.55, 0, 0.35) m, half-length 25 cm, half-width 10 cm, sharpness 2.
  - The near side is straight to 1.4 mm over ±10 cm and 4.4 mm over ±15 cm; s0 is its middle.
  - Tube radius 0.057 (with sharper, narrower turns: 0.017); joint margin 0.57 rad.
- Control:
  - transverse poles −8, beta poles −4, **no integrals** (a held push would wind them up and throw the arm back
    on release);
  - along the rail, `rail-damper` has nu_d = 0 (mu_eta = −k_t eta_dot, k_t = 10), and `rail-spring` has
    eta_d = s0 (k_p,eta 16, k_d,eta 8).
- Limits: |xi| 0.03, flange 0.15 m/s, joint speed 30 %, |beta| 0.1. The robot's collision reflex is at 20 N at
  the end effector.
- Mock with 1× friction. The mock's push is a torque on one joint: joint 1 ≈ along the rail, joint 4 ≈ across.

  | Push | Result |
  |---|---|
  | Across (4 N m on joint 4, ≈ 13 N) | 14 mm off the line (≈ 1 kN/m), back to < 1 mm on release, beta untouched |
  | Along, damper, k_t 4 | 140 mm/s, ran 255 mm; beta stop (see below) |
  | Along, damper, k_t 10 | ≈ 57 mm/s, stops within about 1 cm of release; about 2° attitude left after (wrist friction, no integral) |
  | Along, spring, 4 / 4 | friction held it 41 mm short of the start |
  | Along, spring, 16 / 8 | pushed 42 mm, returns to about 10 mm |

- A torque on joint 1 also drives the self-motion (beta to −0.09, a stop at 0.1 with the soft damper). A hand
  force at the end effector cannot drive beta directly: X^T J^T F = (J X)^T F = 0.

**14. Robot runs: self-motion (2nd) and first rail push.**
- `orbit_selfmotion_hold_2` (beta0 0.5, xi poles −8, beta integral 16 / 1.5): worked.
  - beta 0.5 → −0.003; joints back on the lift (joint 1 +0.000, joint 3 −0.001 rad from q*).
  - End effector off the path up to 4.6 mm / 3.3° during the swing, 0.6 mm / 0.7° at the end.
  - Torque rate peak 470 N m/s (estimated beforehand: about 465). Finished at rest.
- `rail_damper_1`: our flange-speed stop, then the robot's Cartesian reflex. Both safety layers worked; nothing
  was wrong with the controller.
  - The operator pushed the end effector 124 mm along the rail at up to 136 mm/s, then reversed sharply to
    166 mm/s. The joint torque sensors show the hand needed under 1 N m: k_t = 10 is far too light a damper for a
    natural push.
  - At 3.565 s the flange speed passed the 0.15 m/s limit, so the safety stop held the arm (stiff damped hold).
  - The hand kept pushing against the held arm. The robot's own external-torque estimate rose to about 9 N m on
    joints 1 and 3 (about 20 N at the end effector), and libfranka aborted with `cartesian_reflex` at 3.694 s.
  - Afterwards the robot stays in mode Reflex until `robot_driver <ip> error-recovery`.
- Fix: heavier damping along the rail. Tell the operator to let go as soon as the arm stops.
- Applied: `rail-damper` k_t 10 → 30, and the watch text now says to push slowly, avoid sharp reversals and let
  go if the arm stops. On the mock the same push (4 N m on joint 1) now moves the end effector at about 20 mm/s
  (k_t 10: 57 mm/s), and it stops within about 4 mm of release.

**15. `rail_damper_2` (k_t 30): ended by the robot's Cartesian reflex at 11.4 s; no safety stop of ours.**
- The operator pushed the end effector along the rail −116 mm, then +103 mm, then back, at up to 59–89 mm/s.
  The transverse error stayed ≤ 0.012 and the flange under 0.093 m/s.
- At the end it was moving back at 89 mm/s. The robot's external-torque estimate (tau_ext_hat_filtered, last 50
  states) is 8.5 / 8.3 N m on joints 1 and 3, which is **19.1 N along the rail at the flange**: the 20 N
  Cartesian collision threshold.
- That is the damper's own force: 19 N / 0.089 m/s = 214 N s/m, i.e. an effective mass of about 7 kg × k_t 30.
  The operator's hand was supplying it.
- (A quick estimate from tau_J − tau − G gave 1–4 N. That is wrong: it ignores the drive inertia and the arm's
  acceleration. The robot's estimate is the one to trust.)
- So the damper sets a speed–force trade-off: force ≈ 7 kg · k_t · speed.

  | k_t | 0.15 m/s speed stop reached at | 20 N reflex reached at |
  |---|---|---|
  | 10 | 10.7 N | 0.28 m/s |
  | 30 | 32 N | 0.09 m/s |

  k_t 10 hit the speed stop first; k_t 30 hit the reflex first.
- Choice: about k_t 15 (107 N s/m). The speed stop then comes first, at about 16 N; a natural push of about 10 N
  moves it at about 9 cm/s. Alternatively raise the Cartesian reflex threshold for this experiment (a safety
  decision for the operator).
- Applied: `rail-damper` k_t = 15.

**16. Two more runs reviewed, and the good runs copied to `figures/`.**
- `loops_pid_hi_1` (loops-pid-hi, k_t 30, k_i,eta 100, coning circle): full 38 s.
  - Speed error 2.5 mm/s rms (loops_pid_1: 7.8; loops_14 with the friction model: 2.8), pos / att 6.6 mm / 4.3°
    rms (about unchanged).
  - qdd* 50 %, torque rate up to 551 N m/s (close to 600).
  - The stiffer tangential loop removes most of the speed swing that friction causes.
- `rail_damper_3` (k_t 15): pushed about 7 cm along the rail at up to about 42 mm/s, within 3 mm of the line.
  - At 11.5 s a push across the rail took it 18 mm off the line, reaching the |xi| 0.03 stop at 11.7 s, as
    designed.
  - The stop ran its full 2 s without settling (|dq| 0.12 rad/s at the end, 0.05 rad drift): a hand still on
    the arm.
- `figures/` now holds runs 5–11 (see its README); `selfmotion_sim` was renamed `10_selfmotion_simulation`.

### Code and tools changed today (not committed)

- `tools/plot_robot_log.py`:
  - the transverse error is split into position (mm) and attitude (deg);
  - the tangential panel shows the output against its reference;
  - 8 panels.
- `tools/plot_robot_path.py`:
  - paper-style figure: path with body frames, transverse error, tangential tracking (η̇ vs ν_d or η vs η_d),
    elbow deviation;
  - PNG and PDF output.
- `figures/`: the four coning-circle comparison runs, with figures and data (hard links into `logs/`).
- New experiments: `loops-pid`, `loops-pid-hi`, `orbit-velocity`, `orbit-position`, `dome-velocity`,
  `dome-position`, `dome-velocity-lo`, `orbit-selfmotion-hold`, `orbit-selfmotion-moving`, `rail-damper`, `rail-spring`.
- `SampledPath` in `se3pf/path.hpp`: any smooth closed pose curve, reparameterized numerically by arc length.
  Also `inspectionOrbit` and `domeFigureEight`.
- `PathKind` and `PandaExamples`: each experiment names its path, and the driver builds the matching example.
- `StationSchedule`: a position reference that moves through stations with stops.
- Logs carry the desired pose gamma(eta) (columns gx..grz), so the plots work for any path.

### Open questions

- Restore the committed safety limits, or confirm the loosening is intended?
- Torque-rate check: check the command after libfranka's filter (what the robot receives) instead of the raw
  command? This needs the operator's OK.
- For the paper: present PID without a friction model, or add friction compensation (at least on the wrist)?
