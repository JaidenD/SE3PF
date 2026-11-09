# Figures

The good runs, one folder per run. Folders 1 to 4 are the `loops` path (two loops of the coning circle at
50 mm/s), named after the controller and numbered in the order they build on each other. Folders 5 to 11 are
the runs of 2026-11-02 (see below and `labNotes.md`). The files inside keep
the run's name from `logs/` (e.g. `loops_14`). Each folder has the log (`<run>.csv`), the preflight prediction with the robot's model
(`<run>.csv.preflight.csv`) and two figures, each as PNG and PDF:

- `<run>_path`: the paper figure. (a) desired and actual end-effector path with body frames,
  (b) transverse error (position, attitude), (c) tangential tracking (path speed against its
  reference), (d) elbow deviation from the lift.
- `<run>`: the diagnostic figure: errors, beta, path speed, torques, joint margins, accelerations.

| Folder | Run | Controller | Friction model | Drive inertia | Path covered | Speed error RMS | Position / attitude error RMS |
|---|---|---|---|---|---|---|---|
| `1_PD_no-friction-model` | `loops_inertia_1` | PD | no | yes | 15 % (stalls) | 43 mm/s | 5.7 mm / 1.9 deg |
| `2_PD_with-friction-model` | `loops_measured_1` | PD | yes | yes | 76 % | 12.8 mm/s | 4.7 mm / 3.8 deg |
| `3_PID_no-friction-model` | `loops_pid_1` | PID | no | yes | 100 % | 7.8 mm/s | 6.7 mm / 4.5 deg |
| `4_PID_with-friction-model` | `loops_14` | PID | yes (Gaz et al. 2019 levels) | wrist only | 100 % | 2.8 mm/s | 3.3 mm / 0.9 deg |

RMS values are over the cruise at 50 mm/s. Caveat for `loops_pid_1`: the commanded joint
acceleration reached 59 % of the limit (wrist), above the usual 30 % safety margin, which was
relaxed in that build.

Regenerate with:

    python3 tools/plot_robot_path.py figures/<folder>/<run>.csv
    python3 tools/plot_robot_log.py figures/<folder>/<run>.csv

The CSVs are hard links to the files in `logs/` (and `important figs/` for `loops_measured_1`):
the same data stored once, so deleting one copy leaves the other.

## Runs of 2026-11-02

All with the PID controller and the measured drive inertia, no friction model, unless noted.

| Folder | Run | Experiment | Result |
|---|---|---|---|
| `5_coning_PID_stiff-speed-gains` | `loops_pid_hi_1` | loops-pid-hi: k_t 30, k_i,eta 100 | Full 38 s. Speed error 2.5 mm/s rms (k_t 15: 7.8; with the friction model: 2.8); pos / att 6.6 mm / 4.3° rms. |
| `6_orbit_velocity` | `orbit_velocity_1` | inspection orbit, two loops at 50 mm/s | Full run, no stop. Speed 100 %, \|xi\| ≤ 0.022, \|beta\| ≤ 0.091. |
| `7_orbit_position` | `orbit_position_1` | inspection orbit, one loop in six stations | Full run. \|xi\| ≤ 0.018, final 3.6e-3. |
| `8_dome_position` | `dome_position_1` | dome figure-eight, one loop in six stations | Full run. \|xi\| ≤ 0.017, \|beta\| ≤ 0.028, final \|xi\| 1.4e-3. |
| `9_selfmotion_hold_robot` | `orbit_selfmotion_hold_2` | self-motion back to the lift, end effector held (beta0 0.5) | beta 0.5 → −0.003. End effector off the path ≤ 4.6 mm / 3.3° during the swing (wrist friction), 0.6 mm at the end. |
| `10_selfmotion_simulation` | (simulated) | the same, exact model, no friction, PD | beta on the design curve; end effector ≤ 1.2 µm / 0.4 arcsec. See its README. |
| `11_rail_damper_push` | `rail_damper_3` | straight rail, pushed by hand (damper k_t 15) | Pushed about 7 cm along the rail within 3 mm of the line. A sideways push of 18 mm at 11.5 s reached the \|xi\| 0.03 stop, as designed. |

Figures in each folder: `<run>_path` (paper figure), `<run>` (diagnostics); folder 9 also `<run>_selfmotion`
(`tools/plot_self_motion.py`).

Not included:
- `dome_velocity_1`: torque-rate stop at 21.6 s, from velocity noise.
- `dome_velocity_lo_1`: softer gains, worse, and a stop.
- `orbit_selfmotion_hold_1`: beta stalled at 0.04 without the beta integral.
- `rail_damper_1` and `rail_damper_2`: damping too light or too heavy, ending in a speed stop and a reflex.
