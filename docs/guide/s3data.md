## Architecture data for the diagram

Layers bottom-up, every module in exactly one layer; edges from the `#include "..."` lines, then the calls that cross module boundaries. A -> B means A uses B.

### Layers and modules

- foundation
  - `types.hpp`: Eigen fixed-size aliases, `Mat<R, C>`, `Vec<N>` (header only)
  - `se3.hpp` / `se3.cpp`: hat, vee, exp/log SO(3), J_r, Ad, ad
  - `spline.hpp`: periodic cubic Hermite spline template (header only)
- geometry
  - `metric.hpp` / `metric.cpp`: product metric W: Exp_I, Log_I, Γ, distance
  - `path.hpp` / `path.cpp`: `Path` interface; coning and tumbling circles
  - `tube.hpp` / `tube.cpp`: normal frame N, Ψ, projection, J_Ψ, radius
  - `rigid_body.hpp` / `rigid_body.cpp`: rigid-body System; Euler equations, RK4 integrator
  - `panda.hpp` / `panda.cpp`: Panda System: DH kinematics, J_b, limits, hardware bounds
- outputs
  - `redundancy.hpp`: `RedundancyOutput` interface for β, dβ (header only)
  - `lift.hpp` / `lift.cpp`: lift q*, kernel frame X, β outputs, starts
- control
  - `controller.hpp`: `Controller<System>`, gains, references, ramp (header only)
- examples
  - `rigid_body_example.hpp` / `.cpp`: paper's rigid-body setup and closed loop
  - `panda_example.hpp` / `.cpp`: paper's Panda setup and closed loop
- apps / sim (not in the library target `se3pf`)
  - `apps/rigid_body_sim.cpp`: runs rigid-body example, writes CSV
  - `apps/panda_sim.cpp`: runs Panda example, writes CSV
  - `apps/panda_limits.hpp`: writes joint limits CSV (header only)
  - `sim/mujoco_panda.hpp` / `.cpp`: MuJoCo plant, torque-level closed loop (target `se3pf_mujoco`)
  - `apps/panda_mujoco.cpp`: MuJoCo kinematics check (0) and verification runs A to D
  - `apps/mujoco_smoke.cpp`: raw MuJoCo API checks, no se3pf

Placement notes: the two systems sit in geometry because they include only `types.hpp` and use only `se3`, and `lift` (outputs) uses `panda`. `tube` sits in geometry because `lift` uses `Tube::Psi`, so it must be below outputs.

### Header include edges (from the .hpp files)

- `se3.hpp` -> `types.hpp`
- `metric.hpp` -> `types.hpp`
- `path.hpp` -> `types.hpp`, `metric.hpp`
- `tube.hpp` -> `types.hpp`, `metric.hpp`, `path.hpp`, `spline.hpp`
- `rigid_body.hpp` -> `types.hpp`
- `panda.hpp` -> `types.hpp`
- `redundancy.hpp` -> `types.hpp`
- `lift.hpp` -> `types.hpp`, `panda.hpp`, `path.hpp`, `redundancy.hpp`, `spline.hpp`, `tube.hpp`
- `controller.hpp` -> `types.hpp`, `redundancy.hpp`, `tube.hpp`
- `rigid_body_example.hpp` -> `types.hpp`, `controller.hpp`, `metric.hpp`, `path.hpp`, `rigid_body.hpp`, `tube.hpp`
- `panda_example.hpp` -> `types.hpp`, `controller.hpp`, `lift.hpp`, `metric.hpp`, `panda.hpp`, `path.hpp`, `redundancy.hpp`, `tube.hpp`
- `spline.hpp` -> standard library only; `types.hpp` -> `<Eigen/Dense>`

### Source-only include edges (from the .cpp files; implementation detail, not visible to users of the header)

- `metric.cpp`, `path.cpp`, `tube.cpp`, `rigid_body.cpp`, `panda.cpp`, `lift.cpp`, `rigid_body_example.cpp`, `panda_example.cpp` -> `se3.hpp`
- `rigid_body.cpp` -> `<Eigen/Cholesky>`

### Call edges that cross modules

- `metric` -> `se3` (expSO3, logSO3, Jr, JrInv, makePose, inverse)
- `path` -> `se3` (expSO3, makePose); `path` -> `metric` (norm, in the constant-speed check)
- `tube` -> `path` (virtual pose, velocity, acceleration, length)
- `tube` -> `metric` (Exp, Log, Gamma, GammaInv, inner, flat, norm, distance, matrix)
- `tube` -> `spline` (`PeriodicHermite<Mat65>`); `tube` -> `se3` (Ad, inverse, logSO3)
- `rigid_body` -> `se3` (expSO3, JrInv, makePose)
- `panda` -> `se3` (Ad, inverse, makePose)
- `lift` -> `panda` (PandaArm::output, jointMargin, restPosture), `path` (pose, length), `tube` (Tube::Psi in startConfiguration), `spline` (`PeriodicHermite<Vec7>`), `se3` (logSO3, inverse)
- `lift` -> `redundancy` (implements it: `LiftRedundancy`, `CouplingRedundancy`)
- `controller` -> `tube` (outputs, JPsi, wrap, length), `redundancy` (virtual evaluate)
- `controller` -> System template parameter (`System::output`, `System::retract`): `rigid_body` or `panda`, bound at instantiation, no include
- `rigid_body_example` -> `controller` (`Controller<RigidBody>`), `rigid_body` (mass, coriolis, integrate), `tube`, `path` (`TumblingCircle`), `metric`, `se3`
- `panda_example` -> `controller` (`Controller<PandaArm>`), `lift` (Lift, KernelFrame, LiftRedundancy, CouplingRedundancy, startConfiguration), `panda`, `tube`, `path` (`ConingCircle`), `metric`, `se3`
- `apps/rigid_body_sim.cpp` -> `rigid_body_example`
- `apps/panda_sim.cpp` -> `panda_example`, `apps/panda_limits.hpp`
- `apps/panda_limits.hpp` -> `panda`
- `sim/mujoco_panda` -> `types`, `controller`, `panda`, `panda_example`, MuJoCo C API
- `apps/panda_mujoco.cpp` -> `sim/mujoco_panda`, `apps/panda_limits.hpp`, `panda`
- `apps/mujoco_smoke.cpp` -> MuJoCo C API only

### Build targets (CMakeLists.txt)

- `se3pf` (static library): the nine .cpp files of `src/`; header-only modules: `types`, `spline`, `redundancy`, `controller`
- `rigid_body_sim`, `panda_sim` -> `se3pf`
- `tests` (Catch2 v3.5.2) -> `se3pf`
- with `-DSE3PF_WITH_MUJOCO=ON`: `se3pf_mujoco` -> `se3pf` + `mujoco`; `panda_mujoco` -> `se3pf_mujoco`; `tests_mujoco` -> `se3pf_mujoco`; `mujoco_smoke` -> `mujoco` only

### Data that flows offline -> online (for arrows between phases)

- `ConingCircle` -> `Tube` (N, N′ at 2001 nodes) -> online `Tube::frame`
- `ConingCircle` -> `Lift` (q*, q*′ at 2001 nodes) -> `KernelFrame` (X, X′ at the same nodes) -> online `LiftRedundancy::evaluate`
- `Tube` + `LiftRedundancy` -> `Controller<PandaArm>` (references) -> q̈* -> caller (integrator, MuJoCo torque, robot torque)
