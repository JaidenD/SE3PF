# SE(3) Path Following in C++: Design Document

> A snapshot of the design document, taken on 2026-10-13 when milestones 1 to 10 were done. The living
> version is the artifact https://claude.ai/code/artifact/d6354d3d-b74d-42e8-bae0-b8a6e5c707e4. The code
> comments cite its sections ("design doc, Online: projection"). Where the code and this document differ,
> the code and its tests are right: see README.md and TICKETS.md for what was decided along the way (for
> example, the lift uses the projected Newton step, and Metric accepts the product type only).

2026-10-12 · @Jaiden

## Purpose and scope

The project implements the controller of Pose Path Following using Normal Tubular Coordinates in C++17 with Eigen, and reproduces its two simulations: the free rigid body and the seven-joint Franka Panda. Hardware (libfranka 0.9.2) comes last, after both simulations match the paper.

- In scope: the product-type metric W = diag(w_R I, w_p I) with closed-form Exp_I, Log_I and Γ; constant-speed paths; the normal frame, tubular map and outputs (Section IV); the lift q* and redundancy output β for one redundant joint; the controller of Section V; the rigid-body and Panda simulations.
- Optional, later: a general W (Exp_I and Γ from the geodesic ODEs), numerical arc length, the libfranka driver and its safety checks.
- Out of scope: torque-level robustness studies with identified Panda dynamics.

This doc fixes interfaces, conventions, algorithms and tests; every function body is yours to write. Each module section lists what to declare, the paper equations to implement, the pitfalls, and the tests that show it is right. The earlier implementation in ~/Desktop/panda_path_following is a numerical oracle: compare your numbers with it, and read its code only after writing your own.

## Conventions

Fix these once and never mix them: most silent bugs in this kind of code are a swapped twist order or a frame mix-up.

| Quantity | Convention |
| --- | --- |
| Twist ζ | (ω, v), angular first, as in the paper. Murray–Li–Sastry and libfranka use (v, ω). |
| Velocity of a pose | Body coordinates: ζ = (g⁻¹ ġ)^∨, so ω = (RᵀṘ)^∨ and v = Rᵀṗ. |
| Pose g = (R, p) | Eigen::Matrix4d, homogeneous; compose by matrix product; invert in closed form. |
| Body Jacobian J_b | 6 × n, rows (ω, v) in the tool frame: ζ = J_b(q) q̇. |
| Units | SI: m, rad, s. ξ is in “metres-equivalent”: one radian counts as ℓ = √(w_R / w_p) metres. |
| Metric | W = diag(w_R I, w_p I); Panda example w_R = 0.01, w_p = 1 (ℓ = 0.1 m); rigid body w_R = w_p = 1. |
| Path parameters | Curves use θ ∈ [0, 2π), the paper's curve parameter; the path uses arc length s ∈ [0, L) in the metric W; constant speed σ gives s = σθ. λ is the transverse gain (poles at −λ). |
| Tangential output | η wraps into [0, L). The reference η_d is real-valued; e_η is wrapped into (−L/2, L/2]. |
| Indexing | Zero-based in code: ξ₁ is xi[0], joint 1 is q[0], matrix entry a₁₂ is A(0, 1). |
| Panda | Tool frame = flange (Franka frame F). Rest posture q_r = (0, −π/4, 0, −3π/4, 0, π/2, π/4). |
| Names | Namespace se3pf. Fixed-size aliases Vec5, Vec6, Vec7, Mat5, Mat6, Mat65 (N), Mat67 (Panda J_b), Pose, and two alias templates for sizes that depend on the system: template \<int R, int C> using Mat = Eigen::Matrix\<double, R, C>; and template \<int N> using Vec = Eigen::Matrix\<double, N, 1>;. Sizes are int, as in Eigen. The systems are RigidBody and PandaArm. |

### Eigen pitfalls

- No auto for Eigen results. auto y = A * dq; stores an unevaluated product, computed later from whatever A holds then; a function returning auto can hand back an expression pointing into a destroyed local. Write the type: Vec7 y = A * dq;. No compiler warning catches this.
- Nothing starts at zero. Vec7 a; holds garbage. Use ::Zero() or .setZero() for anything you fill conditionally or accumulate into.
- Aliasing. M = 0.5 * (M - M.transpose()); overwrites entries before reading them. Write M = (0.5 * (M - M.transpose())).eval(); or assign to a new variable. Products such as N = N * E are safe.
- Comparisons. == on Eigen objects is exact, and a.isApprox(b, tol) is relative, so it never passes against zero. Use (a - b).cwiseAbs().maxCoeff() < tol, or a.isZero(tol) when the expected value is zero.
- Alignment. With C++17 and Eigen 3.4, fixed-size members and std::vector\<Mat65> need no alignment macros or allocators.

## Architecture

The code is a library of four layers (foundation, geometry, outputs, control) with small programs, the apps, on top; each layer uses only the layers below it, and nothing in the library knows about libfranka, files or the clock.

*(Diagram in the living document: the module layers, five levels and twelve modules.)*

Work is split by when it happens. Everything expensive runs once, offline, before the first control cycle: the normal frame, the lift and the splines that store them. A control step evaluates the outputs three times, at q and at the two points of the central difference for Ȧq̇. Each evaluation is one projection (usually one Newton update), one 6 × 6 factorization of J_Ψ and a few spline lookups; then comes one n × n solve.

```
SE3_pathFollowing/
  CMakeLists.txt
  include/se3pf/   types.hpp se3.hpp metric.hpp path.hpp spline.hpp tube.hpp
                   rigid_body.hpp panda.hpp redundancy.hpp lift.hpp controller.hpp
  src/             one .cpp per non-template header
  apps/            rigid_body_sim.cpp panda_sim.cpp (robot_driver.cpp last)
  tests/           test_<module>.cpp, one per module
  tools/           plotting scripts
```

- Move your header first. It is at se3pf/se3.hpp in the project root, but CMakeLists.txt searches include/, so #include "se3pf/se3.hpp" fails until you run mkdir -p include && mv se3pf include/.
- Growing CMakeLists.txt: add a module's .cpp to add_library(se3pf …) and its test to add_executable(tests …) only when they exist. Each app is add_executable(rigid_body_sim apps/rigid_body_sim.cpp) plus target_link_libraries(rigid_body_sim PRIVATE se3pf). Give tests and the apps -Wall -Wextra too, since templates compile where they are used.
- Two builds: build/ as Debug for unit tests (Eigen's size and range checks on; optionally -fsanitize=address,undefined), and build-rel/ with -DCMAKE_BUILD_TYPE=Release for simulations, the robot driver, offline computations and timing. Unoptimized Eigen is 30–100× slower, so every time in this doc assumes Release. Your current build/ is Debug.
- Git: SE3_pathFollowing is not its own repository yet; it sits untracked inside the one at ~/Desktop. Run git init in it and ignore build/ and build-rel/.

## Paper-to-code map

Every paper object has one home and one signature; the rest of this doc assumes these interfaces. Equation names are the \labels in your .tex, so you can search them. Signatures are suggestions: change them, but update this table when you do.

| Paper object | Where in the paper | Module | Suggested declaration |
| --- | --- | --- | --- |
| hat, vee on R³ | Notation | se3 | Matrix3d hat(const Vector3d&), Vector3d vee(const Matrix3d&) |
| ζ̂ and its inverse on R⁶ | Notation | se3 | Matrix4d hat6(const Vec6&), Vec6 vee6(const Matrix4d&) |
| Ad_g, ad_ζ | Notation | se3 | Mat6 Ad(const Pose&), Mat6 ad(const Vec6&) |
| exp, log on SO(3), J_r | Example 1 | se3 | expSO3, logSO3, Matrix3d Jr(const Vector3d&), JrInv |
| ⟨a, b⟩_W, ‖a‖_W | Notation | metric | double Metric::inner(const Vec6&, const Vec6&) const |
| Exp_I, Log_I, Γ | eq:Gamma, Example 1 | metric | Pose Metric::Exp(const Vec6&) const, bool Metric::Log(const Pose& g, Vec6& z) const, Mat6 Metric::Gamma(const Vec6&) const, Mat6 Metric::GammaInv(const Vec6&) const |
| dist | eq:dist | metric | double Metric::distance(const Pose&, const Pose&) const |
| γ(s), ζ_γ(s), ζ_γ′(s), L | eq:path-tangent | path | Pose pose(double s), Vec6 velocity(double s), Vec6 acceleration(double s), double length() |
| N, q*, X between grid points | Implementation (cubic Hermite) | spline | void PeriodicHermite\<T>::eval(double s, T& y, T* dy) const |
| N(s), N′(s), H, Φ | Lemma 1, eq:transport, eq:holonomy, eq:frame-ode | tube | void Tube::frame(double s, Mat65& N, Mat65* dN) const |
| Ψ(s, ξ) | eq:Psi | tube | Pose Tube::Psi(double s, const Vec5& xi) const |
| e_g(s), φ_g(s), φ_g′(s) | eq:projection-error, eq:projection-condition, Online stage | tube | bool Tube::projectionCondition(const Pose& g, double s, double& phi, double& dphi, Vec6& e) const |
| (η, ξ) = Ψ⁻¹(g) | Proposition 1 | tube | bool Tube::outputs(const Pose& g, double guess, double& eta, Vec5& xi, int* iterations = nullptr) const |
| Initial η (closest point, for starts) | Proposition 1(c) | tube | bool Tube::globalProject(const Pose& g, double& eta, int samples = 720) const |
| J_Ψ(y) | body Jacobian of Ψ (cited as eq:JPsi in Algorithm 1) | tube | Mat6 Tube::JPsi(double eta, const Vec5& xi) const |
| ε̂ | Implementation, Tube radius | tube | double estimateTubeRadius(const Tube&, const RadiusOptions&) |
| Π(q), J_b(q) | Problem formulation, Systems | systems | static Pose System::output(const Config& q, Mat<6, System::kDof>* Jb) |
| M, Cζ, G (rigid body) | eq:EL, Free rigid bodies | rigid_body | Mat6 RigidBody::mass() const, Vec6 RigidBody::coriolis(const Vec6& zeta) const |
| q*(s), q*′(s) | Assumption 1, Offline stage (Lift) | lift | void Lift::eval(double s, Vec7& q, Vec7* dq) const |
| X(s), X′(s) | Lemma 3, Offline stage (Lift) | lift | void KernelFrame::eval(double s, Vec7& X, Vec7* dX) const |
| β, dβ, and the level set {β = value} | eq:beta-flat, eq:dbeta | redundancy (interface, in redundancy.hpp), lift (implementation) | RedundancyOutput::evaluate(const Vec7& q, double eta, const Vec7& deta, double& beta, Vec7& dbeta) const, RedundancyOutput::levelSet(double eta, double value, Vec7& p, double& b) const |
| A(q) | eq:A | controller | filled by Controller::outputs(...) |
| μ∥, μ⋔, μ_β | eq:transverse, eq:tan-pos, eq:tan-vel, eq:beta-law | controller | struct Gains, struct TangentialReference |
| Prediction η̂ | Algorithm 1 | controller | void Controller::reset(double eta0) |
| Feedback, one cycle | eq:fb, Algorithm 1 | controller | bool Controller::step(const Config& q, const VecN& dq, double t, double dt, VecN& qdd, Status& st) |

The draft cites eq:JPsi in Algorithm 1 but no longer defines it; the formula is in the tube section below.

## Module: se3 and metric

se3 holds the group operations, which need no metric; metric holds W and the Riemannian maps Exp_I, Log_I and Γ. Build the product-type metric first: every example in the paper uses it, and its maps are closed forms.

### se3 (se3.hpp, se3.cpp)

Declare hat, vee, hat6, vee6, expSO3, logSO3, Jr, JrInv, makePose, inverse, Ad and ad. The ones you have not written yet:

```math
J_r(\omega) = I - \frac{1-\cos\theta}{\theta^2}\,\hat\omega + \frac{\theta-\sin\theta}{\theta^3}\,\hat\omega^2,\qquad \theta = \|\omega\|
```

```math
J_r(\omega)^{-1} = I + \tfrac12\hat\omega + \Bigl(\frac{1}{\theta^2} - \frac{1+\cos\theta}{2\theta\sin\theta}\Bigr)\hat\omega^2,\qquad \mathrm{ad}_\zeta = \begin{bmatrix}\hat\omega & 0\\ \hat v & \hat\omega\end{bmatrix}
```

- Small angles: below θ ≈ 1e-6 use the series J_r ≈ I − ½ω̂ + ⅙ω̂² and J_r⁻¹ ≈ I + ½ω̂ + (1/12)ω̂².
- Inverse pose: g⁻¹ = (Rᵀ, −Rᵀp) in closed form; Matrix4d::inverse() is only for tests.
- Tests, all to 1e-12 unless stated: hat(w) x ≈ w × x; logSO3(expSO3(w)) ≈ w at tiny angles and at π − 1e-3 and π − 1e-6 (to 1e-9 there); J_r(ω)δ against (exp(ω)ᵀ(exp(ω + hδ) − exp(ω − hδ)) / 2h)^∨ with h = 1e-6, to 1e-8; Ad(g)z ≈ vee6(g · hat6(z) · g⁻¹); ad(a)b ≈ vee6(hat6(a) hat6(b) − hat6(b) hat6(a)).

### metric (metric.hpp, metric.cpp)

A Metric class owns W and answers inner, norm, Exp, Log, Gamma, GammaInv and distance. Make them all const members, even though the product-type closed forms do not depend on w_R and w_p, so a general W can come later behind the same calls; Log returns false when it fails (possible only for a general W). For W = diag(w_R I, w_p I), Example 1 gives, for ζ = (ω, v):

```math
\mathrm{Exp}_I(\hat\zeta) = \begin{bmatrix}\exp(\hat\omega) & v\\ 0 & 1\end{bmatrix},\qquad \mathrm{Log}_I(g)^\vee = (\log(R)^\vee,\ p),\qquad \Gamma(\zeta) = \begin{bmatrix}J_r(\omega) & 0\\ 0 & \exp(\hat\omega)^\top\end{bmatrix}
```

- Not the screw exponential: the translation part of Exp_I is v itself. The matrix exponential of ζ̂ would give J_l(ω)v instead (Remark 1).
- Domain: Log_I needs a rotation angle below π; the injectivity radius is π√w_R.
- Tests: Log(Exp(ζ)) ≈ ζ to 1e-12; Γ(ζ)w against (Exp(ζ)⁻¹(Exp(ζ + hw) − Exp(ζ − hw)) / 2h)^∨ with h = 1e-6, to 1e-8; GammaInv · Gamma ≈ I; dist(I, Exp(ζ)) ≈ ‖ζ‖_W.

Optional: a general W. Exp_I and Γ come from eq:Exp-ode and eq:Gamma-ode. Integrate the body velocity z, the pose c = (R, p), and the six variation pairs (σ_k, δ_k) together over t ∈ [0, 1] with a fixed number of RK4 steps (32 is enough below |ζ| ≈ 1). Store the pairs as two 6 × 6 matrices S and D.

```math
\dot z = W^{-1}\mathrm{ad}_z^\top W z,\quad \dot R = R\hat\omega_z,\quad \dot p = R\,v_z,\quad \dot S = W^{-1}\bigl(\mathrm{ad}_{S}^\top W z + \mathrm{ad}_z^\top W S\bigr),\quad \dot D = S - \mathrm{ad}_z D
```

```math
z(0) = \zeta,\quad c(0) = I,\quad S(0) = I,\quad D(0) = 0,\qquad \mathrm{Exp}_I(\hat\zeta) = c(1),\quad \Gamma(\zeta) = D(1)
```

Here ad_SᵀWz means column k is ad_{σ_k}ᵀWz. RK4 moves R slightly off SO(3), so project it back once at the end (polar decomposition or SVD). Log_I becomes Newton's method, x ← x + Γ(x)⁻¹ log(Exp_I(x)⁻¹g)^∨, started from the product-type log. With a fixed step count, your Exp_I is a smooth map and your Log_I is its exact inverse.

A strong test: with W = M = diag(J, mI), geodesics are the free motions of the rigid body. Exp_I(ζ) must equal one second of your rigid-body integrator, started at the identity pose with initial velocity ζ and zero input. Expect about 56 µs per Exp + Γ in Release: fine for simulation, too slow for 1 kHz with Newton on top.

## Module: path and splines

A path answers four questions at arc length s: L, γ(s), ζ_γ(s) and ζ_γ′(s). Both example paths have constant speed in the product metric, so arc length is a rescaling of the curve parameter θ and everything is closed form.

### path (path.hpp, path.cpp)

An abstract Path with length(), pose(s), velocity(s), acceleration(s), and one subclass per curve. Give the base class virtual ~Path() = default; and mark each subclass's functions override: without the virtual destructor, deleting a subclass through a std::unique_ptr\<Path> is undefined behaviour that -Wall does not report.

Write each curve in θ ∈ [0, 2π) with body velocity ζ̃(θ) = (c⁻¹ dc/dθ)^∨ and its θ-derivative ζ̃′. With speed σ = ‖ζ̃‖_W constant:

```math
s = \sigma\theta,\qquad L = 2\pi\sigma,\qquad \zeta_\gamma(s) = \tilde\zeta(\theta)/\sigma,\qquad \zeta_\gamma'(s) = \tilde\zeta'(\theta)/\sigma^2
```

For any curve with attitude R(θ) and position p(θ), the translational parts follow from ṽ = Rᵀp′ and ṽ′ = Rᵀp″ − ω̃ × ṽ.

| Curve | Attitude R(θ), position p(θ) | ω̃(θ), ω̃′(θ) | Speed σ, length L |
| --- | --- | --- | --- |
| Coning circle (Panda) | R₀ R_z(θ) A R_z(−θ), A = exp(αê₁); p₀ + r(cos θ, sin θ, 0) | R_z(θ) w₀ and R_z(θ)(e₃ × w₀), with w₀ = Aᵀe₃ − e₃ | σ = (4w_R sin²(α/2) + w_p r²)^½; L = 0.770 for ℓ = 0.1 m |
| Tumbling circle (rigid body) | exp(θê₁) exp(θê₃); r(cos θ, sin θ, 0) | (cos θ, −sin θ, 1) and (−sin θ, −cos θ, 0) | σ = (2w_R + w_p r²)^½ = √3; L = 2√3π ≈ 10.883 |

- Pitfall: ζ_γ′ is the derivative of the body velocity, so the −ω̃ × ṽ term is not optional. If it is wrong, the normal frame drifts and stops being periodic.
- Optional, non-constant speed: compute S(θ) = ∫‖ζ̃‖_W by Gauss–Legendre on a grid, invert it by Newton, and interpolate θ(s) with a C² (quintic) spline, since ζ_γ′ needs θ″.
- Tests: ‖ζ_γ(s)‖_W ≈ 1 to 1e-12; ζ_γ against (γ(s)⁻¹(γ(s + h) − γ(s − h)) / 2h)^∨ with h = 1e-6, to 1e-8; ζ_γ′ against a central difference of ζ_γ, to 1e-6; γ(L) ≈ γ(0); the two lengths in Testing.

### spline (spline.hpp, header only)

PeriodicHermite\<T> interpolates values and exact derivatives given on n + 1 uniform nodes s_k = kL/n, where node n repeats node 0. It is a template over T (Mat65 for N, Vec7 for q* and X), so it lives entirely in the header. Declare void build(double L, std::vector\<T> y, std::vector\<T> dy) and void eval(double s, T& y, T* dy = nullptr) const; dy is the derivative of the interpolant, which is what Tube::frame, Lift::eval and KernelFrame::eval return as N′, q*′ and X′. With node spacing h = L/n and t = (s − s_k)/h on interval k:

```math
y(s) = (2t^3 - 3t^2 + 1)\,y_k + (t^3 - 2t^2 + t)\,h\,y_k' + (-2t^3 + 3t^2)\,y_{k+1} + (t^3 - t^2)\,h\,y_{k+1}'
```

```math
y'(s) = \frac{1}{h}\Bigl[(6t^2 - 6t)\,y_k + (3t^2 - 4t + 1)\,h\,y_k' + (-6t^2 + 6t)\,y_{k+1} + (3t^2 - 2t)\,h\,y_{k+1}'\Bigr]
```

- Pitfalls: wrap s into [0, L) with std::fmod plus a correction for negative values; clamp k to n − 1 when rounding lands on L; do not forget the 1/h in y′.
- Tests: interpolate (sin, cos) on a periodic grid; when n doubles, the value error should fall about 16× and the derivative error about 8×.

## Module: tube

The tube turns a pose g into the outputs (η, ξ): offline it builds the periodic normal frame N(s); online it projects g onto the path by Newton's method and evaluates J_Ψ.

class Tube keeps a const Path& and a copy of the Metric, so the path must outlive the tube. The easy way to break this is a temporary: Tube tube(ConingCircle(...), metric); compiles without a warning, builds a correct frame, and then every later call reads a destroyed object. Build paths, tubes, lifts and controllers as named objects that live for the whole run. To make the mistake a compile error, add Tube(const Path&&, const Metric&) = delete;, and the same in every class that stores a const T&.

### Offline: the normal frame (Lemma 1)

1. Initial frame: Gram–Schmidt in ⟨·,·⟩_W on ζ_γ(0), e₁, …, e₆, in that order, skipping any vector whose residual has W-norm below 1e-8. Seven vectors in R⁶ are dependent, so exactly one eᵢ drops out: e₅ on the rigid-body path (ζ_γ(0) ∝ e₁ + e₃ + e₅), e₆ on the Panda path. Check that six vectors survive; the five after ζ_γ are the columns of N₀(0).
1. Transport: integrate eq:transport by RK4 over N_grid = 2000 uniform intervals. The right-hand side needs ζ_γ and ζ_γ′ at the stage points, which the path gives exactly.
1. Holonomy: H = N₀(0)ᵀ W N₀(L). Check HᵀH ≈ I (to about 1e-10) and det H ≈ 1 before going on.
1. Logarithm of H: Φ = log H, a skew 5 × 5 matrix. RK4 leaves H orthogonal only to about 1e-11, so first replace it by the nearest rotation: with the SVD H = UΣVᵀ (Eigen::JacobiSVD\<Mat5> with full U and V), use H_r = UVᵀ. Then take Φ either from the invariant planes (the eigenvectors of (H_r + H_rᵀ)/2 have eigenvalues cos θᵢ, and H_r u = cos θ u + sin θ v defines each plane), or with Eigen's unsupported/Eigen/MatrixFunctions (Mat5 Phi = Hr.log();, Mat5 E = (t * Phi).exp();). With planes, treat a direction as fixed only if 1 − cos θ < 1e-12: the Panda's second plane has θ = 8.0e-5, so 1 − cos θ = 3.2e-9, and it must be kept. Check that Φ is skew to 1e-12 and that two planes are found for these paths. MatrixFunctions allocates memory, which is fine offline.
1. Closed frame: N(s_k) = N₀(s_k) exp(−(s_k / L) Φ), with N′(s_k) from eq:frame-ode. Store both in PeriodicHermite\<Mat65>.

Checks after construction: NᵀWN ≈ I and NᵀWζ_γ ≈ 0 to 1e-9 at random s. Closure: |N(s_n) − N(s_0)| below 1e-9 on the node values from step 5, before they go into the spline (the reference reports 4.8e-11 for the Panda path and 1.8e-13 for the rigid body). Do not measure it through frame, which wraps s = L to 0 and returns node 0 exactly. This node check is the only one that catches a wrong sign in exp(−(s/L)Φ) or a dropped plane of H; such a frame passes every other test in this section.

The frame is not unique, and the test numbers assume this one. N(s)Q, for any constant orthogonal 5 × 5 matrix Q (of either determinant), also satisfies Lemma 1, and a different Gram–Schmidt order, skip rule or sign gives such a frame. The rotation angles of H do not change (H becomes QᵀHQ), so they are a good first comparison. But ξ₀ is a coordinate vector in this frame, so the start poses Ψ(0.1, ξ₀) and most closed-loop numbers in Testing hold only for the construction in step 1. With the columns reversed, the rigid body would start 15.50° and 0.563 m off instead of 25.45° and 0.439 m. The paper says only “Gram–Schmidt”; consider stating the order there.

### Online: projection and outputs (Proposition 1)

```math
e_g(s) = \mathrm{Log}_I\bigl(\gamma(s)^{-1}g\bigr)^\vee,\qquad \varphi_g(s) = \langle e_g(s), \zeta_\gamma(s)\rangle_W,\qquad \varphi_g'(s) = \langle e_g, \zeta_\gamma'\rangle_W - \bigl\langle \Gamma(e_g)^{-1}\mathrm{Ad}_{g_\mathrm{err}^{-1}}\zeta_\gamma,\ \zeta_\gamma\bigr\rangle_W
```

- Newton: from s = guess, repeat at most 20 times: evaluate φ and φ′, update s ← s − φ/φ′, wrap s into [0, L), re-evaluate φ, and stop once |φ| < 1e-10. Always take at least one update. The iteration count is the number of updates (1 to 20), returned through int* iterations; the projections for Ȧq̇ are not counted. Start from the prediction η̂ of the previous cycle; one update suffices.
- Outputs: η is the root; ξ = N(η)ᵀ W e_g(η).
- Failure is a value, not an exception: outputs returns bool. Non-finite values, |φ′| < 1e-12, or no convergence mean the pose left the tube.
- Pitfall: φ has other roots on the circle that are not the closest point. Started from L/2, Newton converges to η = 5.614 for the rigid body's start and η = 0.485 for the Panda's, both wrong (the right answer is 0.1). For an initial condition, use globalProject: evaluate dist(γ(s_k), g) at s_k = kL/720, take the closest s_k, and run Newton from it. It is for starts only, never in the 1 kHz loop.

### Online: J_Ψ

Differentiating Ψ = γ Exp_I(Nξ) gives the body Jacobian, columns ordered (η, ξ₁, …, ξ₅):

```math
J_\Psi(\eta,\xi) = \Bigl[\ \mathrm{Ad}_{\mathrm{Exp}_I(N\xi)^{-1}}\,\zeta_\gamma(\eta) + \Gamma(N\xi)\,N'(\eta)\,\xi\quad\Big|\quad \Gamma(N\xi)\,N(\eta)\ \Bigr]
```

### Offline: tube-radius estimate

Follow the paper's test: ε passes if no point of a fine sampling of the path is closer than ε to any pose Ψ(s, εd). Use 32 values of s, the ten directions ±eᵢ plus 6 random unit directions, 1000 path samples plus a few points just beside each foot point, and a relative tolerance of 1e-7. Double ε from 0.01L until it fails, then bisect 14 times. The result depends on the sampled directions; expected values are in Testing.

### Tests

- Ψ⁻¹∘Ψ round trip: (η, ξ) recovered to 1e-9 for |ξ| up to 0.06 on the Panda path.
- |ξ| equals the minimum distance over 20 000 path samples, to 2e-6 (Proposition 1(c)).
- J_Ψ against a central difference of Ψ along a random direction (h = 1e-6), to 1e-6.
- Newton converges from guesses offset by ±0.05; φ′(η) ≈ −1 for a pose on the path.
- globalProject(Ψ(0.1, ξ₀)) returns 0.1 for both example starts, where Newton from L/2 does not.

## Module: systems

A system is anything with a configuration q, a pose output Π(q) with body Jacobian J_b(q), and dynamics of the form eq:EL. The controller only needs the kinematic part, so that is the interface; the dynamics (M, Cζ and the integrator for the rigid body) are used only by the simulators and the robot driver.

### The system interface

The two systems differ in kind: the rigid body's configuration is a pose (n = 6), the Panda's a vector of joint angles (n = 7). Both n and the configuration type are known at compile time, so the controller is a template over a System class that provides:

| Member | RigidBody | PandaArm |
| --- | --- | --- |
| static constexpr int kDof | 6 | 7 |
| using Config | Pose | Vec7 |
| static Pose output(const Config& q, Mat<6, kDof>* Jb) | Π = id, J_b = I | flange forward kinematics, body Jacobian |
| static Config retract(const Config& q, const Vec\<kDof>& dq, double h) | (R expSO3(hω), p + hRv) | q + h dq |

retract gives a curve through q with velocity dq, which the controller uses for the central difference of A along the motion. Any smooth curve with c(0) = q and c′(0) = dq works, because (A(c(h)) − A(c(−h)))/2h = DA(q)[dq] + O(h²) for any such curve: c″ changes only the size of the O(h²) error, not the limit. For the rigid body this curve is g·Exp_I(hζ) of the product metric, built from functions you already have; the screw exponential would also work but is not needed. A good way to learn templates: write the controller for the rigid body only first, then generalize (pitfalls in the controller section).

### Free rigid body (rigid_body.hpp, .cpp)

```math
M = \mathrm{diag}(J, mI),\qquad C(g,\zeta)\zeta = -\bigl(J\omega\times\omega,\ m\,v\times\omega\bigr),\qquad G = 0,\qquad \dot\zeta = M^{-1}(u - C\zeta)
```

- Integrator: one step of length h with u held constant. Use RK4 for ζ, and for the pose use exponential coordinates around the step's starting attitude R_k: R = R_k exp(φ̂), φ̇ = J_r(φ)⁻¹ω, ṗ = Rv. At the end set R = R_k exp(φ̂), so R never leaves SO(3).
- Sample-and-hold: the controller runs every Δt = 2 ms; u is held over the period while the plant takes 4 substeps of 0.5 ms.
- Test, energy: J = diag(1.0, 1.2, 0.8), m = 2, ζ₀ = (1, −0.7, 0.5, 0.4, 0.1, −0.3), u = 0, h = 0.5 ms for 10 s. The relative drift of E = ½ζᵀMζ must stay below 1e-10 (the reference reaches 6.5e-15), and |RᵀR − I| below 1e-12.
- Test, a closed-form motion: the energy test cannot see a sign error in Cζ, because both terms are orthogonal to the velocity. Use J = diag(1, 1, 2), m = 2, ω₀ = (1, 0, 1), any v₀ (for example (0.3, −0.2, 0.5)), u = 0, and h = 0.5 ms for 10 s; at h = 2 ms, RK4's own error is already about 1e-12. The body rate must follow ω(t) = (cos t, sin t, 1) to 1e-12, and the centre of mass must move in a straight line: Rv stays equal to R₀v₀, and p(t) = p₀ + tR₀v₀, to 1e-10. A flipped sign gives (cos t, −sin t, 1) and a curved track.

### Franka Panda (panda.hpp, .cpp)

- Forward kinematics: modified (Craig) DH, T = ∏ᵢ Tᵢ(qᵢ) · Trans_z(0.107) for the flange, with the parameters in the appendix:

```math
T_i(\theta) = \begin{bmatrix}\cos\theta & -\sin\theta & 0 & a_{i-1}\\ \sin\theta\cos\alpha_{i-1} & \cos\theta\cos\alpha_{i-1} & -\sin\alpha_{i-1} & -d_i\sin\alpha_{i-1}\\ \sin\theta\sin\alpha_{i-1} & \cos\theta\sin\alpha_{i-1} & \cos\alpha_{i-1} & d_i\cos\alpha_{i-1}\\ 0&0&0&1\end{bmatrix}
```

- Body Jacobian: joint i's frame has z-axis zᵢ and origin pᵢ in the base frame. Its spatial column, in (ω, v) order, is (zᵢ, pᵢ × zᵢ); then J_b = Ad_{g⁻¹} J_s with g the flange pose.
- Limits: joint position, velocity, acceleration and torque limits from the Franka docs (appendix), plus a jointMargin(q) helper: the distance to the nearest position limit.
- Dynamics: not needed for the simulations. The paper's Panda results are at the joint-acceleration level, where eq:fb cancels the dynamics exactly. On the robot, franka::Model provides M, C and G, and the robot adds G itself.
- Tests: FK at q = 0 and at the rest posture against the values in Testing; column j of J_b against (g⁻¹(FK(q + heⱼ) − FK(q − heⱼ)) / 2h)^∨ with h = 1e-6, to 1e-8, at random q.

## Module: lift and redundancy output

For the Panda (n = 7, r = 1), the lift q*(s) picks one configuration at each point of the path, and β measures how far q is from it along the self-motion direction X(s). For the rigid body, q* = γ and this module is not used.

### Offline: the lift q* (Assumption 1)

The paper defines q*(s) as the configuration nearest the rest posture q_r on the self-motion manifold over γ(s): it minimizes ‖q − q_r‖ subject to Π(q) = γ(s), “by a Newton-type method started from the solution at the previous grid point”. Two solvers work. The projected step needs only J_b; the square system needs x(q), the unit kernel vector of J_b(q) (Kernel frame, below), which you need anyway for the minimum check and for X(s).

- Projected Newton step (simplest): with J = J_b(q), J⁺ = Jᵀ(JJᵀ)⁻¹ and e the pose error from FK(q) to γ(s) in the flange frame,

```math
\Delta q = J^{+}e + \bigl(I - J^{+}J\bigr)(q_\mathrm r - q)
```

- Stopping rule: |e| < 1e-12 and |(I − J⁺J)(q − q_r)| < 1e-12. The second condition is the Lagrange condition: q − q_r lies in the range of J_bᵀ. It converges linearly: 9 updates per node, 13 at most.
- Square 7 × 7 system (faster): for r = 1 the Lagrange condition is the scalar h(q) = x(q)ᵀ(q − q_r) = 0. Solve [J_b; ∇hᵀ] Δq = (e, −h), 2–3 updates per node. This is what the reference library uses, with |Δq| capped at 0.2 rad; both solvers reach the same q* to about 1e-12, so either reproduces the Testing values. Its gradient is ∇h = Dx(q)ᵀ(q − q_r) + x(q), with Dx the 7 × 7 Jacobian of x, by central differences with step 1e-6. Do not drop the Dxᵀ term: |q* − q_r| is 0.4–1.0 rad on this lift, and with xᵀ alone as the last row you get the projected step back, not Newton.
- Continuation: nodes s_k = kΔs, Δs = L/2000. Solve at the first node s = 0 from q_r. FK(q_r) is 0.255 m and 0.25 rad from γ(0), so the first step is large (about 0.97 rad). The projected step needs no step limit there (11 updates); the square system is safer with |Δq| capped at 0.2 rad. Start every later node from the previous solution; those steps are about 0.002 rad.
- Checks: closure |q*(L) − q*(0)| < 1e-6 (otherwise the branch ended on a different self-motion and Assumption 1 fails); σ_min(J_b(q*)) > 0 at every node; every joint inside its limits with margin.
- Is it a minimum, not a maximum? At each node, find q₊ and q₋ on the same self-motion: solve FK(q) = γ(s_k) with the linear condition xᵀq = xᵀq*(s_k) ± 0.05 (x = x(q*(s_k))) by the 7 × 7 Newton of Start configurations, from q*. Then check f(q₊) − 2f(q*) + f(q₋) > 0 for f = ‖q − q_r‖². Do not use q* ± 0.05x unprojected: f is quadratic, so that difference is always 2·0.05² and tests nothing.
- Derivative q*′: the spline needs it at the nodes. Fourth-order periodic central differences of the node values, (q_{k−2} − 8q_{k−1} + 8q_{k+1} − q_{k+2}) / 12Δs, are accurate enough. The exact value solves [J_b(q*); ∇h(q*)ᵀ] q*′ = (ζ_γ, 0), with the full ∇h. The fourth-order differences agree with it to about 1.4e-9 (largest over all nodes); a second-order difference (q_{k+1} − q_{k−1}) / 2Δs is off by up to 1.9e-5, so an error near 1e-5 means the stencil is wrong.

### Offline: the kernel frame X (Lemma 3)

For r = 1 no integration is needed (Implementation, Lift): X(s) is the normalized generalized cross product of the six rows of J_b(q*(s)),

```math
n_j = (-1)^j\det\bigl(J_b \text{ without column } j\bigr),\quad j = 0,\dots,6,\qquad x = n/\|n\|,\qquad J_b\,n = 0
```

It is smooth in q and depends only on q*(s), so it is periodic with no sign choice to make. Get X′ at the nodes the same way as q*′, and store q*, q*′, X and X′ in PeriodicHermite\<Vec7>. For r > 1 you would transport X like N (Lemma 3) and cancel its holonomy; not needed here.

### Online: β and dβ (eq:beta-flat, eq:dbeta)

```math
\beta(q) = X(\eta)^\top\bigl(q - q^\star(\eta)\bigr),\qquad \mathrm d\beta(q) = X(\eta)^\top + \Bigl(X'(\eta)^\top\bigl(q - q^\star(\eta)\bigr) - X(\eta)^\top q^{\star\prime}(\eta)\Bigr)\,\mathrm d\eta(q)
```

Here dη(q) is the first row of J_Ψ⁻¹J_b, which the controller has already computed. Put β behind an abstract class RedundancyOutput with a virtual destructor and two methods, so a second output can be swapped in:

- evaluate(q, eta, deta, beta, dbeta): the formulas above.
- levelSet(eta, value, p, b): among configurations with η(Π(q)) = η, the set {β = value} is the hyperplane pᵀq = b. Here p = X(η) and b = X(η)ᵀq*(η) + value.

Declare the interface in its own header, redundancy.hpp, in the outputs layer: just the abstract class, which needs only Vec7 from types.hpp and has no .cpp. Write it at milestone 6, because Controller\<RigidBody> already holds a const RedundancyOutput*. controller.hpp and lift.hpp both include it, and the implementations live here in lift. Declaring it in controller.hpp instead would make lift include the layer above it. The joint coupling makes a useful second output on a lift with q₃ = q₁: β = q₃ − q₁, dβ = (e₃ − e₁)ᵀ, p = e₃ − e₁, b = value.

### Start configurations

To start at pose Ψ(s₀, ξ₀) with β = β₀ (s₀ = 0.1 in the examples), get p and b from levelSet(s₀, β₀, p, b), then solve FK(q) = Ψ(s₀, ξ₀) together with pᵀq = b by Newton on the 7 × 7 system [J_b; pᵀ], starting from q*(s₀). This is correct because η = s₀ for that pose.

### Tests

- FK(q*(s)) ≈ γ(s) to 1e-8 at random s from the spline. Closure |q*(s_n) − q*(s_0)| < 1e-6 (reference 5.9e-13), measured on the continuation's node values before they go into the spline: the spline wraps s = L to 0 and returns q*(0) exactly, so it cannot test closure.
- Off the nodes, at s = s_k + Δs/4, Δs/2, 3Δs/4: |X(s) − x(q*(s))|, |J_b(q*(s)) X(s)| and ||X(s)| − 1| all below 1e-12; X′ against a central difference of the X spline, to 1e-9. No sign flips: X(s_k)·X(s_{k+1}) > 0.9999 at every node.
- β(q*(s)) ≈ 0 to 1e-12, and dβ(q*(s)) X(s) = +1 to 1e-12. It is exactly +1, because dη X = 0 on the lift; allowing −1 would hide a sign error.
- The start configuration has the requested pose and β to 1e-10.
- Compare q*(0), the smallest joint margin and the closure gap with the numbers in Testing.

## Module: controller

The controller returns the generalized acceleration q̈* = A(q)⁻¹(μ − Ȧq̇) that makes ȳ̈ = μ; the caller turns it into an input u = M q̈* + C q̇ + G. That split keeps the controller free of dynamics, so one class serves the rigid body, the Panda simulation and the robot.

### Interface

template \<class System> class Controller holds a const Tube&, the Gains, a TangentialReference, the prediction η̂, and a const RedundancyOutput* that must be non-null when n = 7 (throw in the constructor otherwise).

- void reset(double eta0) sets η̂ = eta0 wrapped into [0, L). Call it before the first step of every run; also set η̂ = 0 in the constructor so it is never uninitialized.
- bool outputs(const Config& q, double guess, OutputState& out) const fills η, ξ, β, A and the Newton iteration count.
- bool step(const Config& q, const VecN& dq, double t, double dt, VecN& qdd, Status& st) runs one cycle.
- OutputState and Status are structs nested in the class, so their sizes follow from System. Status holds t, η, η̇, the reference (ν_d or η_d), the tangential error (e_v or e_η), ξ, ξ̇, β, β̇, σ_min(A), the iteration count and an ok flag.
- Gains: K_p = λ²I, K_d = 2λI with λ = 4 rad/s (poles at −λ; K_p = 16, K_d = 8); k_t = 5 for velocity references; k_p = 9, k_d = 6 for the rigid body's pose reference; K_p,β = 16, K_d,β = 8.
- TangentialReference: either a velocity ν_d(t, η) with ∂_tν_d and ∂_ην_d, or a position η_d(t) with η̇_d and η̈_d. std::function members work; calling one does not allocate.

Template pitfalls (C++17), all with confusing error messages:

- Inside the class, write using Config = typename System::Config; and static constexpr int n = System::kDof; using VecN = Vec\<n>;. Without typename, GCC says “need 'typename' before 'System::Config'”.
- Fixed-size member templates on objects whose type depends on n need the template keyword: A.template topRows<6>(). Without it GCC says only “expected primary-expression”. The runtime forms A.topRows(6) and A.block(0, 0, 6, n) need no keyword.
- Put every n = 7 step (β, dβ, the seventh row of A, μ_β) under if constexpr (n == 7). A plain if still compiles that code for Controller\<RigidBody>, where Config is a 4 × 4 Pose, and Eigen stops with YOU_MIXED_MATRICES_OF_DIFFERENT_SIZES.

### outputs (eq:A)

1. g = Π(q) and J_b from System::output.
1. (η, ξ) from Tube::outputs, starting at the guess.
1. Factor J_Ψ(η, ξ) once (partialPivLu()) and solve for the 6 × n block J_Ψ⁻¹J_b; never form the inverse. Return false if the block is not allFinite().
1. If n = 7: β and dβ from the redundancy output, with dη = row 0 of that block. A = [J_Ψ⁻¹J_b; dβ].

### step (Algorithm 1)

1. outputs(q, η̂); ȳ̇ = A q̇, which holds η̇, ξ̇ and β̇.
1. Ȧq̇ ≈ (A(q₊) − A(q₋)) q̇ / 2h with q± = retract(q, dq, ±h), h = 1e-5, both projections started from the current η. When q̇ = 0 you may skip the two projections, but set Ȧq̇ = VecN::Zero(): Eigen leaves it as garbage otherwise.
1. μ∥ from eq:tan-vel (e_v = η̇ − ν_d) or eq:tan-pos, with e_η wrapped into (−L/2, L/2].
1. μ⋔ = −K_pξ − K_dξ̇; if n = 7, μ_β = −K_p,ββ − K_d,ββ̇.
1. q̈* = A⁻¹(μ − Ȧq̇) with partialPivLu().solve(...). That factorization does not check rank: a nearly singular A gives huge but finite values. So compute σ_min(A) with Eigen::JacobiSVD\<Mat\<n, n>>(A).singularValues().minCoeff() (fixed size, no allocation) and return false if it is below 1e-9 or q̈* is not finite.
1. η̂ ← (η + η̇Δt) mod L for the next cycle; fill Status, using the σ_min(A) from step 5.

### Rules for the 1 kHz loop

- No heap allocation in step: fixed-size Eigen types only, no MatrixXd, no std::vector growth, no strings, no new. No console or file output either: log into a buffer allocated before the run and write it afterwards.
- No exceptions in step: failure is false. A simulation stops and reports the time; the robot driver stops the arm with the safety stop (a damped hold) and brings it to rest.
- For hardware only: multiply μ (not Ȧq̇) by the quintic smoothstep r = 10x³ − 15x⁴ + 6x⁵, x = t/T clamped to [0, 1], T = 1–2 s. Its first two derivatives vanish at both ends, so the first command at rest is exactly zero and the torque starts continuously. Use the same ramp in the preflight.

### Tests

- A is the Jacobian of ȳ: compare (ȳ(retract(q, w, h)) − ȳ(retract(q, w, −h))) / 2h with A w for random w, h = 1e-6, to 1e-6, for both systems. For the rigid body w is a body twist; never add a vector to a Pose. Wrap the η difference into (−L/2, L/2] before dividing. Run it off the lift: at the rigid body's start pose, and at the start configurations of both Panda runs. On the lift q − q*(η) = 0, so the X′ term of dβ drops out and a sign error in it goes unseen.
- On the lift: ξ ≈ 0, β ≈ 0 and σ_min(A) > 0 (Proposition 2).
- Closed loops reproduce the designed responses (numbers in Testing).
- No allocation in step: compile with target_compile_definitions(se3pf PUBLIC EIGEN_RUNTIME_NO_MALLOC). Build everything, call Eigen::internal::set_is_malloc_allowed(false), run 1000 Panda steps with q̇ ≠ 0, then allow it again. Any Eigen allocation trips an assert, which exists only in Debug builds. That does not see std::vector or std::function: for those, replace the global operator new in the test file with one that counts calls, read the counter just before and just after the 1000 steps, and REQUIRE the two to be equal after the loop. Keep CAPTURE, INFO and other checks outside that window, because Catch2 allocates for them.
- Step time: below 100 µs for the Panda with the product metric, in the Release build (the paper reports 27 µs). A Debug build takes about 1–3 ms, longer than the 1 ms control period. Tag the test [.timing] so the default Debug run skips it.

## Simulation drivers and the hardware layer

Simulators are small programs in apps/ that wire the library together and write CSV logs; the robot driver is the same loop with libfranka in place of the integrator, and comes last. If step returns false in a simulation, stop and report the time.

### Rigid-body simulation (apps/rigid_body_sim.cpp)

1. Build the tumbling-circle path (r = 1 m, W = I), its tube, and the controller with λ = 4.
1. Start at rest: g₀ = Ψ(0.1, ξ₀) with ξ₀ = (0.4, −0.3, 0.2, −0.1, 0.3), ζ₀ = 0. Call reset(0.1): a pose Ψ(s₀, ξ) has η = s₀.
1. Every Δt = 2 ms: step gives q̈*; u = M q̈* + Cζ is held for the period while the integrator takes 4 substeps.
1. Run twice for 20 s: the pose reference η_d(t) = L/4 + (L/5)(1 − cos(πt/10)) with k_p = 9, k_d = 6, and the velocity reference ν_d(t) = 0.5 + 0.2 sin(0.5t) with k_t = 5.
1. Log t, η, η̇, reference, tangential error, ξ₁…ξ₅, |ξ|, the designed |ξ₀|(1 + 4t)e^{−4t} and u. A short Python script plots the paper's three panels: transverse output, pose tracking, velocity tracking.

### Panda simulation (apps/panda_sim.cpp)

The paper's Panda results are at the joint-acceleration level: with an exact model, eq:fb makes q̈ = q̈*, so no dynamic model is needed.

1. Build the coning-circle path (ℓ = 0.1 m), tube, lift, kernel frame and controller.
1. Start at rest at the configuration with pose Ψ(0.1, ξ₀), ξ₀ = (0.02, −0.015, 0.01, −0.005, 0.015), and β = 0; call reset(0.1).
1. Every Δt = 1 ms: q ← q + Δt q̇ + ½Δt²q̈*, q̇ ← q̇ + Δt q̈*. This is exact for an acceleration held over the period. Velocity reference ν_d = 0.08, 20 s.
1. Log the same columns plus q, q̇ and σ_min(A); report the smallest joint margin, the largest |q̇ᵢ| as a fraction of its limit, and the smallest σ_min(A).
1. Self-motion run (checks the β channel; not in the paper): start on the path, with pose γ(0.1) (ξ₀ = 0) and β₀ = 0.2. Set the initial velocity by solving A(q₀) q̇₀ = (0.08, 0, …, 0), so that η̇ = 0.08 and ξ̇ = β̇ = 0 at t = 0, as the designed response β₀(1 + 4t)e^{−4t} assumes. Step while kΔt < L/0.08: 9627 steps, one loop. Report max |ξ|, max |β − β₀(1 + 4t)e^{−4t}|, and maxᵢ |qᵢ − q*ᵢ(η)| at the end.

### Hardware (last milestone)

> Implemented on 2026-10-15 (HW-1, HW-2): `se3pf/safety.hpp`, `se3pf/experiments.hpp`, `se3pf/preflight.hpp`,
> `apps/robot/robot_driver.cpp`, `mock/franka/`; the lab runbook is docs/HARDWARE.md. Where this section and the
> code differ, the code and docs/HARDWARE.md are right.

Only this executable links libfranka 0.9.2, the newest version that supports the Panda (FER); 0.10 and later support only the FR3. The library itself never includes a libfranka header.

1. Check kinematics: compare FK and J_b with model.pose(franka::Frame::kFlange, state) and model.bodyJacobian(franka::Frame::kFlange, state). Both are column-major arrays; the Jacobian is in the flange frame like yours, but its rows are (v, ω), so swap the two 3-row blocks. Do not use kEndEffector, which includes the gripper offset set in Desk, or zeroJacobian, which is in the base frame. Expect agreement to about 1e-5.
1. Preflight: simulate exactly the motion the robot will run (same start, reference, gains and μ ramp), evaluate torques with franka::Model, and refuse to move if any margin below is exceeded. The ramp is what lets the paper's start pass: without it the first command asks for 49 % of a joint's acceleration limit; with the quintic ramp over 1 s the peak is 11.5 % (22 % of velocity), over 2 s it is 6.8 % (19 %).
1. Approach: move to the start configuration with a slow minimum-jerk joint motion (franka::JointPositions).
1. Torque loop: in robot.control(...) at 1 kHz: step, then τ = M q̈* + C q̇ (the robot adds G), check the margins, send franka::Torques. On any violation, switch to the safety stop until the arm is at rest. (As implemented, after the safety review: a damped hold at the configuration of the violation, τ = −K_s (q − q_stop) − D q̇, `SafetyStop` in se3pf/experiments.hpp; joint damping alone, τ = −D q̇, cannot hold a weight that G does not know.)

Margins used by the earlier implementation: at least 0.10 rad from every joint limit; at most 50 % of joint velocity and 30 % of commanded acceleration; flange speed at most 0.5 m/s; |ξ| ≤ 0.04; |β| ≤ 0.05 rad; σ_min(A) ≥ 1e-3. The torque margin applies to the torque the joint delivers, |τᵢ + Gᵢ(q)| ≤ 0.6 τ_max,i; the torque-rate margin applies to the commanded τ, at most 60 % of 1000 N·m/s. To test the driver without a robot, write a small mock of the few libfranka classes you use; mock/franka in the earlier repo is an example.

## Testing and acceptance numbers

Every derivative is checked against a finite difference, every construction against its defining property, and every simulation against a number from the paper or the reference implementation. The numbers below were computed with the reference library for the paper's definitions, including its rest-posture lift, and each was recomputed independently; the paper's own figures are in brackets.

### Three kinds of test

- Finite differences: J_r, Γ, ζ_γ, ζ_γ′, N′, J_b, J_Ψ, A, q*′ and X′, each against a central difference with h = 1e-6, to 1e-6 or better. Check A at configurations off the lift, where the X′ term of dβ is active.
- Properties: round trips (log∘exp, Log∘Exp, Ψ⁻¹∘Ψ), orthonormality and periodicity of N and X, unit speed, closure of the lift.
- Closed loops: ξ against the designed response |ξ₀|(1 + 4t)e^{−4t}. The deviation is a sample-and-hold effect of order Δt: in the rigid body's velocity-reference run it is 9.5e-4 at 2 ms, 4.75e-4 at 1 ms and 2.37e-4 at 0.5 ms. If halving Δt does not halve your deviation, look for a bug.
- Writing the checks: tolerances here are absolute. For vectors and matrices use REQUIRE((a - b).cwiseAbs().maxCoeff() < tol), or REQUIRE(a.isZero(tol)) against zero. For scalars use REQUIRE_THAT(x, Catch::Matchers::WithinAbs(expected, tol)) from \<catch2/matchers/catch_matchers_floating_point.hpp>. Put CAPTURE(err) before a check so a failure prints the actual error.

### Unit-test anchors

| Quantity | Expected value | Tolerance and notes |
| --- | --- | --- |
| FK(q = 0) | flange at (0.088, 0, 0.926) m, R = diag(1, −1, −1) | 1e-12 |
| FK(q_r) | flange at (0.306891, 0, 0.590282) m, R rows (c, −c, 0), (−c, −c, 0), (0, 0, −1), c = √½ | 1e-6 |
| Panda path, ℓ = 0.1 m | L = 0.77008765, σ = 0.12256326, p(γ(0)) = (0.52, 0, 0.45) m [0.770] | 1e-8 |
| Panda path, ℓ = 1 m | L = 1.73869664 | 1e-8 |
| Rigid-body path | L = 2√3π = 10.88279619, computed in the test as 2 * std::sqrt(3.0) * M_PI; ζ_γ(0) = (1, 0, 1, 0, 1, 0)/√3 | 1e-8 |
| Holonomy angles of H | Panda ℓ = 0.1: 0.161182 and 8.0e-5; Panda ℓ = 1: 0.707417 and 3.09e-4; rigid body: 2.166040 and 0.492067 | 1e-5; independent of the initial frame and the grid |
| Tube radius estimate | Panda ℓ = 0.1: 0.1250 [0.125]; Panda ℓ = 1: 0.307 [0.279]; rigid body: 1.93 [1.89] | ±10 %: depends on the sampled directions (1000 directions give 0.283 and 1.87) |
| logSO3 near π | round trip to 1e-9 at π − 1e-6 with the symmetric-part formula | the reference's near-π branch is only first-order (error ≈ 1.76(π − θ)), so do not compare against it there |

### The paper's lift (Panda, ℓ = 0.1 m)

| Quantity | Expected value |
| --- | --- |
| Closure \|q*(L) − q*(0)\| | 5.9e-13 |
| q*(0) | (0.01457, −0.13186, −0.02614, −2.21552, 0.22458, 2.24779, 0.64710) rad |
| Smallest joint margin | 0.3592 rad, joint 4 at its lower limit, near s = 0.367 |
| Smallest σ_min(J_b) | 0.2070 |
| Largest \|q*′\| | 5.28 rad per unit arc length |
| Newton updates to 1e-12 | projected step: 9 typical, 13 at most, 11 from q_r at s = 0 (14 with a 0.2 rad cap); square 7 × 7 system: 2–3, 8 at s = 0 with a 0.2 rad cap (5 without) |
| Is it a minimum? | yes: with q± on the self-motion at xᵀq = xᵀq* ± 0.05, the second difference of \|q − q_r\|² is positive at all 2000 nodes, at least 4.1e-3 (the unprojected value would be 5e-3 everywhere) |
| Kernel frame | no sign flips: min over k of X(s_k)·X(s_{k+1}) = 0.99999966; off-node errors of the X spline about 3e-13 (see the lift tests) |
| For comparison: coupled lift q₃ = q₁ | closure 5.3e-13, margin 0.3604 rad, σ_min(J_b) 0.2078; differs from the rest-posture lift by up to 0.134 rad |

### Closed loops

These values assume N₀(0) is built exactly as in step 1 of the tube module; with another valid frame the start errors, max |u| and the exact deviations change, although the rigid body's deviation stayed below 1e-3 in every frame tried.

| Run | Expected values |
| --- | --- |
| Rigid body, start | \|ξ₀\| = 0.6245; 25.45° and 0.4389 m from γ(0.1) [25.5°, 0.44 m] |
| Rigid body, pose reference, 20 s | max \|\|ξ\| − design\| = 9.69e-4 [within 1e-3]; final \|ξ\| 1.8e-6; final e_η −1.7e-6; max \|u\|∞ 29.1 (at t = 0) |
| Rigid body, velocity reference, 20 s | max deviation 9.51e-4 [within 1e-3]; final \|ξ\| 2.8e-6; final e_v −3.2e-6; max \|u\|∞ 14.5 (at t = 0) |
| Panda, start | \|ξ₀\| = 0.03122; 15.07° and 16.82 mm from γ(0.1) [15.1°, 16.8 mm] |
| Panda, ν_d = 0.08, 20 s | max deviation 2.43e-5 [2.4e-5]; max \|β\| 1.5e-5 (sampling, O(Δt)); final \|ξ\| 7.1e-7; final η̇ 0.0800025 |
| Panda, along that run | smallest joint margin 0.3592 rad; smallest σ_min(A) 0.0474; largest \|q̇ᵢ\| 26 % of its limit; at most one Newton update per cycle, counted as in the tube section (the reference reports exactly 1) |
| Panda, self-motion (step 5 of the Panda simulation) | max \|ξ\| 1.2e-6; max \|β − β₀(1 + 4t)e^{−4t}\| 1.51e-4 (accept below 2e-4); at the end maxᵢ \|qᵢ − q*ᵢ(η)\| 9.9e-6, 2-norm 1.3e-5 (accept below 2e-5) |
| Panda, step time | 24–28 µs per cycle at -O2 on a Ryzen 5 5600G [27 µs, 2.1 GHz Xeon], measured with the reference's own matrix code; expect the same order of magnitude with Eigen in Release |

## Milestones

The first end-to-end result is the rigid-body figure (milestone 6), before any Panda code. Each milestone ends when its check passes; the numbers come from the Testing section.

- 1. Skeleton: CMake, Eigen, Catch2 and passing tests; headers under include/se3pf/, aliases in types.hpp, everything in namespace se3pf. Check: ./build/tests runs and passes.
- 2. se3: hat, vee, hat6, vee6, expSO3, logSO3, Jr, JrInv, poses, Ad, ad. Check: round trips at tiny angles and near π; Ad and ad against their definitions.
- 3. Product metric: inner, Exp, Log, Γ, Γ⁻¹, distance. Check: Γ against finite differences of Exp to 1e-8.
- 4. Spline and paths: PeriodicHermite\<T>, the coning and tumbling circles. Check: unit speed; derivatives against finite differences; L = 0.77008765 and 2√3π.
- 5. Tube: normal frame, Ψ, projection, outputs, globalProject, J_Ψ. Check: frame properties and node closure; Proposition 1 tests; J_Ψ against finite differences; holonomy angles match Testing.
- 6. Rigid body end to end: dynamics and integrator, redundancy.hpp (the interface only), Controller\<RigidBody>, rigid_body_sim, plot script. Check: the two rigid-body tests; start error 25.45° and 0.4389 m; |ξ| within 1e-3 of the designed response; the paper's three panels.
- 7. Tube radius: the bisection estimate. Check: within 10 % of 1.93 for the rigid body and 0.125 for the Panda path.
- 8. Panda kinematics: DH forward kinematics, body Jacobian, limits. Check: FK at q = 0 and at q_r; J_b against finite differences.
- 9. Lift, kernel frame, β: continuation, checks, splines, redundancy output, start configurations. Check: the lift table in Testing; β = 0 and dβ·X = +1 on the lift.
- 10. Panda end to end: Controller\<PandaArm>, panda_sim with both runs. Check: the Panda rows in Testing; A against finite differences off the lift; no allocation in step; step time below 100 µs in Release.
- 11. Optional, general W: geodesic ODEs, Newton Log, numerical arc length. Check: the W = M rigid-body test; closed forms recovered for a product W.
- 12. Hardware, last: libfranka driver against a mock first, then preflight, approach, torque loop and safety stop on the robot. Check: kinematics agree with franka::Model at the flange; the preflight passes before any motion.

## Open decisions

These are yours to make; each has a default that works, given in brackets. Tick them off as you decide.

- Headers vs. sources: split non-template code into .hpp and .cpp to learn the compile-and-link model, and keep templates (PeriodicHermite, Controller) in headers. [split]
- Metric scope: a product-type-only Metric first, with the general-W path added behind the same interface later, or both from the start. [product first]
- Lift solver: the projected Newton step, or the square 7 × 7 system with h(q) = x(q)ᵀ(q − q_r), which the reference uses (the paper says only “a Newton-type method”). [projected Newton, then compare]
- Node derivatives of q and X:* fourth-order central differences of the node values, or the exact derivative from the implicit function theorem. [differences]
- Polymorphism: virtual Path and RedundancyOutput (chosen at run time), template Controller\<System> (fixed at compile time). [as stated]
- Failure reporting: bool plus out-parameters, or std::optional\<OutputState>, a good C++17 feature to learn. [either]
- Which lift the paper reports: the rest-posture lift of the Implementation section, or also the joint coupling q₃ = q₁ as a second redundancy output for comparison. [rest posture]

## Appendix: example parameters

All numbers the two examples need, in one place. The Panda's kinematics and limits are the values of the earlier implementation's panda.hpp, which takes them from Franka's control parameters for the Panda (FER).

### Rigid body

| Parameter | Value |
| --- | --- |
| Path | tumbling circle, r = 1 m, centred at the origin, R₀ = I |
| Metric | W = I (w_R = w_p = 1, ℓ = 1 m) |
| Inertia, mass | J = diag(1.0, 1.2, 0.8) kg·m², m = 2 kg |
| Start | at rest at Ψ(0.1, ξ₀), ξ₀ = (0.4, −0.3, 0.2, −0.1, 0.3) |
| Gains | λ = 4 rad/s; pose reference k_p = 9, k_d = 6; velocity reference k_t = 5 |
| References | η_d(t) = L/4 + (L/5)(1 − cos(πt/10)); ν_d(t) = 0.5 + 0.2 sin(0.5t) |
| Timing | control period 2 ms with zero-order hold, 4 plant substeps, 20 s |

### Panda

| Parameter | Value |
| --- | --- |
| Path | coning circle, p₀ = (0.40, 0, 0.45) m, r = 0.12 m, α = 0.25 rad, R₀ = flange attitude at q_r (approach axis pointing down) |
| Metric | w_R = 0.01, w_p = 1 (ℓ = 0.1 m); also w_R = 1 for the role-of-the-metric example |
| Rest posture | q_r = (0, −π/4, 0, −3π/4, 0, π/2, π/4) rad |
| Start | at rest at Ψ(0.1, ξ₀), ξ₀ = (0.02, −0.015, 0.01, −0.005, 0.015), β = 0 |
| Gains | λ = 4 rad/s; k_t = 5; K_p,β = 16, K_d,β = 8 |
| Reference, timing | ν_d = 0.08 (about 9.6 s per loop); control period 1 ms; 20 s |

### Panda kinematics: modified DH (Craig)

| Joint i | aᵢ₋₁ (m) | dᵢ (m) | αᵢ₋₁ (rad) |
| --- | --- | --- | --- |
| 1 | 0 | 0.333 | 0 |
| 2 | 0 | 0 | −π/2 |
| 3 | 0 | 0.316 | π/2 |
| 4 | 0.0825 | 0 | π/2 |
| 5 | −0.0825 | 0.384 | −π/2 |
| 6 | 0 | 0 | π/2 |
| 7 | 0.088 | 0 | π/2 |
| flange | 0 | 0.107 | 0 |

### Panda limits

| Joint | q_min (rad) | q_max (rad) | q̇_max (rad/s) | q̈_max (rad/s²) | jerk max (rad/s³) | τ_max (N·m) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | −2.8973 | 2.8973 | 2.1750 | 15 | 7500 | 87 |
| 2 | −1.7628 | 1.7628 | 2.1750 | 7.5 | 3750 | 87 |
| 3 | −2.8973 | 2.8973 | 2.1750 | 10 | 5000 | 87 |
| 4 | −3.0718 | −0.0698 | 2.1750 | 12.5 | 6250 | 87 |
| 5 | −2.8973 | 2.8973 | 2.6100 | 15 | 7500 | 12 |
| 6 | −0.0175 | 3.7525 | 2.6100 | 20 | 10000 | 12 |
| 7 | −2.8973 | 2.8973 | 2.6100 | 20 | 10000 | 12 |

Torque rate is limited to 1000 N·m/s on every joint, and flange speed to 1.7 m/s.
