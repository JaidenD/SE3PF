# Adversarial review of the harmonic-lift hypothesis, and what beats it

Role: adversarial mathematician. Scope: (1) attack the working hypothesis ("the canonical closed lift is the harmonic
section of the self-motion bundle; Lemma 1, Lemma 3 and the lift are one principle"); (2) look for a more elegant
construction. All numbers are for the paper's Panda example (W = diag(0.01 I, I), coning circle, L = 0.77008765,
q_r = (0, -pi/4, 0, -3pi/4, 0, pi/2, pi/4)), N = 500 path nodes, Euclidean joint metric K = I unless "K = M" is said
(K = M uses the MuJoCo-menagerie Panda mass matrix with armature 0.1, through the lead's `num/pandakin.py::MjDyn`).
Where both exist, my numbers reproduce the lead's (`num/run_harmonic.out`: E = 6.669474 for K = I and 3.257517 for
K = M at N = 500, the same minimizers as my A and my near-rest K = M minimizer).
Scripts and raw outputs: `harmonic_lift/adv/` (Section 6). My kinematics are independent of the lead's and agree with
the C++ library to 4.4e-16 (`adv/check_kin.py`).

Status tags used throughout:
**[P]** proven here (a proof or a complete proof sketch is given) or a standard theorem with a reference;
**[PH]** proven under the hypotheses stated next to it;
**[N]** numerical evidence on the paper's path (not a proof);
**[C]** conjecture;
**[F]** false as stated (a counterexample or counter-argument is given).

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict

1. **Theorem A (frames) is right, after three precisions [P].** (i) The transport (eq:transport) is horizontal for the
   *left-trivialization* connection of TSE(3) projected onto zeta_gamma^perp (metric for W, with torsion), not for the
   Levi-Civita connection of W. (ii) The minimizers of (1/2) oint |N^T W N'|^2 over closed W-orthonormal normal frames
   are exactly N_0 exp(-(s/L) Phi) Q_0 with Q_0 in O(5) arbitrary and Phi a minimal-norm logarithm of H; they are unique
   modulo Q_0 iff -1 is not an eigenvalue of H (if it is, there are two, or a pair of 2-spheres, of minimizing families).
   (iii) The "one principle" works for frames because the horizontal part of their Dirichlet energy is gauge invariant
   (it equals oint |zeta_gamma'|_W^2 per column); for the lift the horizontal energy is a potential that varies along
   the fibre, and that is exactly what destroys uniqueness. Numerically, N^T W N' along the paper's frame is constant to
   4.6e-11 with spectrum {+-0.16118 i, +-8.008e-5 i, 0}/L = spec(-Phi)/L [N].
2. **The harmonic lift is not canonical on the paper's own example [N].** With K = I the energy has at least three strict
   local minima in one homotopy class of closed sections: E = 5.128, 6.669, 6.784 (Hessians positive definite). With
   K = M: E = 2.577, 3.258, 3.332. The lowest-energy minimizer is 0.033 rad from joint 2's limit for K = I (below the
   0.10 rad safety margin), and violates joint 2's limit by 0.33 rad for K = M. The minimizer the lead's numerics report
   (6.669 for I, 3.258 for M) is the one in the basin of the rest-posture lift: its "random starts" (|u| <= 0.3 rad in a
   slice chart) never leave that basin. So the rest posture q_r re-enters through the initial guess.
3. **The physics is right but unflattering [P, N].** Euler-Lagrange <=> X^T(M q'' + C(q,q') q') = 0 and the d'Alembert
   reading (gravity-free passive arm dragged at constant speed by a flange wrench) are correct. But every nondegenerate
   minimizer is a *hyperbolic* periodic orbit of that passive motion (proof in 2.5); the computed Floquet multipliers are
   3.81, 5.21 and 4.14 per lap for the three K = I minimizers. The "natural motion" is an inverted-pendulum balance.
4. **Several sub-claims need repair.** The bundle claim fails for full self-motion circles on this path (the circle
   passes within sigma_min = 0.0073 of a singular configuration and its homology class in T^7 changes twice along the
   loop) [N]; with joint limits, the relevant object is a trivial *interval* bundle [N]. The "isometric transport =>
   uniform drift" formula is false as stated [F] and, under the hypothesis that does make it true, gives a whole circle
   of minimizers (non-uniqueness). Aubry-Mather needs circle fibres and is inapplicable here. "Closure is built in" is
   equally true for the pointwise argmin lift (2.2). The k -> infinity statement is right [PH], with error about 19/k
   rad, reached only for k >~ 1e4 [N].
5. **What beats it.** The logarithm of the holonomy diffeomorphism fails (Diff(F) has no locally surjective
   exponential; on this path h is not even a diffeomorphism of the full fibre circle), but *why* it fails points to the
   fix: **reduce the structure group of the self-motion bundle from Diff(F) to the measure-preserving group
   Diff_mu^+(F) = {id} (interval fibres) or U(1) (circle fibres), then cancel the holonomy uniformly (Coulomb gauge)**.
   This is literally what Lemma 1 does with SO(5) and Lemma 3 with SO(r). For the Panda (interval fibres, bounded by
   joint 2 at both ends for every s [N]) it yields a lift that is canonical in the strongest sense: the **midsection**,
   the point bisecting the admissible self-motion arc in joint-space arc length (the unique section fixed by the fibre
   reflection). It is unique and closed by construction, smooth under a checkable transversality hypothesis,
   continuous in all data (no Maxwell walls), joint-limit aware, and parametrization free [PH]. On the paper's path:
   smallest joint margin 0.3623 rad (rest lift 0.3592), sigma_min >= 0.208, within 0.48 rad of the rest lift [N].
6. **Recommendation (Section 4):** midsection lift + the unified "reduce, then Coulomb-gauge" principle for N, X and q*;
   keep the rest-posture lift as the posture-preferring variant, stated as a fibrewise Morse minimum (closure automatic;
   f restricted to every sampled fibre is a perfect Morse function [N]); keep the energy E_k only as a remark that
   connects the pointwise lift (k -> infinity) to the passive-arm motion (k = 0), without claiming canonicity.

---------------------------------------------------------------------------------------------------------------------

## 1. Facts about the paper's path that decide most questions

| Quantity | Value | Script |
|---|---|---|
| Full self-motion fibre through q*_rest(s) on T^7 (limits ignored), Euclidean length | 11.95 to 15.98 rad (20 sampled fibres) | `selfmotion.py`, `fibre_topology.py` |
| Its class in H_1(T^7) = Z^7 | (-1,0,0,0,1,0,0) for s <= 0.239 and s >= 0.616; (-1,0,0,0,0,0,-1) for 0.254 <= s <= 0.539 | same |
| Smallest sigma_min(J_b) on that fibre | 0.0073 at s = 0.2464, at q ~ (-3.09, -2.29, -0.08, -2.34, -1.53, -0.19, -3.79): q2 outside its limit | `fibre_topology.py` |
| Admissible part (open joint box) of the fibre containing q*_rest | for all 500 nodes the admissible arc through q*_rest has **both ends at joint 2's lower limit -1.7628**; on the 11 fully traced circles the fibre leaves the box exactly once, so its admissible part is a single arc; arc length 5.50 to 8.58 rad; q*_rest at offset 0 in (a(s), b(s)), a in [-4.29, -2.76], b in [2.48, 4.37] | `midpoint.py` |
| f = (1/2)\|q - q_r\|^2 on each sampled fibre circle | exactly one local minimum (q*_rest) and one local maximum | `checks.py` |
| Euclidean holonomy h of the full circle over one lap | h(theta) - theta in [-0.649, -0.197] rad for 139 of 140 samples, no fixed point; one sample's horizontal lift passed the near-singular configuration and diverged | `holonomy.py` |
| Pseudoinverse from q*_rest(0) | gap 0.24 rad after one lap, 4.9 rad after 20 laps, no convergence | `../drift/drift` |
| Harmonic minimizers, K = I | A: E = 6.6695, offset +0.64..+0.87 rad from q*_rest, margin 0.3596 (j4); B: E = 5.1283, offset -2.62..-2.12, margin 0.033 (j2); C: E = 6.7838, offset +1.97..+2.40, margin 0.062 (j2). All Hessians positive definite (lambda_min 4.5e-3, 7.0e-3, 5.2e-3) | `run_harmonic.py` |
| Harmonic minimizers, K = M | E = 3.2575 (near rest, margin 0.3614); 2.5772 (q2 down to -2.09: infeasible); 3.3322 (q2 down to -2.31: infeasible) | `harmonic_M.py` |
| Rest-posture lift energy | E_I = 7.3397, E_M = 3.2987 | |
| Midsection lift | E_I = 9.167, margin 0.3623 (j4), sigma_min >= 0.2082, max \|q'\| = 5.74 (rest 5.28), offset from rest in [-0.48, 0.35] | `midpoint.py`, `checks.py` |
| Sensitivity to d5 += 1 mm (path fixed) | max node shift: rest 4.8e-3, harmonic A 6.3e-3, harmonic B 4.6e-3, midsection 4.8e-3 rad | `sensitivity.py` |
| Penalized family E_k = E + (k/2) oint\|q - q_r\|^2, distance to rest lift | k = 10: 0.221; 1e2: 0.0561; 1e3: 0.0140; 1e4: 1.82e-3 (k * dist = 2.2, 5.6, 14.0, 18.2) | `penalty.py` |
| Frame connection form L * N^T W N' | constant along s to 4.6e-11; eigenvalues +-0.16118 i, +-8.008e-5 i, 0 | `checks.py` |

---------------------------------------------------------------------------------------------------------------------

## 2. Critique of the working hypothesis

### 2.1 Is E_gamma = gamma^* Q_reg -> S^1 a fibre bundle with 1-dimensional fibres?

* **[P] Ehresmann.** If Q' is open, Pi|Q' is a submersion and proper, then E_gamma restricted to Q' is a locally trivial
  bundle over S^1. On Q = T^7 (no limits) properness holds, so the statement is true **iff gamma(S^1) avoids every
  critical value of Pi: T^7 -> SE(3)**, i.e. no configuration in *any* component of *any* fibre over the path is
  singular. Checking the lift (sigma_min >= 0.207 along q*) is not enough.
* **[N] It fails on the paper's path.** The fibre circle through q*_rest(s) changes its class in H_1(T^7) between
  s = 0.2387 and 0.2541 (joint 5 winding +1 is replaced by joint 7 winding -1) and back between s = 0.539 and 0.616,
  and its sigma_min drops to 0.0073 at s = 0.2464 (at q6 ~ -0.19). Fibres over a segment free of
  critical values are isotopic, hence homologous, so (if the traces are right) gamma crosses the critical-value
  hypersurface at least twice: the full fibre undergoes a Morse surgery. Either way the fibre passes within 1e-2 of a
  singularity, so "circle bundle" is not a usable model here.
* **[N] With joint limits the useful object is a trivial interval bundle.** The admissible part (open joint box) of the
  fibre component through the lift is a single open arc whose two ends lie on joint 2's lower limit for every node; the
  near-singular configuration lies outside (q2 = -2.29). So E_gamma^adm = {(s, theta): a(s) < theta < b(s)}.
* **[P] In general, with limits, E_gamma is not a bundle.** Pi restricted to the open box is a submersion but not proper.
  Local triviality near s holds when each admissible arc ends at transversal crossings of a single face (implicit
  function theorem); at tangencies of a self-motion curve with a face, arcs are born, die, split or merge (a codimension
  1 event in s for a generic path, by transversality). So "sections of a bundle over a circle" is a hypothesis on the
  path, of the same nature as Assumption 1.
* **[P] Monodromy can permute fibre components.** If a fibre has several circle components (7R arms can have several
  self-motion manifolds over one pose, cf. Burdick 1989), the monodromy around gamma can permute them; components in a
  non-trivial cycle carry no section at all (this is the redundant analogue of "cuspidal" 6R behaviour). The fibres are
  oriented globally on Q_reg by the generalized cross product (because J_b maps E^perp isomorphically onto se(3) = R^6,
  which orients E^perp, and Q is orientable), so the monodromy preserves orientation and the components carrying sections are tori, never Klein
  bottles.

### 2.2 Is closure "built in"?

* Only in a trivial sense. Varying over closed sections presupposes that one exists (that is Assumption 1), and the
  existence of a minimizer needs compact fibres or a barrier (2.9).
* **[P] The pointwise lift is closed by construction too.** If for every s the function f restricted to the admissible
  fibre has a unique global minimizer, nondegenerate and away from the fibre's ends, then s -> argmin is a smooth
  section (implicit function theorem plus uniqueness) and q*(L) = argmin over F_L = argmin over F_0 = q*(0). The
  paper's continuation and closure check are a numerical certificate that the tracked local minimum is the global
  one, not a defect of the definition. On this path f is a perfect Morse function on every sampled fibre (one minimum,
  one maximum) [N] and the second-difference test in the design doc is positive at all 2000 nodes, so the rest-posture
  lift is canonical given q_r and automatically closed.
* **The real advantage** of a variational definition is different: it is insensitive to fibrewise Cerf events of a
  pointwise criterion (a fold where the tracked minimum dies, a Maxwell point where two minima exchange) [PH: the
  energy is smooth on the section space and does not see them]. The price: the lift becomes path-global and
  parametrization dependent (2.3), and the non-uniqueness moves from s to the space of paths (2.4).

### 2.3 Dependence on K, on W and on the parametrization

* **K.** For K = I minimizer A is 0.76 rad (mean) from the rest lift; for K = M the corresponding minimizer is 0.12 rad
  away [N, mine and the lead's agree]. K = M depends on payload, tool, armature and a dynamics model that the paper does
  not otherwise need offline (identified Panda dynamics are out of scope in the repo). A kinematic object (the lift)
  becoming payload dependent is a modelling smell, not a feature.
* **Parametrization [P].** E(q) = (1/2) oint |q'|^2 ds is not invariant under reparametrization of the base circle; a
  non-affine change s = sigma(t) changes the minimizer. Using W-arc length ties the lift to W, i.e. to the tube's
  length scale l: a tube-design parameter leaks into the redundancy resolution. On the coning circle this is invisible
  because the W-speed in theta is constant for every product W, so s is proportional to theta for every l. On a
  generic path it is not.
* The parametrization-free version, the length oint |q'|_K ds, has linear growth (not Tonelli): minimizing sequences
  can develop vertical segments (the base pauses while the arm moves along the fibre), so a length minimizer need not
  be a section [C, standard behaviour of non-superlinear problems].
* Pointwise lifts (rest posture, midsection) depend on the pose gamma(s) only: they are restrictions to C of a local
  inverse-kinematics section Sigma of Pi, independent of parametrization and of W [P].
* **Speed [P].** The passive-arm reading needs constant speed: for q(t) = q*(nu(t)), X^T(M q-ddot + C q-dot) =
  nu-dot X^T M q*' != 0 in general.

### 2.4 Uniqueness; minimum versus saddle

* **[N] Three strict local minima in one homotopy class**, for both metrics (table in Section 1). The class is the same
  because each was reached from the rest lift by a continuous deformation (flowing every node along its self-motion).
  Minimizers B and C (K = I) and both low-energy K = M minimizers sit at or beyond joint 2's lower limit.
* **[PH] Saddles.** Two strict local minima in a connected section space force a third, mountain-pass critical section
  (finite-dimensional discretization; Palais-Smale holds on the compact sublevel sets). So the problem has at least five
  critical sections in this class.
* **[P] Aubry ordering** (crossing lemma, 1-DOF time-periodic Tonelli): two *global* minimizers in the same class do not
  cross. It says nothing about local minimizers, which is what continuation finds.
* **[P] Generic uniqueness, and why it does not help.** For a generic potential the global minimizer is unique (Mane
  1996; Bernard-Contreras 2008). But the map path -> global minimizer is discontinuous across codimension-1 walls in
  path space where two minima have equal energy. Here the global minimizer B is excluded by the 0.10 rad safety margin;
  if the rule is "the lowest-energy minimizer that respects the margin", the competitors are A and C, separated by
  E_C - E_A = 0.114 (1.7 %): a modest change of the path or of K can swap them, and the lift then jumps by 1.3 to 1.5
  rad along the self-motion. This is the Maxwell disease the harmonic construction was meant to cure, moved from s to
  the space of paths.
* **[P] The degenerate case.** If horizontal transport is isometric *and* the horizontal energy is constant along fibres,
  every curve theta(s) = theta_0 - s tau / L is a minimizer: a whole circle of minimizers (2.7). Exactly like the O(5)
  gauge freedom of N, but here nothing fixes theta_0.
* **What the lead's numerics actually compute:** the local minimizer in the basin of the rest-posture lift (random
  starts within |u| <= 0.3 rad of it). That is "the energy-relaxation of the rest-posture lift", i.e. q_r is still the
  selector.

### 2.5 The physical interpretation

* **Euler-Lagrange [P].** Variations of a section are vertical, delta q in ker J_b(q). The first variation of
  E = (1/2) oint q'^T M(q) q' ds is -oint delta q^T (M q'' + C(q, q') q') ds, where only the vector C(q, q') q' (the
  Christoffel term) matters, so any Christoffel-consistent C works. Hence X^T(M q'' + C q') = 0, i.e. the M-orthogonal
  projection of nabla_{q'} q' onto ker J_b vanishes. For K = I this is X^T q'' = 0.
* **d'Alembert [P].** For q(t) = q*(nu t), M q-ddot + C q-dot = nu^2 (M q*'' + C(q*, q*') q*') annihilates ker J_b, so it
  lies in range J_b^T: it is produced by a flange wrench alone (ideal, workless, rheonomic constraint; frictionless,
  gravity-free or gravity-compensated arm, constant nu, any nu by homogeneity). In the paper's controller at steady state
  on the lift this reads: u - G = J_b^T f, zero null-space torque.
* **Stability [P]: nondegenerate minimizers are hyperbolic.** In a fibre chart the second variation is the Hill
  (Sturm-Liouville) form Q[y] = oint (a y'^2 + c y^2) ds with L-periodic a > 0, c. A nondegenerate minimizer has
  Q > 0 on periodic y, i.e. lowest periodic eigenvalue lambda_0 > 0. By Floquet theory lambda_0 is the bottom of the
  spectrum of -(a y')' + c y on L^2(R), so 0 lies strictly below the spectrum: the Jacobi equation has no bounded
  solution and has a positive one, so its monodromy has real positive multipliers mu > 1 > 1/mu. The Jacobi equation is
  the variational equation of the passive dragged-arm dynamics along the orbit, so the orbit is unstable. Computed from
  the discrete Jacobi operator: mu = 3.81 (A), 5.21 (B), 4.14 (C) per lap [N] (K = I; for K = M the same theorem
  applies to the physically meaningful case).
* **The penalized energy [P].** Minimizing oint (T + V) is extremizing the action of the Lagrangian T - (-V): the term
  (k/2)|q - q_r|^2 is a *repulsive* spring in the physical reading. "The arm balanced on a potential hill" again.

### 2.6 Theorem A, checked line by line

Notation: zeta = zeta_gamma, |zeta|_W = 1, N_0' = -zeta zeta'^T W N_0.

1. **[P] Invariants.** d/ds(zeta^T W N_0) = zeta'^T W N_0 - (zeta^T W zeta) zeta'^T W N_0 = 0;
   d/ds(N_0^T W N_0) = -(N_0^T W zeta) zeta'^T W N_0 - (...)^T = 0. So W-orthonormality and normality persist.
2. **[P] Horizontality.** N_0^T W N_0' = -(N_0^T W zeta) zeta'^T W N_0 = 0. Precisely: N_0' is parallel to zeta, so the
   W-normal component of N_0' vanishes: N_0 is parallel for the connection nabla^perp n = P^perp n' on the normal
   bundle zeta^perp, where ' is the derivative of body coordinates. That is the left-trivialization (Cartan-Schouten
   minus) connection of TSE(3), which is W-metric but has torsion; "rotation minimizing" is with respect to it, not to
   the Levi-Civita connection of W (which would add (1/2)(ad_zeta y - W^{-1}(ad_zeta^T W y + ad_y^T W zeta)) and give a
   different frame).
3. **[P] Energy reduction.** Closed frames are N = N_0 Q, Q in H^1([0, L], O(5)), with N(L) = N(0) <=> Q(L) = H^T Q(0)
   (H = N_0(0)^T W N_0(L)). Then N^T W N' = Q^T Q' and E(N) = (1/2) oint |Q'|_F^2.
4. **[P] Minimizers.** For curves with Q(L) = H^T Q(0), E >= d(Q(0), H^T Q(0))^2 / (2L) with equality iff Q is a
   constant-speed minimizing geodesic of the bi-invariant metric; by bi-invariance d(Q_0, H^T Q_0) = |Phi|_F for a
   minimal-norm logarithm Phi of H, independent of Q_0. Geodesics are Q(s) = Q_0 exp(sB) with exp(LB) = Q_0^T H^T Q_0,
   so Q(s) = exp(-(s/L) Phi) Q_0 and N = N_0 exp(-(s/L) Phi) Q_0, Q_0 in O(5) free, energy |Phi|_F^2/(2L) =
   (theta_1^2 + theta_2^2)/L. All other logarithms Phi' of H give the remaining closed critical frames (the boundary
   term Q^T Q'(L) = Q^T Q'(0) holds automatically for one-parameter subgroups), with larger energy.
5. **[P] The eigenvalue -1 case.** H in SO(5) has the eigenvalue -1 with even multiplicity 2m. Minimal logarithms
   restricted to that eigenspace are pi J with J an orthogonal complex structure on it: two choices for m = 1, the
   space O(4)/U(2) (two 2-spheres) for m = 2. So the minimizing frame is unique modulo Q_0 iff -1 is not in spec H; the
   set {H has eigenvalue -1} is the cut locus of I in SO(5), a codimension-1 set, and across it the minimal frame jumps
   (as a function of the path). The paper's principal logarithm is fine for its examples (angles 0.161, 8.0e-5 for the
   Panda; 2.166, 0.492 for the rigid body). Rotation angles exactly pi would need a stated tie-break.
6. **[P] Coulomb gauge.** Connections on the trivial principal G-bundle over S^1, modulo based gauge, are classified by
   their holonomy; A = N^T W N' constant (d^*A = 0) is the Coulomb gauge; constant gauges correspond to logarithms of
   the holonomy, minimal |A| to minimal logarithms. Numerically A is constant along the paper's frame to 4.6e-11 with
   spec(L A) = {+-0.16118 i, +-8.008e-5 i, 0} [N].
7. **[P] Lemma 3 is the same statement.** With P the K-projector onto E^perp, (I - P) P' (I - P) = 0 implies that
   D_s x = -(D_s P) x keeps x in E and has (I - P) D_s x = 0 (horizontal); |D_s X|^2 = |(D_s P) X|^2 + |X^T K D_s X|^2 and
   the first term depends only on X X^T K (the projector onto E), so it is gauge invariant; minimizers are
   X^0 exp(-(s/L) Phi) Q_0, unique modulo O(r) iff -1 is not in spec H.
8. **[P] What "one principle" means, precisely.** The full Dirichlet energy of N as a map into (R^{6x5}, W-Frobenius) is
   oint (|zeta'|_W^2 + |N^T W N'|^2): minimizing it is minimizing the vertical (Wood-type) energy, because the horizontal
   part is gauge invariant. For the lift the full energy is oint (U(s, q) + |q'^V|^2) with U = zeta^T (J_b K^{-1} J_b^T)^{-1}
   zeta, which varies along the fibre. Hence: the frame problem is a geodesic problem on a compact Lie group with an
   explicit, essentially unique solution; the lift problem is a non-convex 1-DOF periodic Lagrangian problem with
   several critical points. The analogy is formal, not structural. Also note the terminology: in C. M. Wood's sense a
   "harmonic section" is a critical point of the *vertical* energy (1/2) oint |q'^V|^2; for the lift that is a different
   functional, whose minimum is 0 exactly when the holonomy of the K-orthogonal connection has a fixed point (it does
   not on this path: 2.1 and (a)).

### 2.7 "If horizontal transport maps fibres isometrically, the harmonic lift is the horizontal lift plus a uniform drift"

* **[F] as stated.** In isometrically transported fibre coordinates theta, E = (1/2) oint (theta'^2 + U(s, theta)) ds and
  the Euler-Lagrange equation is theta'' = (1/2) U_theta. The uniform drift theta_0 - s tau / L solves it only where
  U_theta = 0 along it. Isometric transport (implied by totally geodesic fibres: for the K-orthogonal connection,
  L_{X^h} K on vertical vectors equals -2 K(X^h, T(.,.)) with T the fibres' second fundamental form) says nothing
  about U.
* **[PH] Corrected statement.** If transport is isometric and U is constant on fibres (Pi is a Riemannian submersion with
  totally geodesic fibres onto (SE(3), g_U)), the minimizers are theta(s) = theta_0 - s tau / L with tau the minimal
  representative of the holonomy shift and **every** theta_0: a circle of minimizers, not one.
* **[N] Neither hypothesis holds for the Panda:** the fibre lengths vary from 11.95 to 15.98 rad along the path, which
  is impossible under isometric transport. The lead's numerics nevertheless find the uniform-drift curve within 0.045
  rad of minimizer A; that is a quantitative coincidence on this path, not a theorem.

### 2.8 Aubry-Mather

* **[P] Where it applies** (circle fibres, e.g. continuous-rotation joints and a path free of critical values): in a
  global trivialization of the torus E_gamma the reduced Lagrangian (1/2)(a theta'^2 + 2 b theta' + c) is Tonelli and
  L-periodic in s, so Aubry-Mather theory applies to it directly (Mather's variational theory; Bangert 1988 for the
  1-DOF periodic case): minimizers exist in every rotation class, global ones are ordered, and the variational form of
  Poincare-Birkhoff gives at least two closed critical sections per class (a minimizer and a mountain pass). Its
  time-L map is a monotone twist map when no conjugate points occur within time L (Moser 1986 showed every twist map
  arises this way); in that case minimizers have Greene residue <= 0 (MacKay-Meiss). Nondegenerate minimizers are
  hyperbolic in any case (2.5).
* **[P] "Winding number" is not intrinsic to the bundle.** Homotopy classes of sections of a torus bundle over S^1 form
  a Z-torsor; an integer label needs a reference section. On T^7 there is a canonical label, the joint winding vector in
  H_1(T^7) = Z^7, and inside the joint box every lift has winding 0.
* **[N] It does not apply here.** The admissible fibres are intervals, and the full circles undergo surgery at a
  critical value (2.1). With interval fibres there is one homotopy class and no rotation number; what survives is the
  1-DOF ordering of global minimizers and the mountain-pass structure of 2.4. So Aubry-Mather explains the
  multiplicity of closed lifts rather than selecting one.

### 2.9 Joint limits and singularities

* **[P]** A barrier (or a complete metric K_B = exp(2B) K on the admissible set) gives existence of minimizers in each
  class. But a complete metric changes the energy everywhere, not only near the limits, and the barrier's profile then
  decides which local minimum is global: on this path the global minimizers without barrier (B for K = I, both
  low-energy ones for K = M) live at or beyond joint 2's limit [N], so the barrier is doing the selecting [C: the
  constrained global minimizer hugs joint 2's barrier].

### 2.10 The rest-posture lift as k -> infinity

* **[PH]** If f restricted to each fibre has a unique, nondegenerate global minimum q_f(s) away from the fibre's ends,
  minimizers of E_k near q_f exist for large k and q_k = q_f + (1/k) (Hess^V f)^{-1} (nabla_{q_f'} q_f')^V + O(k^{-2})
  (from the Euler-Lagrange equation grad^V f = (1/k) (nabla_{q'} q')^V; regular perturbation: the linearization
  -d^2/ds^2 + k Hess^V f + O(1) is positive and invertible with inverse O(1/k) for large k, so a contraction argument
  gives the solution near q_f).
  Numerically k * dist = 2.2, 5.6, 14.0, 18.2 for k = 10, 1e2, 1e3, 1e4: the 1/k regime (constant about 19 rad) is
  reached only for k >~ 1e4 [N].
* **[C]** At fibrewise Maxwell or fold points of f the convergence is non-uniform, with interior layers of width
  O(k^{-1/2}) (heteroclinic "kinks"); not proven here.
* The pseudoinverse with null-space pull (gap 0.36 / 0.030 / 1.3e-3 rad for k = 1 / 10 / 100) is a different,
  dissipative family. For r = 1 its lap map is an increasing diffeomorphism of an interval onto its image, so every
  bounded orbit converges monotonically to a fixed point, i.e. to a closed lift (Massera 1950) [P]. It needs no
  variational principle to close.

### 2.11 Cost and robustness against the pointwise lift

| | Rest-posture lift | Harmonic, K = I | Harmonic, K = M | Midsection (Section 3, h1) |
|---|---|---|---|---|
| Offline work | one Newton solve per node (2-3 updates, square system) | periodic BVP; ~5 Newton steps with a cyclic tridiagonal Hessian | same, plus M and its derivatives | two exit searches per node along the self-motion (8 s in batched Python for 500 nodes) |
| Uniqueness | unique global minimum per fibre [N] | >= 3 local minima [N] | >= 3, two infeasible [N] | unique by construction [P] |
| Smallest joint margin | 0.3592 | 0.3596 (A), 0.033 (B), 0.062 (C) | 0.3614 (near rest) | 0.3623 |
| Energy E_I | 7.340 | 6.669 (A), 5.128 (B) | | 9.167 |
| max \|q'\| | 5.28 | 5.02 (A) | 5.41 (lead) | 5.74 |
| Shift for d5 + 1 mm | 4.8e-3 | 6.3e-3 (A) | | 4.8e-3 |
| Depends on | q_r | K, parametrization (W), initial guess | + payload, dynamics model | joint box, fibre measure |

Conditioning is not the issue: the harmonic minimizer is about as sensitive to a kinematic perturbation as the
pointwise lifts. Non-uniqueness and the hidden selector are the issue.

---------------------------------------------------------------------------------------------------------------------

## 3. Alternatives

Throughout r = 1 unless stated; F_s is the admissible fibre over gamma(s); x(q) the unit kernel vector (generalized
cross product), which orients all fibres.

### (a) Logarithm of the holonomy diffeomorphism

**Statement [P].** Let H be an Ehresmann connection on E_gamma with complete horizontal lifts, Phi_s: F_0 -> F_s its
transport and h = Phi_L in Diff^+(F_0) its holonomy. If Z is a vector field on F_0 with exp(Z) = h (time-1 map), put
G_s = Phi_s o exp(-(s/L) Z): F_0 -> F_s. Then G_L = id, so s -> G_s(p) is a closed lift for every p, (s, p) -> G_s(p)
trivializes E_gamma = S^1 x F_0, and the closed lifts are the integral curves of hor - (1/L) Z_s with
Z_s = (Phi_s)_* Z (using (exp(tZ))_* Z = Z). For a linear holonomy H = exp(Phi) on R^5 this is exactly
N_0 exp(-(s/L) Phi).

**When does a logarithm exist? [P, literature]** For F = S^1:
* exp: Vect(S^1) -> Diff(S^1) is not locally surjective (Hamilton 1982, Milnor 1984): there are diffeomorphisms
  arbitrarily close to id on no one-parameter subgroup. Diffeomorphisms that embed in a flow are meagre (Palis 1974).
* Irrational rotation number: h embeds in a smooth flow iff h is C^infinity-conjugate to a rotation (the flow of a
  nowhere-zero field is conjugate to a rotation flow). Herman (1979) and Yoccoz (1984): true for every smooth h with
  Diophantine rotation number; false for some Liouville numbers. The logarithm is then unique up to an integer winding
  k (the centralizer of an irrational rotation is the rotations).
* Rational rotation number p/q, q >= 2: h embeds only if h^q = id.
* Fixed points: on each complementary interval, Szekeres (1958) gives a C^1 vector field (unique among C^1 fields by
  Kopell's lemma), C^infinity inside, smooth at hyperbolic fixed points (Sternberg), with Takens-type obstructions at
  degenerate ones; but the fields obtained from the two ends of an interval agree only if Mather's invariant is trivial,
  which fails generically. So generic h with fixed points does not embed.
* In a one-parameter family of paths the rotation number moves continuously and locks to rationals on open sets (Arnold
  tongues): **existence of a smooth logarithm is a nowhere-continuous property of the path.**

**Panda [N].** On the full fibre circle (limits ignored, K = I) h has no fixed point (displacement -0.649 to -0.197 rad
per lap, rotation number about -0.027 turn per lap), and one sample's horizontal lift diverged because the full circle
passes the near-singular configuration: h is not even a diffeomorphism of the full circle on this path. With joint
limits the admissible arc is not invariant (points drift by 0.2-0.65 rad per lap towards and through joint 2's limit),
so h is only a partial map and "log h" is meaningless.

**Even when it exists** it yields a trivialization, i.e. an F_0-family of closed lifts; a selector (a point p) is still
needed, the analogue of Q_0.

**Verdict.** Most elegant formally, the exact mirror of Lemma 1, and unusable: not robust, not canonical, fails on the
paper's path. **But the failure is instructive**: Lemma 1 works because its structure group SO(5) is a compact Lie group
with surjective exponential. Reducing the self-motion bundle's structure group to such a group is the fix: see (h1).

### (b) Flat connections, transverse foliations, global redundancy coordinates

* **[P] In E_gamma (a surface) every transverse line field is integrable**; Frobenius is vacuous there. A flat
  connection is the same thing as a transverse foliation, and "closed lift" means "leaf with trivial holonomy".
* **[P] Global redundancy coordinates.** If phi: U -> R (U open in Q_adm) has d phi(x) != 0 on U, then each level set
  Sigma_c = phi^{-1}(c) is a 6-manifold transverse to the self-motions, Pi|Sigma_c is a local diffeomorphism, and
  q*(s) = (Pi|Sigma_c)^{-1}(gamma(s)) is closed by construction whenever it is single valued (e.g. phi monotone on each
  admissible arc). Examples: the coupling phi = q3 - q1 (the repo's coupled lift), phi = q7 (the parameter of He-Liu's
  analytical Panda IK), the swivel/arm angle of SRS arms (Kreutz-Delgado, Long, Seraji 1992), and **phi = d f(x) =
  x(q)^T (q - q_r): the rest-posture lift is the level set {x(q)^T (q - q_r) = 0}** (the augmented/extended Jacobian of
  Baillieul 1985). Algorithmic singularities are where d phi(x) = 0, i.e. fibrewise folds.
* **[P] Cyclicity theorems.** A velocity-level resolution is cyclic for all loops iff its horizontal distribution in Q is
  integrable (Shamir-Yomdin 1988), iff it is the derivative of an inverse-kinematic function (Baker-Wampler 1988): the
  cyclic schemes are exactly the foliations of Q transverse to the self-motions, i.e. the families Sigma_c. The 0.24 rad
  drift shows the K = I pseudoinverse distribution is not integrable [N].
* **Controller bonus [P, cost C].** The redundancy output can be the defining function itself, beta = phi(q) (as the
  coupling output already is): no eta-dependence, no q*(eta), X(eta) splines, no X' term in d beta, and Proposition 2
  needs only d phi(x) != 0 on the lift. For phi = x^T(q - q_r), d phi needs the derivative of the generalized cross
  product (cofactor derivatives); whether this fits the 24 us budget is untested [C].
* On this path q2 is not a global coordinate on the admissible arcs: both ends are at q2's *lower* limit, so q2 rises
  and falls along each arc [N].
* **Verdict.** Not a new lift, but the right *form* for every pointwise lift: "a lift is the trace on C of a transversal
  Sigma to the self-motion foliation". Robust, cheap, closure automatic. Elegance moderate; practicality high.

### (c) Fibrewise Morse / Morse-Bott functions and Cerf theory

* **[P] Statement.** For generic f the fibrewise critical set Sigma_f = {q: d f(x(q)) = 0} is a smooth 6-manifold
  (transversality of the vertical differential to the zero section), and Pi|Sigma_f is a Lagrangian-type map whose fold
  values (caustic) and Maxwell set (two global minima of equal value) are codimension-1 in SE(3). By Cerf (1970) a
  generic one-parameter family of functions on a 1-manifold has only birth-death (A2) and Maxwell events. **The global
  argmin lift is smooth and closed iff gamma avoids the Maxwell set and the caustic of the global branch; a continuation
  lift closes iff gamma avoids its branch's fold and the branch has trivial monodromy.**
* **[N] Panda:** f = (1/2)|q - q_r|^2 restricted to each sampled fibre has exactly one minimum and one maximum (a perfect
  Morse function on the circle), so gamma is far from both bad sets.
* **Verdict.** The correct diagnostic and failure criterion for every pointwise lift, and a clean way to state
  Assumption 1 for the rest-posture lift. Not itself a new construction.

### (d) Harmonic map heat flow / curve shortening on sections

* **[P]** The L^2 gradient flow of E_k on sections, partial_t q = (nabla_s q_s)^V - k (q - q_r)^V, exists for all time
  on compact fibres (semilinear parabolic on S^1, no bubbling in dimension 1; cf. Eells-Sampson, Ottarsson 1985). The
  Panda FK and the path are real analytic, so the discretized energy is analytic in the fibre chart and the
  Lojasiewicz gradient inequality makes every bounded trajectory converge to a single critical section (Lojasiewicz;
  Simon 1983 in the infinite-dimensional case) [PH: analyticity of the chart, compact fibres or barrier].
* **[N]** The limit depends on the start: from the rest lift it reaches A, from shifted starts B or C.
* **Verdict.** A method and an existence proof; as a selector it is "the critical section downhill from the initial
  guess". Not canonical.

### (e) Geometric phases and gauge theory

* **[P]** The per-lap drift of a velocity-level scheme is the holonomy of its Ehresmann connection on E_gamma (a
  geometric phase with values in Diff(F)); its curvature is Omega(X, Y) = -[X^h, Y^h]^V, which vanishes iff the scheme is
  cyclic (Shamir-Yomdin). Since gamma is contractible in SE(3), for small loops the phase is to leading order the
  curvature integrated over a spanning disk (Montgomery's falling cat, Shapere-Wilczek swimmers), but here the fibre is
  not a group orbit, so the connection is not principal and the phase is non-linear.
* **[P] Hannay-Berry averaging** (Marsden-Montgomery-Ratiu 1990) makes it principal: give each circle fibre its
  normalized arc length (a U(1)-torsor, oriented by x), average the connection over the U(1) action; the averaged
  connection's holonomy tau in U(1) is a "Hannay angle of the self-motion", and its Coulomb gauge
  theta(s) = theta_0 - (s/L) tau gives a U(1)-family of closed lifts. Exactly Lemma 1 with SO(5) replaced by U(1). See
  (h1), where the same thing arises from optimal transport.
* **Coulomb gauge for Diff(F)** is not defined (no bi-invariant metric on Diff(S^1); the right-invariant L^2 metric gives
  Burgers geodesics).
* **Verdict.** The right language for the drift (0.24 rad per lap for K = I, 0.39 for K = M per the lead) and, after
  averaging, a canonical U(1) Coulomb gauge for circle fibres. Irrelevant for interval fibres (no holonomy at all).

### (f) Aubry-Mather and twist maps

See 2.8. **Verdict:** beautiful and explanatory (why several closed lifts must exist, why minimizers are ordered and
hyperbolic), but it does not select a lift, and it needs circle fibres that this path does not have.

### (g) Symplectic or Hamiltonian reduction by the self-motion

* **[P]** Marsden-Weinstein or Routh reduction needs a symmetry group acting on Q and preserving the Lagrangian; the
  self-motion foliation is not the orbit foliation of any such action (the flow of x does not preserve M or K in
  general, and the fibres have different lengths). What remains is reduction *by the constraint*: the 1-DOF
  time-periodic Lagrangian on E_gamma used in 2.5-2.8.
* **[P] Equivariant special case.** If gamma is the orbit of a one-parameter group of base rotations that is a symmetry of
  the arm (a path circle centred on the joint-1 axis), any lift defined by an invariant criterion is equivariant, and a
  unique minimizer is equivariant. The coning circle is centred at (0.40, 0, 0.45), off that axis: no symmetry.
* **Verdict.** Inapplicable.

### (h) Other constructions

#### (h1) Reduce, then Coulomb-gauge: the measure-preserving (optimal-transport) connection and the midsection  **(recommended)**

Choose a smooth positive measure mu_s on each admissible fibre arc, depending smoothly on s. Canonical candidates, all
determined by K alone: the K-arc length ds_K (recommended); the Leray (coarea) measure ds_K / NJ, NJ = det(J_b K^{-1}
J_b^T)^{1/2} the manipulability, which is the measure of the disintegration of vol_K along Pi (canonical, but it pulls
towards singularities); or ds_K * NJ (pulls away from them).

**Theorem B (interval fibres) [PH].** *Hypothesis:* each admissible fibre F_s is a single open arc of finite mu-mass
whose ends are transversal crossings of the boundary of the joint box (on the Panda path: both ends on joint 2's lower
limit, for all 500 nodes [N]). Let u_s: F_s -> (0, 1) be the normalized cumulative mu_s-mass from the lower end
(orientation from x). Then:
1. (s, q) -> (s, u_s(q)) is a diffeomorphism E_gamma^adm -> S^1 x (0, 1). Proof: the ends move smoothly in s (implicit
   function theorem at transversal crossings), the cumulative mass is smooth in (s, q) with derivative the density
   (> 0) along the fibre, and u depends only on gamma(s), so s = L and s = 0 give the same map.
2. For each c in (0, 1), s -> u_s^{-1}(c) is a smooth closed lift. **Closure is automatic** because the construction is
   pointwise in the pose.
3. The orientation-preserving measure-preserving diffeomorphisms of an interval form the trivial group
   (F_mu(phi(x)) = F_mu(x) forces phi = id). Hence there is exactly one connection on E_gamma^adm whose transport pushes
   mu_s to mu_s' (the monotone rearrangement, which in any orientation-preserving coordinates is the 1-D
   optimal-transport map for every strictly convex cost), it is flat with **trivial holonomy**, and its leaves are the
   sections u = c.
4. Adding orientation reversal, the residual group is Z/2 = {u -> u, u -> 1 - u}; its unique fixed section is
   **the midsection u = 1/2**: the point that splits the admissible self-motion arc into two halves of equal mu-mass.
   With mu = ds_K, the midsection maximizes the smaller of the two self-motion distances to the joint limits.
5. Continuity in all data: the midsection depends continuously (smoothly, under the transversality hypothesis) on the
   path, the kinematics and the box; there are no Maxwell walls, because nothing is minimized: it is the centre of a
   one-dimensional convex set.
6. It is also a global fibre coordinate: u extends to every q in Q_adm whose pose has a single transversal admissible
   arc, and (Pi, u) is a diffeomorphism onto its image (a global inverse-kinematics chart with redundancy parameter u)
   [PH].

**Theorem C (circle fibres) [P].** If the admissible fibres are full circles (no limits on the self-motion, path free of
critical values), the orientation-preserving measure-preserving diffeomorphisms are the rotations of the normalized
coordinate u in R/Z, so E_gamma is a principal U(1)-bundle over S^1 (trivial, not canonically). Among connections whose
transport pushes mu_s to mu_s', the one of least vertical kinetic energy oint int_F |w|_K^2 d mu_s ds (w the vertical
correction to the K-horizontal lift) rotates each fibre at a weighted fibre-average of the horizontal rotation rate
(Hannay-Berry averaging); its holonomy is an angle tau, and the Coulomb gauge u(s) = u_0 - (s/L) tau~, tau~ the
representative in (-1/2, 1/2], gives a U(1)-family of closed lifts, unique modulo u_0 unless tau~ = 1/2. This is
Lemma 1 with SO(5) -> U(1) and Phi -> 2 pi tau~; the residual u_0 is the analogue of Q_0 and needs a gauge-fixing rule
(e.g. minimize oint f(q_{u_0}(s)) ds over u_0 in U(1), a one-dimensional problem).

**The unifying principle [P for the three instances].** *Reduce the structure group of the natural bundle along the loop
to a compact Lie group G using the geometry of its fibres, then cancel the holonomy uniformly in arc length (Coulomb
gauge).* Normal frame: G = SO(5) via W-orthonormality (Lemma 1). Kernel frame: G = SO(r) via K-orthonormality (Lemma 3).
Lift: G = Diff_mu^+(F), which is {1} for interval fibres (no holonomy to cancel: the midsection) and U(1) for circle
fibres (Theorem C). This also explains (a): Diff(F) itself is not a Lie group with a usable exponential, so its holonomy
cannot be "logged"; the measure reduction removes exactly the part that cannot. And Theorem A becomes the statement that
each Coulomb gauge is the minimal-energy closed trivialization; for r = 1 interval fibres the minimal-energy statement
is vacuous (G is trivial).

**Panda numbers [N].** Ends a(s) in [-4.29, -2.76], b(s) in [2.48, 4.37] rad from the rest lift, both on joint 2's lower
limit for all nodes; finite differences of a, b and of the midsection are smooth at the grid scale (max |3rd difference|
/ h^3 about 1e3, consistent with a''' of order a'' / (L / 2 pi), no kinks). Midsection: offset from the rest lift
-0.48 to +0.35 rad; smallest joint margin 0.3623 rad (joint 4; rest lift 0.3592; joint 4's margin is fixed by the pose,
as CLAUDE.md notes); sigma_min >= 0.2082; E_I = 9.167 (rest 7.340); max |q'| = 5.74 (rest 5.28); 1 mm d5 change moves it
4.8e-3 rad (rest 4.8e-3).

**Weaknesses (honest).**
* Hypothesis of a single transversal admissible arc: fails when the excursion outside the box closes (fibre becomes a
  full admissible circle: switch to Theorem C) or the arc splits (a self-motion curve touching a face from inside).
  Both are codimension-1 events in s; the paper would state the hypothesis as its Assumption 1 and check it numerically,
  as it now checks closure.
* The measure is a choice (arc length in the joint metric is natural because the limits are a box in joint
  coordinates and beta is measured in that metric, but it is a choice).
* Not energy optimal: E_I 25 % above the rest lift; max |q*'| 9 % higher.
* Offline cost: tracing each fibre to both limits (seconds in C++), versus milliseconds for the rest lift.
* r >= 2: Diff_mu(F) is infinite-dimensional for dim F >= 2, so the reduction trick no longer gives a Lie group; the
  natural generalization is the Riemannian centre of mass (Karcher mean) of mu_s on the r-dimensional admissible
  region, unique when that region lies in a strongly convex ball (Karcher 1977) [C: hypotheses not checked for any arm].
* Elegance is by canonicity and simplicity, not by a deep variational principle; a referee may find "the midpoint of the
  self-motion arc" plain. I consider that a feature.

#### (h2) Fisher-Rao / H-dot^1 interpolation on Diff(S^1)/Rot

For circle fibres one can also close the holonomy without a logarithm: Diff(S^1)/Rot with the right-invariant
H-dot^1 metric is isometric to an open subset of an L^2 sphere via phi -> 2 sqrt(phi') (Lenells 2007; Khesin, Lenells,
Misiolek, Preston 2013), so the minimal-energy path from id to h^{-1} is an explicit great circle of square-root
densities, always defined [P]. Elegant and always defined, but it still leaves a rotation (U(1)) gauge and is heavier
than Theorem C. **Verdict:** a curiosity for circle fibres.

#### (h3) The attracting cycle of a dissipative velocity resolution

For r = 1 any velocity-level scheme q' = hor + V(s, q) whose lap map keeps an interval invariant has a closed lift, and
every orbit converges monotonically to one (Massera 1950; increasing interval maps) [P]. The null-space-pull
pseudoinverse is such a scheme (closed after a few laps, gaps 0.36 / 0.030 / 1.3e-3 rad for k = 1 / 10 / 100 [N]).
Canonical given the scheme and its gain, stable (unlike the harmonic orbit), but the scheme is the choice. **Verdict:**
a good numerical method for any lift defined as a fixed point; not a definition.

#### (h4) Max-min joint margin

q*(s) = argmax over F_s of min_i (distance of q_i to its limits). Nonsmooth (a max of a min), and degenerate on this path:
the binding joint (4) barely changes along the self-motion (the midsection, 2.4 rad away in places, improves joint 4's
margin by only 0.003 rad). **Verdict:** reject.

---------------------------------------------------------------------------------------------------------------------

## 4. Ranking and recommendation

| Rank | Construction | Canonical? | Smooth closed lift? | Joint limits | Practical | Elegance |
|---|---|---|---|---|---|---|
| 1 | (h1) Reduce, then Coulomb-gauge; midsection for interval fibres | yes: unique fixed section of the residual Z/2; continuous in all data | yes, by construction, under the single-transversal-arc hypothesis | built in | yes (seconds offline) | high: one principle for N, X, q* |
| 2 | (b)+(c) Pointwise lift as level set of a global fibrewise function (rest posture) | given q_r; unique on this path | yes while gamma avoids Maxwell set and caustic | indirect | best (already implemented) | moderate |
| 3 | E_k family (penalized harmonic) | no; interpolates rest (k -> infinity) and passive arm (k = 0) | yes for large k (unique near q_f) | no | yes | moderate; good as a remark |
| 4 | Pure harmonic section (working hypothesis) | no: >= 3 minima, selection by initial guess | yes, but which one | no (global minimizer at or beyond joint 2) | yes for K = I; needs dynamics for K = M | high on paper, fails in practice |
| 5 | (e)/Theorem C, U(1) Coulomb gauge for circle fibres | up to U(1) | yes | no | where circle fibres exist (not here) | high |
| 6 | (h3) attracting cycle of a velocity scheme | given the scheme | yes | partly | yes | low |
| 7 | (d) heat flow | no (initial-guess dependent) | yes | via barrier | yes | method only |
| 8 | (f) Aubry-Mather | no | | | | explanatory |
| 9 | (a) log of holonomy in Diff(F) | no; generically does not exist | | | no | formally highest |
| 10 | (g) reduction | inapplicable | | | | |

**Recommended construction for the paper.**
1. *Lemma 1 and Lemma 3:* keep the constructions, add Theorem A as corrected (minimal-energy closed frames = Coulomb
   gauge; unique modulo O(5)/O(r) iff -1 is not an eigenvalue of H; connection = left-trivialization connection
   projected onto the normal space).
2. *Lift:* define the admissible self-motion bundle over the path; replace Assumption 1 by "every admissible fibre is a
   single open arc whose ends cross the joint-limit box transversally" (checked numerically, like closure now); define
   q*(s) as the midsection (Theorem B). State the unifying principle: the three objects are Coulomb gauges of the
   natural bundles after reducing the structure group to SO(5), SO(r) and Diff_mu^+(F) = {1}.
3. *Alternative lift:* the rest-posture lift, stated as the trace of the fibrewise-critical manifold
   {x(q)^T (q - q_r) = 0} on the path, closed by construction while gamma avoids its Maxwell set and caustic (Cerf);
   say that on the example f is a perfect Morse function on every fibre.
4. *Remark only:* the energy E_k, the passive-arm reading (with "constant speed", "gravity compensated", and
   "the orbit is hyperbolic"), and k -> infinity. Do not call the harmonic lift canonical; the paper's own example has
   several minimizers.

If the authors prefer to keep the rest-posture lift (no code change at all), item 3 alone already removes the
"closure is only checked" objection, and item 1 gives the paper its geometric unification for N and X.

---------------------------------------------------------------------------------------------------------------------

## 5. Gaps and open points (not papered over)

* The trace-based facts (fibre class change, sigma_min 0.0073) are numerical; a branch jump of the tracer near the
  near-singular point cannot be fully excluded, but either reading makes the full-circle bundle model unusable here.
* Whether the fibres over the path have further components (other self-motion manifolds) was not examined; the
  statements concern the component through the lift.
* Multi-start used 9 starts (K = I) and 5 (K = M) along the self-motion; there may be more critical sections.
* Smoothness of a(s), b(s) and of the midsection is verified at the 500-node grid scale and by the transversality
  argument, not between nodes.
* K = M results depend on the MuJoCo-menagerie model (not an identified Panda); only the qualitative conclusion
  (several minima, the low-energy ones outside the limits) is claimed.
* The Floquet multipliers come from the discrete Jacobi operator (N = 500); the hyperbolicity theorem itself is exact.
* The interior-layer claim for E_k at Maxwell points and the Karcher-mean generalization to r >= 2 are conjectures.

## 6. Files (all in `harmonic_lift/adv/`)

* `kin.py`, `kinb.py`, `check_kin.py`: independent kinematics (single and batched), validated against the C++ dump.
* `selfmotion.py`, `fibre_topology.py`: full self-motion circles, admissible arcs, classes in H_1(T^7), sigma_min.
* `holonomy.py` (+ `hol_dbg.py`): Euclidean holonomy map of the full circle over one lap.
* `harmonic.py`, `run_harmonic.py`: K = I closed-section energy minimization (Riemannian Newton, cyclic tridiagonal
  Hessian with fibre curvature), multi-start; outputs `harm_euclid_500*.npy`.
* `harmonic_M.py`: K = M multi-start (run with `/home/jaiden/Desktop/SE3_pathFollowing/.venv/bin/python`); outputs
  `harmM_500_start*.npy`.
* `midpoint.py`: admissible arcs a(s), b(s), bounding joints, midsection; output `midpoint_500.npz`.
* `checks.py`: Morse count of f on fibres, fibre winding, Theorem A check on the paper's frame, Floquet multipliers.
* `sensitivity.py`: 1 mm d5 perturbation. `penalty.py` / `penalty.out`: the E_k family. `q2check.py`: q2 along arcs.
