# Harmonic lifts: closed lifts of a path as harmonic sections of the self-motion bundle

*Companion note to "Pose Path Following using Normal Tubular Coordinates" (Dahlke, Nielsen). Notation follows the
paper: twists $\zeta=(\omega,v)$, body Jacobian $J_b$, left-invariant metric $W$ on $SE(3)$, unit-speed path
$\gamma:\mathbb R/L\mathbb Z\to SE(3)$ with body velocity $\zeta_\gamma=(\gamma^{-1}\gamma')^\vee$, $\|\zeta_\gamma\|_W=1$.
To avoid a clash with the paper's $E_s=\ker J_b(q^\star(s))$, the self-motion bundle is called $\mathcal S$ and the
energy $\mathcal E$.*

**Status tags.** Every claim carries one of:

| tag | meaning |
|---|---|
| **[P]** | proven in this note, with a complete proof (in the text or in Appendix A) |
| **[PH]** | proven under the hypotheses written in the statement |
| **[PS]** | standard argument, proof sketched; the omitted details are routine but not written out |
| **[K]** | known result, cited; not re-proven here |
| **[C]** | conjecture or open question |
| **[N]** | numerical observation on the paper's Panda example (Euclidean joint metric), produced by the scripts in this folder |

Scripts: `pandageo.py` (kinematics, path; checked against the design-doc anchors), `lifts.py` (lifts, discrete
harmonic sections). Each experiment script writes the log of the same name:

| script → log | content |
|---|---|
| `run1.py` | rest-posture lift; pseudoinverse holonomy from $q^\star_{\rm rest}(0)$ |
| `run_thmA.py` | Theorem A check |
| `run_holo.py` | self-motion arc over $\gamma(0)$; holonomy map $T$ on 41 starts |
| `run_harm.py 1000` → `harm_1000.log` | harmonic lift, $m=1000$ |
| `run_multi.py` | multi-start from constant offsets (summarised in Section 5) |
| `run2.py mp/k/v` → `run2_*.log` | index-1 point; $k$-continuation; vertical energy |
| `run3.py a/d` → `run3_*.log` | second minimiser; corrector constant; stretch rate; arcs at 7 more $s$ |
| `run4.py` | vertically harmonic lift details (Cor. 4.4 check) |
| `run_shoot.py` | one lap of the acceleration-level ODE (4.2) |

Every number tagged [N] below comes from these logs (Section 10).

---

## 0. Summary

The working hypothesis ("closed lifts are sections of a bundle over the loop; the canonical one is a harmonic
section; the paper's frames $N$ and $X$ are harmonic sections too") holds. It needs three corrections and several
sharpenings.

1. **Bundle and holonomy.** Closed lifts are exactly the sections of $\mathcal S=\gamma^*\mathsf Q_{\rm adm}\to S^1$
   [P]. $\mathcal S$ is a fibre bundle in the compact case (Ehresmann) and, for $r=1$, whenever the admissible
   self-motion pieces are single arcs (Lemma 1.4) [P]. Pseudoinverse drift is the holonomy of the $K$-orthogonal
   Ehresmann connection, and closed pseudoinverse lifts are its fixed points [P]. On the Panda path the holonomy has
   **no fixed point on the whole admissible self-motion arc** over $\gamma(0)$: $T(t)-t\in[-0.65,-0.22]$ rad of joint
   arclength. So no velocity-level pseudoinverse lift closes on this component [N].
2. **Theorem A (frames = Coulomb gauge).** The closed $W$-orthonormal normal frames that minimise
   $\tfrac12\oint|N^\top WN'|_F^2$ are exactly $N_0e^{-(s/L)\Phi}Q_0$, with $\Phi$ a principal logarithm of the
   holonomy and $Q_0\in O(5)$ arbitrary. The minimum is $\sum\theta_i^2/L$. The minimiser is unique modulo $Q_0$ iff
   $-1\notin\operatorname{spec}H$. Critical points correspond to *all* logarithms; they are the Coulomb gauges (constant
   connection form) [P]. The same holds for Lemma 3's kernel frame [P].
3. **Theorem B (harmonic lift).** The Euler–Lagrange equation is $(\nabla_{q'}q')^V=0$ [P]. For $K=M$ it reads
   $X^\top(Mq''+Cq')=0$. Equivalently, the gravity-compensated torque that executes the lift at any constant speed has
   **zero dynamically consistent null-space component**, and the lift is a periodic motion of the passive,
   gravity-free arm dragged by its flange [P]. For $K=I$ it reads $q''=J_b^+(\zeta_\gamma'-\dot J_bq')$: **the harmonic
   lift is a periodic orbit of the acceleration-level pseudoinverse** [P]. Numerically, one lap of that ODE from the
   lift's initial data closes to $7\cdot10^{-6}$ [N].
4. **Correction 1 (two energies).** The full energy $\mathcal E$ and Wood's vertical energy $\mathcal E^v$ differ by
   $\tfrac12\oint\zeta_\gamma^\top\Lambda_K(q)\zeta_\gamma$, the operational-space inertia along the path. They agree
   up to a constant in the principal case of Theorem A, but not for the self-motion bundle [P]. "Horizontal lift plus
   uniform vertical drift" (Prop. D) holds for $\mathcal E^v$. For $\mathcal E$ it holds only if $\Lambda_K$ is
   constant along the fibres; as stated in the hypothesis it is false in general [P].
5. **Correction 2 (no canonical winding 0).** For circle fibres the homotopy classes of closed lifts form a
   $\mathbb Z$-*torsor*: the integer depends on a trivialisation, and no class is distinguished [P]. The canonical
   choice is the **global** minimiser over all classes, which also selects the class [P]. With joint limits (arc
   fibres) there is exactly one class [P].
6. **Correction 3 (barrier metrics).** A conformal change $K=\varphi^2K_0$ does **not** change the connection, so
   pseudoinverse transport still runs into the limits [P]. It does make the full energy coercive, because it inflates
   $\Lambda_K$ [P]. The vertical energy stays non-coercive on arc fibres.
7. **Theorem C ($r=1$).** Existence [PH], Aubry's crossing lemma and total ordering of minimisers [P], canonical
   extremal minimisers on arc fibres [P], uniqueness and nondegeneracy open and dense in potentials [P]. A
   nondegenerate local minimiser is a **hyperbolic** periodic orbit with positive multipliers; more generally
   $\operatorname{sign}(\operatorname{tr}DP-2)=(-1)^{\text{Morse index}}$ [P, from Hill's oscillation theorem K].
   Minimising harmonic lifts are therefore unstable orbits of the passive dynamics. They cannot be found by forward
   simulation; a boundary-value problem is needed.
8. **Prop. E.** The rest-posture lift is the $k\to\infty$ limit of minimisers of $\mathcal E+\tfrac k2\oint|q-q_r|^2$
   [PH]. The asymptotic rate is $O(1/k)$ [PS]. On the Panda the asymptotic regime starts only at $k\gtrsim10^3$, with
   $k\cdot\mathrm{dist}$ approaching the predicted corrector constant 18.9 (18.2 at $k=10^4$) after a pre-asymptotic regime [N].
9. **Panda, $K=I$ [N].** The minimising harmonic lift exists inside the joint limits. It lies $0.64$–$0.87$ rad from the
   rest-posture lift along the self-motion, with 9.1 % lower energy (6.6696 vs 7.3398). Its joint-4 margin is the same
   (0.3596 vs 0.3592 rad: this margin is set by the path); its smallest $\sigma_{\min}(J_b)$ is lower (0.159 vs 0.207).
   It is nondegenerate (lowest Hill eigenvalue 2.94) and hyperbolic (multipliers 3.81, 0.262). Newton from the
   rest-posture lift converges in 4 iterations. A second, higher local minimiser ($\mathcal E=6.783$) and an index-1,
   elliptic harmonic lift between the two ($\mathcal E=7.460$, trace $-1.33$) also exist, exactly as the
   ordering/minimax picture predicts.

**Recommendation.** Define the paper's lift as the **minimising harmonic section of the full energy** (with $K=M$ when
a dynamic model is available, otherwise $K=I$, and a barrier term for guaranteed existence). Present Lemma 1, Lemma 3
and the lift as three instances of one principle: harmonic sections of natural bundles over the loop. The first two
are the principal ("gauge") case and have closed forms; the third does not.

---

## 1. Setting

### 1.1 Spaces and the vertical bundle

$\mathsf Q$ is the joint space: $\mathbb R^n$, or $T^n=(\mathbb R/2\pi\mathbb Z)^n$ for unlimited revolute joints; in
both cases $\mathsf Q$ is oriented. $n=6+r$, $r\ge1$, and $\Pi:\mathsf Q\to SE(3)$ is the (real-analytic) forward
kinematics with body Jacobian $J_b(q)\in\mathbb R^{6\times n}$. Let

$$
\mathsf Q_{\rm reg}=\{q:\operatorname{rank}J_b(q)=6\},\qquad \Sigma=\mathsf Q\setminus\mathsf Q_{\rm reg},\qquad
\mathsf Q_{\rm adm}=\mathsf Q_{\rm reg}\cap\operatorname{int}B ,
$$

with $B=\prod_i[q_i^-,q_i^+]$ the joint box ($B=\mathsf Q$ without limits). $\Pi|_{\mathsf Q_{\rm reg}}$ is a
submersion, and the **vertical bundle** $V=\ker d\Pi\subset T\mathsf Q_{\rm reg}$ has rank $r$; in coordinates
$V_q=\ker J_b(q)$. Its leaves (the fibres $\Pi^{-1}(g)\cap\mathsf Q_{\rm reg}$) are the self-motion manifolds of Burdick
(1989). Since $T\mathsf Q/V\cong\Pi^*TSE(3)$ is trivial and $\mathsf Q$ is oriented, $V$ is oriented [P]. For $r=1$ it
is trivialised by the unit generalised cross product $x(q)$ of the rows of $J_b(q)$ (Euclidean), as in the paper.

### 1.2 The metric and the Ehresmann connection

Let $K$ be a Riemannian metric on $\mathsf Q_{\rm adm}$, written as a matrix field $K(q)$. Examples are the Euclidean
metric $K=I$, the kinetic-energy metric $K=M(q)$, and the barrier metrics of Section 8. The **horizontal bundle**
$H=V^{\perp_K}$ is an Ehresmann connection on $\Pi|_{\mathsf Q_{\rm adm}}$. Write $\mathcal V,\mathcal H$ for the
$K$-orthogonal projectors onto $V,H$. In coordinates,

$$
H_q=\operatorname{range}K^{-1}J_b^\top,\qquad
\Lambda_K:=(J_bK^{-1}J_b^\top)^{-1},\qquad
J_K^{\#}:=K^{-1}J_b^\top\Lambda_K,\qquad \mathcal H=J_K^\#J_b,\quad \mathcal V=I-J_K^\#J_b .
$$

The horizontal lift of a body velocity $\zeta$ at $q$ is $J^\#_K(q)\zeta$. For $K=I$ this is the Moore–Penrose
pseudoinverse $J_b^+$. For $K=M$ it is Khatib's **dynamically consistent generalised inverse**
$\bar J=M^{-1}J_b^\top(J_bM^{-1}J_b^\top)^{-1}$, and $\Lambda_M$ is the operational-space inertia (Khatib 1987). So the
"dynamically consistent connection" is $H=V^{\perp_M}$ [P, by definition]. Two identities are used repeatedly [P]:

$$
|J^\#_K\zeta|_K^2=\zeta^\top\Lambda_K\zeta,\qquad |w|_K^2=(J_bw)^\top\Lambda_K(J_bw)+|\mathcal Vw|_K^2 .
\tag{1.1}
$$

*Proof.* $(J^\#\zeta)^\top KJ^\#\zeta=\zeta^\top\Lambda J_bK^{-1}KK^{-1}J_b^\top\Lambda\zeta=\zeta^\top\Lambda\zeta$.
For the second identity, $w=J^\#J_bw+\mathcal Vw$ is a $K$-orthogonal splitting. $\square$

**Conformal invariance [P].** If $\tilde K=\varphi^2K$ with $\varphi>0$, then $\tilde H=H$, $J^\#_{\tilde K}=J^\#_K$
and $\Lambda_{\tilde K}=\varphi^2\Lambda_K$. Orthogonal complements do not change under conformal rescaling. This is
used in Section 8.

**Curvature.** For vector fields $\xi,\eta$ on $SE(3)$ with horizontal lifts $\xi^h,\eta^h$, set
$F(\xi,\eta):=-\mathcal V[\xi^h,\eta^h]$. $F$ is a $V$-valued tensor, the curvature of $H$ (sign conventions vary;
Kolář–Michor–Slovák 1993). Frobenius: $H$ is integrable iff $F\equiv0$ [K].

### 1.3 The self-motion bundle; closed lifts are sections

$$
\mathcal S:=\{(s,q)\in\mathbb R/L\mathbb Z\times\mathsf Q_{\rm adm}:\ \Pi(q)=\gamma(s)\},\qquad p(s,q)=s,\qquad
\mathcal S_s:=p^{-1}(s)\cong\Pi^{-1}(\gamma(s))\cap\mathsf Q_{\rm adm}.
$$

**Lemma 1.1 [P].** $\mathcal S$ is a smooth $(1+r)$-manifold and $p$ is a submersion. A smooth $L$-periodic map
$q^\star$ satisfies Assumption 1 of the paper (inside $\mathsf Q_{\rm adm}$) iff $s\mapsto(s,q^\star(s))$ is a smooth
section of $p$. The map $(s,q)\mapsto q$ embeds $\mathcal S$ onto $\mathsf Z:=\Pi^{-1}(\mathcal C)\cap\mathsf Q_{\rm adm}$
(the "self-motion tube" over the path $\mathcal C=\gamma(\mathbb R/L\mathbb Z)$).

*Proof.* $(s,q)\mapsto\gamma(s)^{-1}\Pi(q)$ is a submersion at every point, because $d\Pi$ is onto. So
$\mathcal S=\{\gamma(s)^{-1}\Pi(q)=I\}$ is a submanifold of dimension $1+n-6$. Given $\dot s$, the equation
$J_b\dot q=\zeta_\gamma\dot s$ has solutions, so $dp$ is onto. $\gamma$ is injective (embedded circle), so
$(s,q)\mapsto q$ is injective. It is an immersion: a tangent vector $(\dot s,\dot q)$ with $\dot q=0$ has
$\zeta_\gamma\dot s=J_b\dot q=0$, so $\dot s=0$. It is
proper onto its image $\mathsf Z$ (closed in $\mathsf Q_{\rm adm}$), hence an embedding. $\square$

The connection induces on $\mathcal S$ the vector field
$\tilde Y(s,q)=(\partial_s,\ Y_s(q))$, $Y_s(q):=J^\#_K(q)\zeta_\gamma(s)$ (the horizontal lift of $\partial_s$). Its
local flow $\Theta$ is **horizontal transport**: $\Theta_s$ maps (part of) $\mathcal S_0$ into $\mathcal S_s$. For
$K=I$ it is velocity-level pseudoinverse redundancy resolution, and $T:=\Theta_L|_{\mathcal S_0}$ is the **holonomy
(monodromy)**.

### 1.4 When is $\mathcal S$ a fibre bundle?

**Proposition 1.2 (compact case) [P].** Suppose $\mathsf Q=T^n$ without joint limits, and let $\mathcal S'$ be a union
of connected components of $\gamma^*\mathsf Q_{\rm reg}$ that is closed in $\mathbb R/L\mathbb Z\times T^n$ (it does
not accumulate on $\Sigma$). For instance, take all of $\gamma^*\mathsf Q_{\rm reg}$ when $\mathcal C$ avoids the
critical values $\Pi(\Sigma)$. Then $p:\mathcal S'\to S^1$ is a smooth fibre bundle with compact fibre. For $r=1$ the
fibre is a finite disjoint union of circles.

*Proof.* $\mathcal S'$ is compact, so $p$ is a proper surjective submersion. Ehresmann's fibration theorem (Ehresmann
1951) [K] applies. $\square$

**Proposition 1.3 (complete transport ⇒ mapping torus) [P].** If the flow of some lift $\tilde Y$ of $\partial_s$ to
$\mathcal S$ is complete (for instance the $K$-horizontal one), then $\Theta_t:\mathcal S_s\to\mathcal S_{s+t}$ are
diffeomorphisms, and $\mathcal S$ is diffeomorphic to the **mapping torus**
$\mathcal M_T=[0,L]\times F/\{(L,u)\sim(0,T(u))\}$, with $F=\mathcal S_0$ and $T=\Theta_L|_F$.

*Proof.* $(t,u)\mapsto\Theta_t(u)$ from $\mathbb R\times F$ is a surjective local diffeomorphism. Its only
identifications are $(t+L,u)\sim(t,T(u))$, because $\Theta_{t+L}=\Theta_t\circ\Theta_L$ on $F$. $\square$

When $\mathcal S$ is compact, every lift of $\partial_s$ is complete. With joint limits, the $K$-horizontal lift is
usually **not** complete. Pseudoinverse transport reaches a joint limit in finite $s$, and the Panda example does this
after a few laps (Section 2.4). For $r=1$, the bundle structure can be obtained without any connection.

**Lemma 1.4 (arc fibres, $r=1$) [P].** Let $\mathcal S^\circ\subset\mathcal S$ be open, such that every fibre
$\mathcal S^\circ_s$ is connected and non-compact (an open arc). Then $p:\mathcal S^\circ\to S^1$ is a smooth fibre
bundle with fibre $\mathbb R$. Since $V$ is oriented, it is trivial, $\mathcal S^\circ\cong S^1\times\mathbb R$.

*Proof (Appendix A.1).* Take a complete metric $g$ on $\mathcal S^\circ$ (Nomizu–Ozeki 1961) and the $g$-unit oriented
vertical field $\hat\nu$. Its flow is complete, because fibres are closed in $\mathcal S^\circ$. Through a local
section $\sigma$, $(s,t)\mapsto\Phi^{\hat\nu}_t(\sigma(s))$ is a fibre-preserving diffeomorphism
$U\times\mathbb R\to p^{-1}(U)$. $\square$

A natural $\mathcal S^\circ$ is the union, over $s$, of the components of $\mathcal S_s$ that contain $q_{\rm ref}(s)$
for a given continuous closed lift $q_{\rm ref}$, such as the rest-posture lift. This set is open [P]: a compact arc
inside the open set $\mathcal S$ persists under small changes of $s$. Lemma 1.4 applies if none of these components is
a whole circle. On the Panda path this was checked at eight values of $s$: the component through the rest-posture lift
is a single arc (joint arclength 5.5–8.6 rad), ended at both sides by the limits of joint 2 [N, Section 10].

---

## 2. Topology: existence and classification of closed lifts; why pseudoinverse control drifts

### 2.1 Sections of a mapping torus

Let $\mathcal S\cong\mathcal M_T$ (Prop. 1.3). In the transport picture a section is $s\mapsto\Theta_s(c(s))$ with
$c:[0,L]\to F$ and

$$
c(0)=T(c(L))\qquad(\text{and } c'(0)=dT\,c'(L)\text{ for } C^1\text{ closure}). \tag{2.1}
$$

**Proposition 2.1 [P].** (a) A closed lift through the component $F_i$ of $F$ exists iff $T(F_i)=F_i$. (b) If $T$
permutes the components cyclically with period $m\ge2$, $T^m(F_i)=F_i$, there are no closed lifts through $F_i$. There
are lifts that close after $m$ laps (sections of the $m$-fold cover), so the arm returns to its configuration every
$m$ laps. (c) If $F$ is connected, closed lifts exist.

*Proof.* A continuous $c$ stays in one component $F_i$, and (2.1) requires $T(c(L))\in F_i$, i.e. $T(F_i)=F_i$. If
$F_i$ is path connected and $T(F_i)=F_i$, any path from $T(y)$ to $y$ in $F_i$ satisfies (2.1). It can be chosen
smooth and with matching derivatives at the ends. For (b), apply (a) to $T^m$ on the $m$-fold cover. $\square$

### 2.2 Classification of closed lifts up to homotopy

Two closed lifts are *equivalent* if they are homotopic through closed lifts (sections).

**Proposition 2.2.** (a) [K, standard; not used below] For connected $F$, equivalence classes of sections of
$\mathcal M_T$ correspond bijectively to the Reidemeister (twisted-conjugacy) classes of $T_\#$ on $\pi_1(F)$:
$\alpha\sim\beta\,\alpha\,T_\#(\beta)^{-1}$ (with $T_\#$ defined by a base path from $x_0$ to $T(x_0)$). Via
$\pi_1(\mathcal M_T)\cong\pi_1(F)\rtimes_{T_\#}\mathbb Z$, sections are degree-one loops, and their classes are the
twisted-conjugacy classes; this is the mapping-torus picture of Nielsen fixed-point theory (e.g. Jiang 1983). Only the
explicit cases (b)–(d), which are proven here, are used below. (b) [P] $F=S^1$ and $T$ orientation preserving (always so here: $V$ is
oriented and transport preserves the orientation): the classes form a $\mathbb Z$-torsor. A trivialisation
$\mathcal S\cong S^1\times S^1$ identifies them with the winding number $k$ of the section around the fibre per lap.
Changing the trivialisation shifts every $k$ by the same integer, so **no class is distinguished**. (c) [P]
$F=S^1$, $T$ orientation reversing (Klein bottle; does not occur for oriented $\mathsf Q$): two classes. (d) [P]
$F\cong\mathbb R$: a single class.

*Proof of (b)–(d).* For (b), $T$ lifts to $\tilde T:\mathbb R\to\mathbb R$ with $\tilde T(u+\ell)=\tilde T(u)+\ell$.
Sections correspond to $\tilde c:[0,L]\to\mathbb R$ with $\tilde c(0)=\tilde T(\tilde c(L))+k\ell$, and $k$ is
invariant under homotopy. For fixed $k$ the space of such $\tilde c$ is contractible: it fibres over $\tilde c(L)\in\mathbb R$ with
contractible path-space fibres. Replacing $\tilde T$ by $\tilde T+j\ell$ (another lift) replaces $k$ by $k-j$. For
(c), $\tilde T(u+\ell)=\tilde T(u)-\ell$, and $\tilde c\mapsto\tilde c+\ell$ changes $k$ by 2, so $k$ mod 2 classifies.
For (d), the fibre is contractible, so in $\mathcal S^\circ\cong S^1\times\mathbb R$ sections are functions
$S^1\to\mathbb R$, joined by convex combination. $\square$

**Consequence (Correction 2).** "The winding-0 minimiser" is meaningful only after a reference section has been
fixed, for example the rest-posture lift. Without one, the canonical object is the global minimiser over all classes
(Theorem C7). With joint limits ($F$ an arc) the question disappears: there is one class.

### 2.3 Holonomy, curvature and repeatability

**Proposition 2.3 [P].** Velocity-level resolution $q'=J^\#_K(q)\zeta_\gamma(s)$ (pseudoinverse control for $K=I$)
started at $q_0\in F$ returns to $q_0$ after one lap iff $T(q_0)=q_0$. Closed horizontal lifts are the fixed points
of the holonomy.

This is a tautology given Section 1.3, but it places the classical results in the right setting.

- **Klein and Huang (1983)** observed that pseudoinverse control is not cyclic. In this language, $T\neq\mathrm{id}$ in
  general.
- **Shamir and Yomdin (1988)** showed that a resolution of the form $\dot q=G(q)\dot x$ is repeatable for all closed
  task paths near a configuration iff the distribution spanned by the columns of $G$ is integrable (Frobenius). For
  $G=J^\#_K$ this says: *repeatable for all small loops* $\iff F\equiv0$, i.e. **flat connection** [K]. The
  infinitesimal reason is that the holonomy of a small loop spanning the plane of $\xi,\eta$ with area $\varepsilon^2$
  displaces by $\varepsilon^2F(\xi,\eta)+o(\varepsilon^2)$ (sign depending on convention) [PS].
- **Baker and Wampler (1988)** showed that a cyclic resolution on a region is the same as an inverse-kinematic
  function (a section of $\Pi$ over the region), with topological obstructions to its existence [K]. In connection
  terms, a flat connection with trivial holonomy has integral manifolds that are sections.

The point for the paper is different. **Closure along one given loop does not require flatness**, because a section
of $\mathcal S$ over the circle always exists (Prop. 2.1). What fails is only that a *Markovian*, velocity-level rule
produces it. Harmonic sections are non-local in $s$: they are determined by the whole loop.

### 2.4 The Panda path [N]

With $K=I$, one lap of $q'=J_b^+(\zeta_\gamma+50e)$ ($e$ is a small pose-error feedback that keeps $q$ on the fibre)
from $q^\star_{\rm rest}(0)$ ends $0.2423$ rad (2-norm) away, almost entirely along the self-motion:
$x^\top(q(L)-q(0))=-0.2423$. This reproduces the lead's 0.24 rad. The admissible component of $\mathcal S_0$ through
$q^\star_{\rm rest}(0)$ is an arc of joint arclength $t\in[-4.28,4.29]$, ended at both sides by joint 2's limits, with
$\sigma_{\min}(J_b)\ge0.103$ on it. Of 41 starting points spread along the arc, 35 stay in the joint box for the whole lap, and for them
$T(t)-t\in[-0.65,-0.22]$; between neighbouring samples it changes by at most $0.094$. The other 6 (near the two ends)
leave the box during the lap, so they are not admissible. **The holonomy has no fixed point on the admissible arc.** Every start drifts monotonically
towards joint 2's lower-$t$ limit, and no velocity-level pseudoinverse lift closes on this component.

---

## 3. Theorem A: the paper's frames are harmonic sections (Coulomb gauge on the circle)

### 3.1 The normal connection and the paper's transport

Let $\nu_s:=\zeta_\gamma(s)^{\perp_W}\subset\mathbb R^6$. This is the normal bundle of the path in body coordinates, of
rank 5, with fibre metric $W$. With the $W$-orthogonal projector $P_\nu=I-\zeta_\gamma\zeta_\gamma^\top W$, define the
**normal connection** $\nabla^\nu_sn:=P_\nu n'$ for sections $n$ of $\nu$. This is the connection induced on the
subbundle $\nu$ by the flat connection $d/ds$ of the left trivialisation. Differentiating $\zeta_\gamma^\top Wn=0$
gives $\zeta_\gamma^\top Wn'=-\zeta_\gamma'^\top Wn$, hence

$$
\nabla^\nu_s n=n'+\zeta_\gamma\,\zeta_\gamma'^\top Wn .\tag{3.1}
$$

**Lemma 3.1 [P].** The solution of the paper's (eq:transport), $N_0'=-\zeta_\gamma\zeta_\gamma'^\top WN_0$, with
$N_0(0)^\top WN_0(0)=I$ and $N_0(0)^\top W\zeta_\gamma(0)=0$, satisfies for all $s$:
$\zeta_\gamma^\top WN_0=0$, $N_0^\top WN_0=I$, $\nabla^\nu N_0=0$, and

$$
N_0^\top W N_0'=0 .\tag{3.2}
$$

$N_0$ is the **rotation-minimising (Bishop) frame** of $\nu$, i.e. parallel for $\nabla^\nu$. Its holonomy is
$H=N_0(0)^\top WN_0(L)\in SO(5)$, and $N_0(s+L)=N_0(s)H$.

*Proof.* $(\zeta_\gamma^\top WN_0)'=\zeta_\gamma'^\top WN_0-(\zeta_\gamma^\top W\zeta_\gamma)\zeta_\gamma'^\top WN_0
=0$, using $\|\zeta_\gamma\|_W=1$. Then $(N_0^\top WN_0)'=-N_0^\top W\zeta_\gamma'\zeta_\gamma^\top WN_0-N_0^\top
W\zeta_\gamma\zeta_\gamma'^\top WN_0=0$. Equation (3.2) follows from $N_0^\top W\zeta_\gamma=0$, and $\nabla^\nu N_0=0$
from (3.1). $\det H=1$ because $\nu\oplus\operatorname{span}\zeta_\gamma=\mathbb R^6$, so $\nu$ is oriented and
transport preserves orientation. $N_0(s+L)$ and $N_0(s)H$ solve the same $L$-periodic linear ODE with the same value at
$s=0$. $\square$

### 3.2 Closed frames, connection form, energy

A **closed frame** is an $L$-periodic $H^1$ map $N$ with $N^\top WN=I$ and $N^\top W\zeta_\gamma=0$. Every closed
frame is $N=N_0Q$ with $Q:[0,L]\to O(5)$ in $H^1$ satisfying

$$
Q(L)=H^{-1}Q(0). \tag{3.3}
$$

Its **connection form** is $A_N:=N^\top WN'\in\mathfrak{so}(5)$. By (3.2), $A_N=Q^\top Q'$, so $|A_N|_F=|Q'|_F$. Define

$$
\mathcal E^v[N]:=\tfrac12\oint|N^\top WN'|_F^2\,ds=\tfrac12\int_0^L|Q'|_F^2\,ds .
$$

This is the energy of the curve $Q$ in $O(5)$ for the bi-invariant metric $\langle A,B\rangle=\operatorname{tr}(A^\top B)$.

**Theorem A [P].** Let $\Phi=\log H$ be a principal logarithm (rotation angles $\theta_i\in[0,\pi]$ in the invariant
planes of $H$).

1. $\min\mathcal E^v=\dfrac{|\Phi|_F^2}{2L}=\dfrac{1}{L}\sum_i\theta_i^2$. The minimisers are exactly
   $$N(s)=N_0(s)\,e^{-(s/L)\Phi}\,Q_0,\qquad Q_0\in O(5),\ \Phi\text{ a principal logarithm of }H .$$
   These are smooth and $L$-periodic on $\mathbb R$, and include the paper's frame ($Q_0=I$).
2. The minimiser is unique modulo right multiplication by a constant $Q_0\in O(5)$ iff $-1\notin\operatorname{spec}H$.
3. The **critical points** of $\mathcal E^v$ (variations $N\mapsto Ne^{\epsilon\Xi}$, $\Xi$ periodic) are exactly the
   closed frames with **constant connection form**, i.e. the Coulomb gauge $\partial_sA_N=0$. They are
   $N=N_0e^{-(s/L)Y}Q_0$ with $Y$ *any* real logarithm of $H$, and have energy $|Y|^2_F/(2L)$.
4. The full energy satisfies $\tfrac12\oint|N'|_W^2=\tfrac12\oint\|\zeta_\gamma'\|_W^2+\mathcal E^v[N]$, where
   $|B|_W^2:=\operatorname{tr}B^\top WB$. The full and vertical energies therefore have the same critical points and
   minimisers.

*Proof.* (1) For $Q$ satisfying (3.3), Cauchy–Schwarz gives
$\int_0^L|Q'|^2\ge\frac1L(\int_0^L|Q'|)^2\ge\frac1Ld(Q(0),H^{-1}Q(0))^2$, with $d$ the bi-invariant distance. By
bi-invariance $d(Q_0,H^{-1}Q_0)=d(I,Q_0^\top H^{-1}Q_0)$. In $SO(m)$ with this metric, $d(I,R)^2=2\sum\theta_i^2(R)$
with angles in $[0,\pi]$ (Appendix A.2) [P]. Conjugation preserves the angles, and $H^{-1}$ has the same angles as $H$,
so the bound is $|\Phi|_F^2/(2L)$ for every $Q_0$. Equality holds iff $Q$ is a constant-speed minimising geodesic from
$Q_0$ to $H^{-1}Q_0$. Those geodesics are $Q_0e^{(s/L)Z}$ with $e^Z=Q_0^\top H^{-1}Q_0$ and $|Z|_F$ minimal, i.e.
$Z=-Q_0^\top\Phi Q_0$ for a principal $\Phi$. Then $Q(s)=e^{-(s/L)\Phi}Q_0$. Periodicity: $N(s+L)=N_0(s)He^{-(s/L)\Phi}
e^{-\Phi}Q_0=N(s)$, since $H=e^\Phi$ commutes with $e^{-(s/L)\Phi}$.

(2) The principal logarithm is unique iff no angle equals $\pi$. An angle $\pi$ gives several (or a continuum of)
minimal logarithms, i.e. inequivalent minimiser families. Their ratios $e^{(s/L)Y_1}e^{-(s/L)Y_2}$ are not constant.

(3) Under $N\mapsto Ne^{\epsilon\Xi}$, $\delta A=\Xi'+[A,\Xi]$, so
$\delta\mathcal E^v=\oint\langle A,\Xi'+[A,\Xi]\rangle=\oint\langle A,\Xi'\rangle=-\oint\langle A',\Xi\rangle$, using
$\langle A,[A,\Xi]\rangle=0$ (ad-invariance). Critical iff $A'\equiv0$. If $A$ is constant, then $Q=Q_0e^{sA}$ and
(3.3) give $e^{LA}=Q_0^\top H^{-1}Q_0$, so $LA=-Q_0^\top YQ_0$ with $e^Y=H$.

(4) $N'=N_0'Q+N_0Q'$. $N_0'Q$ has range in $\operatorname{span}\zeta_\gamma$ and $N_0Q'$ has range in $\nu$, so the
cross term vanishes. Moreover $|N_0'Q|_W^2=\operatorname{tr}(N_0^\top W\zeta_\gamma'\zeta_\gamma'^\top WN_0)=\|P_\nu
\zeta_\gamma'\|_W^2=\|\zeta_\gamma'\|_W^2$, since $\zeta_\gamma^\top W\zeta_\gamma'=0$. $\square$

**Interpretation.** On the frame bundle $O(\nu)\to S^1$, with $\nabla^\nu$ as connection, a closed frame is a gauge.
$\mathcal E^v$ is the $L^2$ norm of the connection form in that gauge. On $S^1$ the Coulomb condition $d^*A=0$ means
$A$ constant, and the gauge-equivalence classes of connections are conjugacy classes of holonomy (standard gauge theory
on the circle). Theorem A says that **the paper's holonomy cancellation is the minimal-energy Coulomb gauge**. For
space curves ($SO(2)$, a single angle $\theta$) this is the familiar "rotation-minimising frame plus uniformly
distributed closing twist $\theta/L$". It appears for closed curves in computer graphics (the double-reflection RMF of
Wang, Jüttler, Zheng and Liu 2008 is corrected this way) and as the zero-material-twist equilibrium of closed elastic
rods (Bergou, Wardetzky, Robinson, Audoly and Grinspun 2008). The rotation-minimising frame is Bishop's (1975).

**Numerical check [N].** On the Panda path ($\ell=0.1$ m, 2000 RK4 steps) the angles are $0.161182$ and
$8.0075\cdot10^{-5}$, as in the design doc. $|\Phi|_F^2/(2L)=0.0337361403=\sum\theta_i^2/L$. Quadrature of
$\mathcal E^v$ for the paper's frame gives $0.0337361403$. Three closed perturbations
$e^{-(s/L)\Phi}e^{0.05\sin(2\pi ks/L)\Xi}$ give $0.67$, $4.46$ and $4.47$. Identity (4) holds to $5\cdot10^{-12}$
relative. $\max|N_0^\top WN_0'|=4.6\cdot10^{-11}$.

### 3.3 Lemma 3 (kernel frame) is the same theorem

Along a closed lift $q^\star$, $E_s=\ker J_b(q^\star(s))=V_{q^\star(s)}$ is a rank-$r$ subbundle of $q^{\star*}T\mathsf Q$.
It has the fibre metric $K$ and the induced connection $\nabla^E_sx=\mathcal V D_sx$, where $D_s$ is the $K$-Levi-Civita
derivative along $q^\star$.

**Lemma 3.2 [P].** The paper's transport $D_sx^0=-(D_sP)x^0$ ($P=I-\mathcal V$ along $q^\star$) is exactly
$\nabla^E$-parallel transport.

*Proof.* From $P^2=P$, $(I-P)(D_sP)(I-P)=0$. With $y:=Px^0$, $D_sy=(D_sP)x^0+PD_sx^0=(I-P)(D_sP)x^0=(I-P)(D_sP)y$, a
linear ODE with $y(0)=0$, so $y\equiv0$ and $x^0\in E$. Then $\mathcal VD_sx^0=-(I-P)(D_sP)(I-P)x^0=0$. $\square$

Theorem A therefore holds verbatim for $X$, with $O(r)$, the holonomy $H\in SO(r)$ ($\det H=1$ by orientability of
$V$, i.e. the paper's Assumption 2) and $\mathcal E^v[X]=\frac12\oint|X^\top KD_sX|_F^2$ [P]. For $r=1$, $SO(1)$ is
trivial, $\Phi=0$, and every unit section of the oriented line bundle is parallel:
$\nabla^Ex=K(D_sx,x)x=\tfrac12(|x|^2)'x=0$. This is why the generalised cross product needs no integration [P].

---

## 4. Theorem B: harmonic lifts

### 4.1 Two energies

Let $\Gamma^1(\mathcal S)$ be the $H^1$ sections, i.e. $L$-periodic $H^1$ maps $q$ into $\mathsf Q_{\rm adm}$ with
$\Pi(q(s))=\gamma(s)$. Its tangent space at $q$ consists of the $H^1$ vertical fields $w$ along $q$, $w(s)\in
V_{q(s)}$. Every section satisfies $J_b(q)q'=\zeta_\gamma(s)$. Define

$$
\mathcal E_K(q)=\tfrac12\oint|q'|_K^2\,ds,\qquad
\mathcal E^v_K(q)=\tfrac12\oint|\mathcal Vq'|_K^2\,ds,\qquad
U(s,q):=\zeta_\gamma(s)^\top\Lambda_K(q)\zeta_\gamma(s).
$$

By (1.1), with $J_bq'=\zeta_\gamma$ [P]:

$$
\boxed{\ \mathcal E_K(q)=\tfrac12\oint U(s,q(s))\,ds+\mathcal E^v_K(q)\ }\tag{4.1}
$$

**Definition.** A *harmonic lift* is a critical point of $\mathcal E_K$ on $\Gamma^1(\mathcal S)$. A *minimising
harmonic lift* minimises $\mathcal E_K$ in its homotopy class of sections. The *canonical lift* is the global
minimiser, when it is unique. A *vertically harmonic lift* is a critical point of $\mathcal E^v_K$. This is the notion
of harmonic section of C. M. Wood (1987), defined there for Riemannian submersions.

**Lemma 4.1 (regularity) [PS].** An $H^1$ critical point of $\mathcal E_K$ or $\mathcal E^v_K$ is $C^\infty$, and
real-analytic if $\Pi,\gamma,K$ are. *Sketch.* In a local trivialisation the weak Euler–Lagrange equation reads
$\int(L_p\phi'+L_u\phi)=0$. By du Bois-Reymond, $L_p(s,u,u')$ is absolutely continuous with derivative $L_u$. Strict
convexity in $p$ gives $u'$ continuous, and then bootstrapping gives smoothness (Buttazzo, Giaquinta and Hildebrandt
1998, ch. 4). The same argument works for vector $u$ ($r>1$).

### 4.2 First variation and Euler–Lagrange equation

Let $\Gamma_K(q)(v,v):=K^{-1}\big[(dK(q)[v])v-\tfrac12\nabla_q(v^\top K(q)v)\big]$ be the Christoffel term, so that
the $K$-covariant derivative along a curve is $D_sq'=q''+\Gamma_K(q)(q',q')$. For $K=M$ one has
$M\Gamma_M(q)(v,v)=C(q,v)v$, the Coriolis/centrifugal vector; it is independent of the factorisation chosen for $C$.
Write $\dot J_b:=dJ_b(q)[q']$, and let $X(s)$ be any $n\times r$ matrix whose columns span $\ker J_b(q(s))$.

**Theorem B [P].** For a $C^2$ section $q$ the following are equivalent.

- (a) $q$ is a harmonic lift (critical point of $\mathcal E_K$).
- (b) $\mathcal V\,D_sq'=0$: the covariant acceleration has no vertical (self-motion) component.
- (c) $X^\top\big(K(q)q''+K\Gamma_K(q)(q',q')\big)=0$. For $K=M$ this is $X^\top\big(Mq''+C(q,q')q'\big)=0$. For
  $K=I$ it is $X^\top q''=0$, i.e. $q''\in\operatorname{range}J_b^\top$.
- (d) $Kq''+K\Gamma_K(q',q')=J_b(q)^\top\lambda$ for some continuous $\lambda:S^1\to\mathbb R^6$. Then
  $\lambda=\Lambda_K(\zeta_\gamma'-\dot J_bq'+J_b\Gamma_K(q',q'))$.
- (e) $q$ solves the second-order ODE
  $$q''=J^\#_K(q)\big(\zeta_\gamma'(s)-\dot J_bq'\big)-\big(I-J_K^\#J_b\big)\Gamma_K(q)(q',q').\tag{4.2}$$
  For $K=I$ this is $q''=J_b^+(\zeta_\gamma'-\dot J_bq')$, the **acceleration-level pseudoinverse
  (minimum-acceleration-norm) resolution**. For $K=M$ it is the **dynamically consistent acceleration-level
  resolution with zero null-space torque**.
- (f) ($K=M$) For one, equivalently every, $\nu>0$, $q_\nu(t):=q(\nu t)$ is a motion of the **gravity-free arm under
  the ideal rheonomic constraint** $\Pi(q(t))=\gamma(\nu t)$:
  $$M(q_\nu)\ddot q_\nu+C(q_\nu,\dot q_\nu)\dot q_\nu=J_b(q_\nu)^\top f(t),\qquad f(t)=\nu^2\lambda(\nu t).$$
  In words, the flange is dragged along $\mathcal C$ at constant speed by an external wrench $f$, and the arm follows
  passively.
- (g) ($K=M$) The torque that executes $q_\nu$ on the real arm, $\tau=M\ddot q_\nu+C\dot q_\nu+G$, satisfies
  $\big(I-J_b^\top\bar J^\top\big)(\tau-G)=0$. Beyond gravity compensation, it has **no component in Khatib's dynamically
  consistent null space**.

*Proof.* (a)⇔(b): variations through sections have vertical variation fields, and every smooth vertical $w$ comes
from one (retract $q+\epsilon w$ onto the fibres). Then
$\delta\mathcal E_K=\oint K(q',D_sw)=-\oint K(D_sq',w)$, and $w$ ranges over all vertical fields.
(b)⇔(c): $K(D_sq',w)=w^\top(Kq''+K\Gamma_K)$, with $w=Xa$.
(c)⇔(d): $(\operatorname{range}X)^{\perp}=\operatorname{range}J_b^\top$ in the Euclidean pairing.
(d)⇔(e): substitute $q''=K^{-1}J_b^\top\lambda-\Gamma_K$ into the differentiated constraint
$J_bq''+\dot J_bq'=\zeta_\gamma'$, which gives $\lambda$; then $q''=J^\#_K(\zeta'_\gamma-\dot J_bq')+J^\#_KJ_b
\Gamma_K-\Gamma_K$.
(f): $\dot q_\nu=\nu q'$ and $\ddot q_\nu=\nu^2q''$, and $C(q,v)v$ is quadratic in $v$, so the equation is $\nu^2$
times (d). A force $J_b^\top f$ is the generalised force of a body wrench $f$ at the flange, and an ideal constraint
force does no virtual work along $V$ (d'Alembert).
(g): $\tau-G=J_b^\top f\iff(I-J_b^\top\Lambda J_bM^{-1})(\tau-G)=0$. If $\tau-G=J_b^\top f$, the projector annihilates
it; conversely, $\tau-G=J_b^\top(\bar J^\top(\tau-G))$. $\square$

**Remarks.**

1. **Speed independence [P].** By (f), harmonic lifts do not depend on the speed at which the path is traversed.
   Gravity breaks this: the action of the dragged arm with gravity,
   $\int(T-V_G)\,dt=\nu\cdot\tfrac12\oint|q'|_M^2ds-\nu^{-1}\oint V_G\,ds$, depends on $\nu$, and the gravity-free
   harmonic lift is its $\nu\to\infty$ limit [PS]. Gravity-free is the right model, because the controller cancels $G$.
2. **Why forward integration cannot find it.** By (e) the harmonic lift is a periodic orbit of acceleration-level
   resolution. Resolving redundancy at the acceleration level with a pseudoinverse leaves the self-motion undamped and
   can let it grow; acceleration-level schemes add null-space damping for this reason (see e.g. Hsu, Hauser and
   Sastry 1989). Theorem C5 adds that *minimising* periodic orbits are
   hyperbolic, hence unstable. The canonical lift is a saddle of the passive dynamics and must be computed as a
   periodic boundary-value problem (Section 9). Numerically [N], one lap of (4.2) with $K=I$ (plus horizontal
   pose/velocity feedback, which does not affect the vertical part), started from the computed lift's
   $(q^\star(0),q^{\star\prime}(0))$, closes to $|q(L)-q^\star(0)|=7.0\cdot10^{-6}$ and
   $|q'(L)-q^{\star\prime}(0)|=1.5\cdot10^{-5}$. From the rest-posture lift's data it misses by 1.24 rad.
3. **Literature.** The Euler–Lagrange equation of the kinetic-energy integral under a task constraint appears in global
   redundancy resolution (Kazerounian and Wang 1988; Nakamura and Hanafusa 1987 via the maximum principle), there with
   fixed end points. New here are the periodic setting (sections over the loop), the bundle viewpoint, and the
   existence, ordering and hyperbolicity theory of Section 5.

### 4.3 Operational-space reading of the energy

By (4.1), $\mathcal E_K=\tfrac12\oint\zeta_\gamma^\top\Lambda_K\zeta_\gamma+\tfrac12\oint|\mathcal Vq'|_K^2$. The first
term is the operational-space kinetic energy that unit path speed costs in the current posture. For $K=I$ it is
$|J_b^+\zeta_\gamma|^2$, the squared joint speed needed per unit path speed: the inverse square of Chiu's (1988)
velocity transmission ratio along $\zeta_\gamma$. The second term is the self-motion energy. **The harmonic lift trades
task compatibility along the path against self-motion.** It is drawn towards postures in which the path direction is
cheap, and it pays for getting there in self-motion. On the Panda [N]:

| lift | $\tfrac12\oint U$ | $\mathcal E^v$ | $\mathcal E$ |
|---|---|---|---|
| rest posture | 7.2552 | 0.0846 | 7.3398 |
| harmonic (full) | 6.6005 | 0.0691 | 6.6696 |
| vertically harmonic | 7.529 | 0.0328 | 7.5621 |

### 4.4 The vertical energy, O'Neill's tensors, and how the two Euler–Lagrange equations differ

For the splitting $T\mathsf Q=V\oplus H$, O'Neill's (1966) tensor $T$ is $T_vw=\mathcal H\nabla_vw=\mathrm{II}(v,w)$
(the second fundamental form of the self-motion manifolds) and $T_vY=\mathcal V\nabla_vY$ for $v,w\in V$, $Y\in H$.
$T_v$ is $K$-skew. These properties need only the metric and $H=V^\perp$, not a Riemannian submersion.

**Lemma 4.2 [P].** Let $Y$ be *basic* (horizontal and $\Pi$-related to a vector field on $SE(3)$) and $w$ vertical.
Then $K(\nabla_YY,w)=-\tfrac12\,w\big(|Y|_K^2\big)$. For a Riemannian submersion $|Y|^2$ is constant on fibres and the
right side vanishes.

*Proof.* $K(\nabla_YY,w)=-K(Y,\nabla_Yw)=-K(Y,\nabla_wY)-K(Y,[Y,w])$. Since $Y$ is projectable and $w$ vertical,
$[Y,w]$ is vertical, so the last term is 0, and $K(Y,\nabla_wY)=\tfrac12w|Y|^2$. $\square$

Write $q'=v+Y_s(q)$ with $v=\mathcal Vq'$ and $Y_s=J^\#_K\zeta_\gamma(s)$. This is basic for frozen $s$: it is the
horizontal lift of the left-invariant field $g\mapsto g\hat\zeta_\gamma(s)$.

**Proposition 4.3 [P].** The Euler–Lagrange equations are

$$
\text{vertical energy:}\quad \mathcal V D_sv+T_vY_s=0,\qquad\qquad
\text{full energy:}\quad \mathcal V D_sv+T_vY_s=\tfrac12\operatorname{grad}^VU(s,\cdot).
$$

For $r=1$, let $x$ be the oriented $K$-unit vertical field, $\kappa:=\nabla_xx=\mathrm{II}(x,x)$ the curvature vector
of the self-motion curve (horizontal), and $v=K(q',x)$ the signed self-motion speed. Then

$$
\text{vertical:}\quad v'=K(Y_s,\kappa)\,v,\qquad\qquad \text{full:}\quad v'=K(Y_s,\kappa)\,v+\tfrac12\,\partial_xU .\tag{4.3}
$$

*Proof.* Full energy: by Theorem B, $\mathcal VD_sq'=0$, and
$D_sq'=D_sv+(\partial_sY_s)+\nabla_vY_s+\nabla_{Y}Y_s$. Here $\partial_sY_s=J^\#\zeta'_\gamma$ is horizontal,
$\mathcal V\nabla_vY=T_vY$, and $\mathcal V\nabla_YY=-\tfrac12\operatorname{grad}^VU$ by Lemma 4.2, since $|Y_s|^2=U$.

Vertical energy: $\delta\mathcal E^v=\oint K(v,(D_w\mathcal V)q'+\mathcal VD_sw)$. For horizontal $Y$, $\mathcal VY=0$
gives $(D_w\mathcal V)Y=-T_wY$, and $D_w\mathcal V$ maps $V$ into $H$. Hence
$K(v,(D_w\mathcal V)q')=-K(v,T_wY)=K(T_wv,Y)=K(\mathrm{II}(w,v),Y)=-K(w,T_vY)$. Integrating
$\oint K(v,D_sw)$ by parts gives $\delta\mathcal E^v=-\oint K(w,\mathcal VD_sv+T_vY)$.

For $r=1$: $\mathcal VD_s(vx)=v'x$ because $K(D_sx,x)=0$, and $K(T_xY,x)=-K(Y,T_xx)=-K(Y,\kappa)$. $\square$

**Corollary 4.4 [P].** A vertically harmonic lift ($r=1$) has $v(s)=v(0)\exp\int_0^sK(Y_\sigma,\kappa)\,d\sigma$
(evaluated along the lift). It is either a closed horizontal lift ($v\equiv0$) or its self-motion speed never vanishes.
In particular, when the holonomy has no fixed point, as on the Panda, every vertically harmonic lift moves strictly
monotonically along the self-motion.

**Which energy reduces to Theorem A?** In the principal case of Theorem A the "potential" is
$\tfrac12\|\zeta_\gamma'\|_W^2$, which is the same at every point of the fibre (Theorem A(4)): the structure group acts
by isometries that preserve the horizontal lift. So the full and vertical energies differ by a constant and have the
same critical points, and the paper's $N$ is harmonic in both senses [P]. For the self-motion bundle there is no group
action. The two coincide (up to a constant) iff $U(s,\cdot)$ is constant on each fibre $\mathcal S_s$, i.e. iff
$\Lambda_K$ is constant along the self-motion manifolds over $\mathcal C$ [P]. That holds when $\Pi$ restricted over
$\mathcal C$ is a Riemannian submersion for some left-invariant metric; a serial arm does not satisfy this in general.
On the Panda, $\tfrac12\oint U$ differs by 9 % between the two lifts in the table of Section 4.3 (7.255 vs 6.600), so $U$ is far from constant on the fibres.

**Which one to use.** Both are geometrically natural. The comparison favours the **full** energy for the paper:

| | full $\mathcal E_K$ | vertical $\mathcal E^v_K$ |
|---|---|---|
| physical meaning | passive dragged arm; zero null-space torque (Thm B) | least self-motion |
| contains closed pseudoinverse lifts | no (they are not critical in general) | yes, as zero-energy minimisers |
| coercive with joint-limit barrier (Sec. 8) | yes | no (arc fibres) |
| stiffness on the Panda [N] | lowest Hill eigenvalue 2.94; multipliers 3.81, 0.262 | 0.072; multipliers 1.23, 0.81 (nearly degenerate) |
| distance from rest lift on the Panda [N] | 0.64–0.87 rad | 0.14–0.24 rad |

The near-degeneracy of $\mathcal E^v$ on the Panda is explained by Prop. D below. If transport were isometric, the
vertically harmonic lifts would form a continuous family (exactly like $Q_0$ in Theorem A). The Panda's transport is
not isometric ($|\kappa|$ up to 0.32), but $\mathcal E^v$ is only weakly pinned.

### 4.5 A reparametrisation-invariant variant (remark)

$\mathcal E_K$ depends on the arc-length parametrisation of $\gamma$, hence on $W$. Lemma 1.1 identifies sections with
closed curves in the self-motion tube $\mathsf Z\subset\mathsf Q$ that cover $\mathcal C$ once. A parametrisation-free
alternative is a **closed geodesic of $(\mathsf Z,K|_{\mathsf Z})$** in the class of a section. For $K=M$, these are
the periodic motions of the passive arm whose flange slides freely on a frictionless rail of shape $\mathcal C$
(Jacobi). Reparametrised by $\eta$, such a geodesic is a harmonic lift for the parametrisation of $\mathcal C$
proportional to its own $K$-arclength [P, since $(D_sq')^V=\tau''\dot c^V$ under $s\mapsto\tau(s)$]. Existence of a
length minimiser in the class is standard in the compact or barrier-complete case [PS]. However, it is **not known to
be a section**: it might be tangent to a fibre somewhere, where the flange would stop on the rail [C]. Harmonic
sections are sections by construction. This variant is not pursued further.

### 4.6 Second variation; nondegeneracy

**Proposition 4.5 [P].** At a harmonic lift $q$, for vertical $w$ along $q$,

$$
d^2\mathcal E_K(w,w)=\oint\Big(|D_sw|_K^2-K\big(R(w,q')q',w\big)-K\big(\mathrm{II}(w,w),D_sq'\big)\Big)ds .
$$

For $r=1$ and $w=\phi x$ this is a **periodic Sturm–Liouville (Hill) form**,

$$
d^2\mathcal E_K(\phi x,\phi x)=\oint\big(\phi'^2+\mathcal R\,\phi^2\big)ds,\qquad
\mathcal R=|D_sx|_K^2-K(R(x,q')q',x)-K(\kappa,D_sq').
$$

For $K=I$: $\mathcal J\phi=-\phi''+(|x'|^2-\kappa\cdot q'')\phi$.

*Proof.* Take a variation $q_\epsilon$ through sections. Then
$\frac{d^2}{d\epsilon^2}\mathcal E=\oint|D_sw|^2+K(D_\epsilon D_s\partial_\epsilon q,q')$ and
$D_\epsilon D_s\partial_\epsilon q=D_sD_\epsilon\partial_\epsilon q+R(\partial_\epsilon q,q')\partial_\epsilon q$.
Integrating by parts, $\oint K(D_sD_\epsilon\partial_\epsilon q,q')=-\oint K(D_\epsilon\partial_\epsilon q,D_sq')$. At a
harmonic lift $D_sq'$ is horizontal, and the horizontal part of the acceleration $D_\epsilon\partial_\epsilon q$ of a
curve inside a fibre is $\mathrm{II}(w,w)$. For $r=1$, $|D_s(\phi x)|^2=\phi'^2+\phi^2|D_sx|^2$ because
$K(x,D_sx)=0$. $\square$

$q$ is **nondegenerate** if $0$ is not an eigenvalue of the periodic Jacobi operator, and its **Morse index** is the
number of negative periodic eigenvalues. Consequences of nondegeneracy:

- **Smooth dependence [PS].** The harmonic lift persists, and is locally unique and smooth, under $C^\infty$-small
  perturbations of $\gamma$, $K$ or added potentials. *Sketch:* the Euler–Lagrange operator is a smooth Fredholm map
  of index 0 between $C^{2,\alpha}$ and $C^{0,\alpha}$ sections; nearby bundles are identified by a fibre-preserving
  diffeomorphism close to the identity, and the implicit function theorem applies.
- **Numerics.** Newton's method on the discretisation converges quadratically (Section 9).

---

## 5. Theorem C ($r=1$): existence, ordering, canonical choice

**Setting C.** $r=1$. $\mathcal S^\circ\subset\mathcal S$ is an open union of fibre components that is a fibre bundle
over $S^1$ with fibre $F=S^1$ (Prop. 1.2) or $F=\mathbb R$ (Lemma 1.4). It is oriented by $V$, and therefore trivial
[P]: $\operatorname{Diff}^+(S^1)$ and $\operatorname{Diff}^+(\mathbb R)$ are connected, so the monodromy is isotopic to
the identity. Fix a trivialisation $\Psi:S^1\times F\to\mathcal S^\circ$, write $q(s)=\Psi(s,u(s))$, and for
$F=\mathbb R/\ell\mathbb Z$ use lifts $\tilde u$ with $\tilde u(s+L)=\tilde u(s)+k\ell$ (class $k$).

**(C1) Reduction [P].** $\mathcal E_K(q)=\oint\mathcal L(s,u,u')\,ds$ with

$$
\mathcal L(s,u,p)=\tfrac12a(s,u)p^2+b(s,u)p+\tfrac12c(s,u),\quad a=|\partial_u\Psi|_K^2>0,\ b=K(\partial_s\Psi,\partial_u\Psi),\ c=|\partial_s\Psi|_K^2 .
$$

This is a smooth Lagrangian, $L$-periodic in $s$ and strictly convex and quadratic in $p$: a time-periodic Tonelli
Lagrangian with one degree of freedom. In transport coordinates (when the connection is complete),
$b=0$, $a=\rho^2$ (the stretch factor, Section 6) and $c=U$. For $\mathcal E^v$, replace $c$ by $b^2/a$; it is the
quadratic form $\tfrac12a(p+b/a)^2$.

**Theorem C.** Assume Setting C.

1. **Existence [PH].** (a) $F=S^1$: every class $k$ contains a minimiser of $\mathcal E_K$ (and of $\mathcal E^v_K$).
   (b) $F=\mathbb R$: assume (B): every sublevel set $\{\mathcal E_K\le C\}$ of sections of $\mathcal S^\circ$ lies in
   a compact subset of $\mathcal S^\circ$. Then a minimiser exists. (B) holds for the complete barrier metrics and the
   barrier potentials of Section 8. In all cases minimisers are $C^\infty$ harmonic lifts (Lemma 4.1).
2. **Crossing lemma [P]** (Aubry and Le Daeron 1983; Mather 1982; Bangert 1988). Two distinct minimisers of the same
   class do not intersect. For $F=\mathbb R$ they are strictly ordered, $u_1<u_2$ everywhere. For $F=S^1$ the lifts of
   all minimisers of class $k$, together with all their deck translates $\tilde u+j\ell$, are pairwise strictly
   ordered.
3. **Order structure [P].** The set $\mathcal M_k$ of minimisers of class $k$ is compact in $C^2$ and totally
   ordered. For $F=\mathbb R$ it has a least and a greatest element $u_-\le u_+$. These are **canonical even without
   uniqueness**, once the orientation of $V$ is fixed (here by the generalised cross product). For $F=S^1$ there is no
   least element: the minimisers form a cyclically ordered family of disjoint curves on the torus, and choosing one
   needs a reference.
4. **Genericity [P].** For a fixed class, and under (B) uniformly for the perturbed problems, the set of potentials
   $P\in C^2_b(\mathcal S^\circ)$ for which $\mathcal E_K+\oint P$ has a *unique and nondegenerate* minimiser is open
   and dense. For $F=S^1$, intersecting over $k$ gives a residual set. Whether a generic *path* $\gamma$ (rather than a
   generic potential) has a unique minimiser is not proven [C].
5. **Hyperbolicity [P, using Hill's oscillation theorem K].** Let $u^\star$ be a nondegenerate harmonic lift with
   Morse index $\iota$, and let $M$ be the monodromy of its Jacobi equation, i.e. the linearised Poincaré (period) map
   of the Euler–Lagrange flow. Then $\det M=1$ and
   $$\operatorname{sign}(\operatorname{tr}M-2)=(-1)^{\iota}.$$
   A nondegenerate (local) minimiser is therefore a **hyperbolic** periodic orbit with real positive multipliers
   $\mu,1/\mu$, $\mu>1$. A degenerate minimiser is parabolic ($\operatorname{tr}M=2$). An index-1 harmonic lift is
   elliptic or inversely hyperbolic.
6. **Minimax lifts [PS; K for twist maps].** If $u_1<u_2$ are strict local minimisers of the same class with no
   critical point strictly between, there is a mountain-pass harmonic lift $u_1<u_3<u_2$ (index $\le1$). These are
   Birkhoff's minimax orbits (for twist maps see e.g. Katok 1982).
7. **Global canonical choice [P].** If all classes have minimisers, the class minima satisfy
   $m_k\ge\frac{a_0}{4L}k^2\ell^2-C\to\infty$, so the global minimum over all classes is attained, in finitely many
   classes. The **canonical lift** is the global minimiser. It exists, it is generically unique (item 4), and it
   selects its own winding class, which removes the arbitrariness in Prop. 2.2(b).

*Proofs.*

1. Direct method. For $F=S^1$, $a\ge a_0>0$, $|b|\le b_0$ and $c\ge0$ on the compact $S^1\times F$, so
   $\mathcal E\ge\frac{a_0}{4}\oint u'^2-\frac{b_0^2L}{a_0}$. Normalise $\tilde u_n(0)\in[0,\ell)$ by a deck
   translation. A minimising sequence is then bounded in $H^1$; a subsequence converges weakly in $H^1$ and uniformly,
   which preserves $k$. $\mathcal L$ is convex in $p$ and continuous, so $\mathcal E$ is weakly lower semicontinuous
   (Tonelli). For $F=\mathbb R$, (B) places a minimising sequence in a compact subset of $\mathcal S^\circ$, on which
   the same estimates hold. The uniform limit stays in $\mathcal S^\circ$: each $\mathcal S^\circ_s$ is closed in
   $\mathcal S_s$, being a union of components.
2. Appendix A.4.
3. A minimiser satisfies the Euler–Lagrange equation. With the energy bound this gives a uniform $C^2$ bound on a
   compact set (by (B)). By Arzelà–Ascoli and lower semicontinuity, limits of minimisers are minimisers. For
   $F=\mathbb R$, total order and compactness make the minimiser with the least value at $s=0$ lie below all others
   everywhere.
4. *Density.* Let $u_1$ be a minimiser for $P$. Add $P_\varepsilon(s,u)=\varepsilon\psi(u-u_1(s))$, with
   $\psi(t)=t^2/(1+t^2)$ for $F=\mathbb R$ or $\psi(t)=1-\cos(2\pi t/\ell)$ for $F=S^1$. Then
   $(\mathcal E+\oint P+\oint P_\varepsilon)(u)>m$ for every $u\neq u_1$ in the class, so $u_1$ is the unique
   minimiser. Its Jacobi operator gains $\varepsilon\psi''(0)>0$, so its lowest eigenvalue becomes positive.
   *Openness.* If $u_1$ is the unique minimiser and nondegenerate, minimisers of $C^2$-close problems converge to $u_1$
   (compactness and uniqueness). Near $u_1$ the implicit function theorem gives a unique critical point, which is
   nondegenerate.
5. Appendix A.5.
6. *Sketch.* The order interval $[u_1,u_2]$ is positively invariant for the gradient (heat) flow of Section 9, by the
   parabolic comparison principle (Angenent 1988 uses this mechanism for twist maps). Palais–Smale holds on it by
   compactness. The mountain-pass theorem (Ambrosetti and Rabinowitz 1973) inside the invariant order interval gives
   $u_3$.
7. $\int_0^L\tilde u'^2\ge(k\ell)^2/L$ by Cauchy–Schwarz, combined with the lower bound in item 1. $\square$

**What needs circle fibres and what does not.**

| result | $F=S^1$ (no limits, compact) | $F=\mathbb R$ (arc fibres, joint limits) |
|---|---|---|
| homotopy classes | $\mathbb Z$-torsor (Prop. 2.2) | one class |
| existence | automatic, every class, both energies | needs (B): barrier metric or potential; full energy only |
| crossing lemma, total order | yes (on lifts, deck-invariant) | yes |
| canonical extremal minimiser | no (cyclic order) | yes ($u_\pm$) |
| genericity, hyperbolicity, parity | yes | yes |
| Aubry–Mather sets of irrational rotation number (quasi-periodic, non-closing lifts of minimal energy per lap) | yes (Mather 1982) | not defined |

**For $r\ge2$** only the following survive: existence in the compact case or under (B), by the same direct method [PH];
regularity; the second-variation formula; and smooth dependence. The crossing lemma and the order structure are
one-dimensional phenomena. Positivity of the second variation of a minimiser does not imply hyperbolicity for $r\ge2$
in general [C: we know of no general result for this periodic, non-autonomous setting].

**Panda [N].** The full-energy minimisation ($K=I$) finds exactly the structure of items 2–6 along the single class:

- global minimiser: offset $0.64$–$0.87$ rad above the rest-posture lift, $\mathcal E=6.6691$ ($m=250$);
- index-1 harmonic lift: offset $1.39$–$1.96$, $\mathcal E=7.4596$, $\operatorname{tr}M=-1.33$ (elliptic);
- second local minimiser: offset $1.97$–$2.40$, $\mathcal E=6.7835$, multipliers 4.14 and 0.242 (hyperbolic), only
  0.062 rad from joint 2's limit.

The three are strictly ordered. Every index-0 point has $\operatorname{tr}M>2$ and the index-1 point has
$|\operatorname{tr}M|<2$, as item 5 predicts. Multi-start from constant offsets $-1,0,1,2$ found no other minimiser.
Global minimality of the first is therefore established only numerically.

---

## 6. Proposition D: isometric transport and the drift formula

Assume the $K$-horizontal transport $\Theta_s$ is complete (Prop. 1.3), $F=\mathcal S_0$ with the metric induced by
$K$, $T=\Theta_L|_F$. Sections are $q(s)=\Theta_s(c(s))$ with (2.1).

**Lemma 6.1 [P].** For the horizontal lift $Y$ of the path and vertical $v,w$,
$(\mathcal L_YK)(v,w)=-2K(Y,\mathrm{II}(v,w))$. Hence $\Theta_s$ is a fibre isometry for all $s$ iff
$K(Y_s,\mathrm{II})\equiv0$ on $\mathcal S$: the self-motion manifolds are totally geodesic **in the direction of the
path**. This holds in particular if they are totally geodesic, which is Hermann's (1960) situation for Riemannian
submersions [K]; see also Vilms (1970).

*Proof.* $(\mathcal L_YK)(v,w)=K(\nabla_vY,w)+K(v,\nabla_wY)$, and
$K(\nabla_vY,w)=v\,K(Y,w)-K(Y,\nabla_vw)=-K(Y,\mathrm{II}(v,w))$, because $K(Y,w)\equiv0$. The fibre metric pulled
back by the flow of $\partial_s+Y_s$ has $s$-derivative $\mathcal L_YK|_V$. $\square$

**Proposition D.**

1. **[P] Isometric transport, vertical energy.** If $\Theta_s$ are fibre isometries, then
   $\mathcal E^v(q)=\tfrac12\int_0^L|c'|_F^2ds$. The vertically harmonic lifts are exactly
   $q(s)=\Theta_s(c(s))$ with $c$ a constant-speed geodesic of $F$ satisfying $c(0)=T(c(L))$ and $c'(0)=dT\,c'(L)$:
   **horizontal lift plus uniform drift along a geodesic of the fibre**. The minimisers in a class are the geodesics
   realising the minimal displacement $\min_y d_F(y,T(y))$ in that class, with energy $d^2/(2L)$.
   - $r=1$, $F$ a circle of length $\ell$, $T$ = rotation by $\tau$: $c(s)=u_0-\frac sL(\tau+k\ell)$ with $u_0$
     *arbitrary* and energy $(\tau+k\ell)^2/(2L)$. The global minimiser has $|\tau+k\ell|$ minimal, i.e. the
     *principal* drift. It is unique modulo $u_0$ unless $\tau\equiv\ell/2$ mod $\ell$, the analogue of
     $-1\in\operatorname{spec}H$.
   - $F=\mathbb R$, $T$ = translation by $\tau$: $c(s)=u_0-\frac sL\tau$.
   - $F$ a bounded open interval: $T=\mathrm{id}$, and closed horizontal lifts exist.
   - **Theorem A is the case** $F=O(5)$ with the bi-invariant metric, $T(Q)=HQ$. Every point has displacement
     $|\Phi|_F$, so every $Q_0$ gives a minimiser.
2. **[P] Isometric transport, full energy.** $\mathcal E=\mathcal E^v+\tfrac12\oint U$. If $U(s,\cdot)$ is constant on
   fibres, the conclusions of item 1 hold for $\mathcal E$. Otherwise the Euler–Lagrange equation is
   $D_sc'=\tfrac12\operatorname{grad}_FU(s,c)$: Newtonian motion in the time-periodic potential $-U/2$ on $F$, with
   **non-uniform** drift. The lead's formulation of Prop. D ("the harmonic lift is horizontal lift plus uniform drift")
   is therefore **false for the full energy** unless $\Lambda_K$ is constant along the fibres.
3. **[P] General transport, $r=1$ (drift formula).** Let $u$ be $K$-arclength on $F$ and
   $\rho(s,u):=|d\Theta_s(\partial_u)|_K$ the stretch factor. Then
   $$\partial_s\log\rho(s,u)=-K(Y_s,\kappa)\big(\Theta_s(u)\big),\qquad
   \rho(s,u)=\exp\Big(-\int_0^sK(Y_\sigma,\kappa)(\Theta_\sigma(u))\,d\sigma\Big),$$
   $$\mathcal L(s,u,u')=\tfrac12\rho(s,u)^2u'^2+\tfrac12U(s,u),\qquad (\rho^2u')'=\rho\rho_uu'^2+\tfrac12U_u .$$
   In terms of $v=\rho u'$ this is (4.3). If $\rho=\rho(s)$ does not depend on $u$ (conformal transport) and $U_u=0$,
   the vertically harmonic lift is explicit:
   $$u(s)=u_0-(\tau+k\ell)\,\frac{\int_0^s\rho^{-2}}{\int_0^L\rho^{-2}} .$$
   The drift is distributed with density $\propto\rho^{-2}$: it moves along the self-motion where transport has
   stretched the fibre least.

*Proof.* (1) $\mathcal Vq'=d\Theta_s(c')$, because $\Theta_s$ maps fibres to fibres and $Y_s$ is horizontal; isometry
gives $|\mathcal Vq'|=|c'|_F$. The closure conditions are (2.1). The Euler–Lagrange equation of
$\tfrac12\int|c'|^2$ with the end condition $c(0)=T(c(L))$ is the geodesic equation, together with the natural
boundary condition $c'(0)=dT\,c'(L)$ obtained from varying $c(L)$. Minimising over $c(L)$ gives the
minimal-displacement statement. On a circle $T$ is a rotation, on $\mathbb R$ a translation, and on a bounded interval
the only orientation-preserving isometry is the identity.
(2) By (4.1). (3) $\partial_s\rho^2=(\mathcal L_YK)(\hat\partial_u,\hat\partial_u)\rho^2=-2K(Y,\kappa)\rho^2$ by
Lemma 6.1 with $\mathrm{II}(x,x)=\kappa$. The Lagrangian follows from $|q'|^2=|Y|^2+\rho^2u'^2$, using $b=0$ in
transport coordinates. If $\rho_u=0$ and $U_u=0$, then $\rho^2u'$ is constant, and closure fixes the constant. $\square$

**Panda [N].** Transport is not isometric: along the rest-posture lift $|\kappa|\in[0.065,0.32]$ and
$-K(Y,\kappa)\in[-1.24,0.92]$. Item 1 still explains the qualitative picture. The pseudoinverse displacement $T(t)-t$
is smallest near the rest-posture lift ($\approx-0.24$), and the vertically harmonic lift sits there (offset
$-0.24..-0.14$). Its energy $\mathcal E^v=0.0328$ is close to the isometric prediction $0.2423^2/(2L)=0.0381$. It is
nearly degenerate (Section 4.4).

---

## 7. Proposition E: the rest-posture lift as a singular limit; the lift in the controller

### 7.1 The family $\mathcal E_k$

Let $f\ge0$ be a smooth function on $\mathcal S^\circ$. The paper's choice is $f(s,q)=\tfrac12|q-q_r|^2$. Assume:

- (E1) for each $s$, $f(s,\cdot)$ has a unique minimiser $q_f(s)$ on $\mathcal S^\circ_s$, and it is nondegenerate:
  $\partial_x^2f>0$ along the fibre. By the implicit function theorem $q_f$ is then a smooth section, the
  **pointwise lift**; with $f=\frac12|q-q_r|^2$ it is the paper's rest-posture lift.
- (E2) (B) holds for $\mathcal E_k:=\mathcal E_K+k\oint f$ uniformly in $k\ge0$ (e.g. with a barrier, Section 8).

**Proposition E.**

1. **[PH]** For every $k\ge0$, $\mathcal E_k$ has minimisers $q_k$ (Theorem C1).
2. **[PH] (singular limit)** $q_k\to q_f$ uniformly as $k\to\infty$, for every choice of minimisers. *Proof:*
   Appendix A.3.
3. **[PS] (rate)** With $h(s)=\partial_x^2f(s,q_f(s))$ (for $f=\tfrac12|q-q_r|^2$ and $K=I$,
   $h=1+(q_f-q_r)\cdot\kappa$),
   $$q_k=q_f+\frac1k\,\frac{K(D_sq_f',x)}{h}\,x+O(k^{-2})\qquad\text{uniformly}. \tag{7.1}$$
   *Sketch.* Write $q=\Psi(s,u_f+v)$. The Euler–Lagrange equation is $f_u(s,u_f+v)=\frac1kG(s,u,u',u'')$, with $G$
   the (smooth) Euler–Lagrange expression of $\mathcal E_K$. Linearise: $(-\frac1k\partial_s P\partial_s+h)v=
   \frac1kG_0+\dots$. The operator is invertible with norm $\le1/\min h$ in $C^0$ (maximum principle) and in the
   $k$-scaled $C^2$ norm. A contraction argument on the circle (no boundary layers) gives $v=\frac{G_0}{kh}+O(k^{-2})$.
   The $O(k^{-2})$ constant involves $(G_0/h)''$, so the asymptotic regime requires $k^{-1/2}$ to be small compared
   with the scale on which $q_f''$ varies.
4. **[P] (continuation)** Along any branch of nondegenerate minimisers, $k\mapsto q_k$ is smooth (implicit function
   theorem). If the branch from $k=\infty$ stays nondegenerate down to $k=0$, it ends at a harmonic lift. That this lift
   is the *global* minimiser is not guaranteed [C].

The two endpoints are canonical: the paper's lift at $k=\infty$ (pointwise optimal, ignores the motion) and the
harmonic lift at $k=0$ (dynamically natural, ignores the rest posture). For finite $k$, $\mathcal E_k$ minimisers
correspond to the dragged passive arm with a joint-space spring of stiffness $k\nu^2$ towards $q_r$.

**Comparison with the null-space-pull family.** The lead's closed orbits of $q'=J^+\zeta_\gamma+k(I-J^+J)(q_r-q)$ are
the attracting periodic orbits of a dissipative first-order system. They form a different family. They also tend to
$q_f$ as $k\to\infty$; the analogous corrector is $\frac1k(x\cdot q_f')/h$, a velocity term with inner scale $1/k$
rather than $1/\sqrt k$ [PS]. As $k\to0$ they cannot converge to a closed lift on the Panda [PS]: their Poincaré maps
$T_k\to T$ uniformly on compact parts of the admissible arc, and $T(t)-t\le-0.22$ there, so for $k$ below some threshold
no admissible closed orbit remains. The variational family $\mathcal E_k$ has a
closed limit at both ends.

**Panda [N]** ($K=I$, $m=250$, continuation from the rest-posture lift):

| $k$ | $10^4$ | $10^3$ | 300 | 100 | 30 | 10 | 3 | 1 | 0.3 | 0 |
|---|---|---|---|---|---|---|---|---|---|---|
| $\max_s\lvert q_k-q_f\rvert$ (rad) | 1.82e-3 | 1.40e-2 | 3.14e-2 | 5.61e-2 | 0.109 | 0.221 | 0.450 | 0.668 | 0.800 | 0.871 |
| $k\cdot$dist | 18.2 | 14.0 | 9.4 | 5.6 | 3.3 | 2.2 | 1.4 | 0.67 | 0.24 | — |
| $\mathcal E_K(q_k)$ | 7.333 | 7.298 | 7.261 | 7.225 | 7.161 | 7.033 | 6.825 | 6.706 | 6.674 | 6.669 |
| lowest Hill eigenvalue | 8519 | 897 | 274 | 93 | 29.3 | 11.2 | 5.07 | 3.55 | 3.11 | 2.94 |

The predicted constant in (7.1) is $\max_s|x\cdot q_f''/h|=18.9$; the measured $k\cdot$dist approaches it from below
(18.2 at $k=10^4$). The branch stays nondegenerate and ends at the global minimiser. For comparison, the lead's
null-space-pull orbits were 0.36 / 0.030 / 1.3e-3 rad from $q_f$ at $k=1/10/100$.

### 7.2 The harmonic lift in the controller

Nothing in Sections IV–V of the paper changes [P].

- **Assumption 1.** A minimising harmonic lift in $\mathcal S^\circ\subset\mathsf Q_{\rm adm}$ is $C^\infty$ (Lemma 4.1),
  closed (by construction), regular and inside the limits (by admissibility, guaranteed with a barrier).
- **Lemma 2 and Lemma 3.** $X$ is built exactly as before. For $r=1$ it is the normalised generalised cross product
  along $q^\star$, periodic with no integration. $\beta=X(\eta)^\top(q-q^\star(\eta))$, $d\beta$ (eq:dbeta), $A(q)$ and
  the feedback (eq:fb) are unchanged.
- **Derivatives.** $q^{\star\prime}$ comes from the discrete solution, and $q^{\star\prime\prime}$ comes exactly from
  (4.2): $q^{\star\prime\prime}=J_b^+(\zeta_\gamma'-\dot J_bq^{\star\prime})$ for $K=I$. So Hermite data of higher
  order are available without numerical differentiation.
- **New interpretation.** For $K=M$, at constant path speed and with $\xi\equiv0$, $\beta\equiv0$, the torque that
  realises the motion is $J_b^\top f+G$ (Theorem B(g)). The redundancy channel requires **no null-space torque in
  steady state**; the $\beta$-loop only rejects deviations. For $K=I$ the analogous statement is zero self-motion
  acceleration, $X^\top\ddot q=0$.
- **Stability.** The lift is an unstable (hyperbolic) orbit of the *passive* dynamics. The feedback (eq:beta-law)
  stabilises $\beta=0$ regardless, so Theorem 2 of the paper is unaffected.

---

## 8. Joint limits and singularities

Let $\Omega=\mathsf Q_{\rm adm}=\mathsf Q_{\rm reg}\cap\operatorname{int}B$ (bounded) and $d(q)=\operatorname{dist}
(q,\partial\Omega)$ (Euclidean), with $\partial\Omega\subset\partial B\cup(\Sigma\cap B)$.

**Proposition 8.1 (complete conformal barrier metric) [P].** Let

$$
\varphi(q)=1+\sum_{i=1}^n\frac{c_i}{(q_i-q_i^-)(q_i^+-q_i)}+\frac{c_\Sigma}{\sqrt{\det\big(J_b(q)J_b(q)^\top\big)}},\qquad c_i,c_\Sigma>0,
$$

and $K=\varphi^2K_0$ with $K_0\ge\underline k I$ smooth on $\bar B$ (e.g. $I$ or $M$). Then:

- (a) $\varphi\ge c/d$ on $\Omega$, so $K$ is a **complete** metric on $\Omega$.
- (b) $H_K=H_{K_0}$: the connection, and hence pseudoinverse transport, is unchanged.
- (c) $U_K=\varphi^2U_{K_0}\to\infty$ at $\partial\Omega$, and hypothesis (B) holds for $\mathcal E_K$. So Theorem C1(b)
  applies, and minimising harmonic lifts exist and are admissible.
- (d) $\mathcal E^v_K$ is in general **not** coercive: there is no confining term along arc fibres. This is the reason
  the full energy is preferred (Correction 3).

*Proof.* (a) If the nearest boundary point lies on the face $q_i=q_i^\pm$, then $d=\min(q_i-q_i^-,q_i^+-q_i)$ and the
$i$-th term is $\ge c_i/((q_i^+-q_i^-)d)$. If the nearest boundary point $\bar q$ lies in $\Sigma$, then
$\sigma_{\min}(J_b(\bar q))=0$. $\sigma_{\min}$ is 1-Lipschitz in the operator norm (Weyl) and $J_b$ is
$\ell_J$-Lipschitz on $B$, so $\sigma_{\min}(J_b(q))\le\ell_Jd$ and
$\sqrt{\det J_bJ_b^\top}=\prod\sigma_i\le\sigma_{\max}^5\sigma_{\min}\le C\ell_Jd$. Completeness: a curve leaving
every compact subset of $\Omega$ has $d\to0$, and its $K$-length is $\ge\int\varphi|dq|\ge\int c\,|d(d)|/d=\infty$
($|\nabla d|\le1$). By Hopf–Rinow, $(\Omega,K)$ is complete.
(b) Section 1.2.
(c) $\Lambda_{K_0}\ge\lambda_{\max}(J_bK_0^{-1}J_b^\top)^{-1}I\ge c_0I$ on the bounded $B$, and
$\|\zeta_\gamma\|\ge c\|\zeta_\gamma\|_W=c$, so $U_{K_0}\ge u_0>0$ and $U_K\ge u_0\varphi^2$. If $\mathcal E_K(q)\le C$,
then some $q(s_0)$ lies in the compact set $\{u_0\varphi^2\le2C/L\}$, and the $K$-length of the loop is
$\le\sqrt{2LC}$. The loop therefore lies in a closed $K$-ball around a compact set, which is compact by completeness.
Uniform limits stay in $\mathcal S^\circ$ (proof of C1).
(d) An example suffices. Take a bundle $S^1\times\mathbb R$ with complete fibre metric $du^2$ and horizontal transport
$u'=-\epsilon(u)$, $\epsilon(u)=1/(1+u^2)$. Then $\mathcal E^v(u)=\frac12\oint(u'+\epsilon(u))^2ds$. Constant sections
$u\equiv c$ give $\frac L2\epsilon(c)^2\to0$ as $c\to\infty$, so $\inf\mathcal E^v=0$. It is not attained: zero energy
would need $u'=-\epsilon(u)<0$ throughout, which is impossible for a periodic $u$. (This is an abstract bundle with
connection; that a serial arm realises such behaviour is not claimed.) $\square$

**Proposition 8.2 (barrier potential) [P].** Let $K$ be smooth on $\bar B$ (e.g. $I$ or $M$) and $\mathsf B:\Omega\to
[0,\infty)$ smooth with $\mathsf B\ge c/d^2$. If $\mathcal E_K(q)+\oint\mathsf B(q)\le C$, then
$\min_sd(q(s))\ge\delta(C)>0$. Hence (B) holds for $\mathcal E_K+\oint\mathsf B$, and its minimisers exist and are
interior critical points: their Euler–Lagrange equation is Theorem B with $+\operatorname{grad}^V\mathsf B$. *Proof:*
Appendix A.6 (a logarithmic divergence of $\oint\mathsf B$ near a near-contact).

**Proposition 8.3 [P; PS for the global part].** If $\mathcal E_K$ without barrier has a nondegenerate minimiser $q^\star$
in the interior of $\Omega$, then for barrier strengths $\varepsilon\to0$ ($\varepsilon\mathsf B$, or
$\varepsilon c_i,\varepsilon c_\Sigma$) there are nondegenerate critical points $q^\star_\varepsilon\to q^\star$
(implicit function theorem). If moreover $q^\star$ is the strict global minimiser and the infimum is not approached at
$\partial\Omega$, the barrier global minimisers converge to $q^\star$ [PS].

**Panda [N].** The barrier-free minimiser ($K=I$) is interior, with smallest joint margin 0.3596 rad (joint 4) and
$\sigma_{\min}(J_b)\ge0.159$. So the paper's example does not need a barrier, and a weak barrier changes the lift only
slightly (Prop. 8.3). The safety margin of the paper's `SafetyMonitor` (0.10 rad) can be built in by putting the barrier
at the shrunk box $q_i^\pm\mp0.10$. The second local minimiser (margin 0.062 rad at joint 2) would then be excluded or
pushed away.

**Not proven.** (i) Existence of a minimiser for $K=I$ on $\Omega$ **without** a barrier: a minimising sequence could
approach $\partial\Omega$, and only the numerics show it does not on the Panda. (ii) Quantitative bounds on how far the
barrier moves the lift for given $c_i$. (iii) Avoidance of the *algorithmic* singularities of the paper's
$A(q)$ (Proposition 2) by the lift: these concern $\sigma_{\min}(A)$, not $\sigma_{\min}(J_b)$, and are not addressed by
the barrier on $J_b$.

---

## 9. Computation

### 9.1 Discrete problem

Use the grid $s_j=jh$, $h=L/m$, $j\in\mathbb Z/m\mathbb Z$ (the paper's grid), with unknowns $q_j\in\mathcal S_{s_j}$ and
the discrete Lagrangian (variational-integrator style; Marsden and West 2001)

$$
\mathcal E_h=\sum_j\frac1{2h}\,|q_{j+1}-q_j|^2_{K(q_{j+1/2})}\quad\big(+\,h\textstyle\sum_jP(s_j,q_j)\ \text{for a potential or barrier}\big),
$$

which is exactly $\sum_j|q_{j+1}-q_j|^2/(2h)$ for $K=I$. There is **one unknown per node** in the chart
$u_j\mapsto Q_j(u_j)$, defined by $\Pi(Q_j)=\gamma(s_j)$ and $x_j^\top(Q_j-\bar q_j)=u_j$, with $x_j=x(\bar q_j)$. It is
solved by Newton on the $7\times7$ system $[J_b;x_j^\top]$, which is the paper's own level-set / start-configuration
solver. The chart is re-centred at the current iterate $\bar q$ in every Newton step, so $u_j$ is the paper's
$\beta$-offset from the current iterate. Re-centring avoids any validity problem of the hyperplane chart for large
offsets (on the Panda the solution is 0.87 rad from the rest-posture lift).

### 9.2 Exact derivatives ($K=I$) [P]

At $u=0$: $Q_j'=x_j$ and $Q_j''=\kappa_j=-J_b^+\big(dJ_b[x_j]\,x_j\big)$, the curvature vector of the self-motion
curve. Differentiate $J_b(Q)Q'=0$ and $x_j^\top Q''=0$; $J_b^+$ maps into $x_j^\perp$. Hence

$$
g_j=\tfrac1h\,(2q_j-q_{j-1}-q_{j+1})\cdot x_j,\qquad
H_{jj}=\tfrac1h\big[2+(2q_j-q_{j-1}-q_{j+1})\cdot\kappa_j\big],\qquad H_{j,j\pm1}=-\tfrac1h\,x_j\cdot x_{j\pm1}.
$$

$g=0$ is the discrete form of $X^\top q''=0$. $H$ is **cyclic tridiagonal**, with $H\approx h\,\mathcal J$ (the Hill
operator of Prop. 4.5). Each Newton step costs $m$ retractions (2–3 $7\times7$ solves each), $m$ kernel and curvature
evaluations, and an $O(m)$ cyclic tridiagonal solve (Sherman–Morrison on the two corner entries, or cyclic Thomas). For
$K=M$ the formulas acquire $dM$ terms and $K(q_{j+1/2})$ weights, but the structure is the same.

### 9.3 Newton, globalisation, verification

Take $\Delta u=-(H+\mu I)^{-1}g$ ($\mu=0$ when $H\succ0$), backtrack on $\mathcal E_h$, and retract. Convergence is
quadratic near a nondegenerate minimiser (standard). On the Panda, starting from the rest-posture lift, it takes 4
iterations to $|g|_\infty<2\cdot10^{-12}$ for both $m=500$ and $m=1000$ [N]. Verify:

1. $H\succ0$ (lowest eigenvalue $>0$): nondegenerate minimiser;
2. the **discrete Floquet test**: the monodromy of the discrete Jacobi equation
   $b_{j-1}\phi_{j-1}+a_j\phi_j+b_j\phi_{j+1}=0$ is the product of the transfer matrices
   $\begin{bmatrix}-a_j/b_j&-b_{j-1}/b_j\\1&0\end{bmatrix}$, with determinant 1. Expect $\operatorname{tr}>2$ at a
   minimiser (discrete C5; compare MacKay and Meiss 1983 for the relation between $\det H$ and the residue);
3. multi-start from constant offsets, to probe global minimality (C2–C7 say what can be found);
4. closure of one lap of the ODE (4.2) from $(q_0,q_0')$, which must close up to $|\mu|\cdot$(discretisation error).

### 9.4 Frenkel–Kontorova structure [P; K for the consequences]

$\mathcal E_h=\sum_j\mathsf H_j(u_j,u_{j+1})$ with $\partial^2\mathsf H_j/\partial u_j\partial u_{j+1}=-\frac1h
x_j\cdot x_{j+1}<0$ whenever consecutive kernel directions are positively aligned. This is the **twist condition** of
Aubry and Le Daeron's generalised Frenkel–Kontorova chain, here on a ring of $m$ sites. Aubry's discrete crossing lemma,
ordering and minimax theory therefore apply to the discrete problem [K; Aubry's arguments for periodic
configurations carry over to site-dependent, $m$-periodic $\mathsf H_j$, which we have not written out]. Caveat: twist is verified along the iterates
(small $|u_{j+1}-u_j|$, which bounded energy enforces for small $h$), not globally on $F\times F$. On the Panda
$\min_jx_j\cdot x_{j+1}=0.99998$ ($m=250$) [N].

### 9.5 Harmonic-map heat flow (globally convergent alternative)

The negative $L^2$-gradient flow of $\mathcal E_K$ on sections (Eells and Sampson 1964 for maps) is

$$
\partial_tq=\mathcal V\,D_sq'\qquad\Big(\tfrac{d}{dt}\mathcal E_K=-\oint|\mathcal VD_sq'|^2\Big). \tag{9.1}
$$

In a trivialisation ($r=1$) it is the scalar, uniformly parabolic, quasilinear equation
$a\,\partial_tu=\partial_s(\mathcal L_p)-\mathcal L_u$ on the circle [P]. Consequences:

- energy decreases [P];
- the **comparison principle** for scalar parabolic equations makes the flow order-preserving, so order intervals
  between a sub- and a supersolution are invariant [K]. The discrete flow $\dot u=-D^{-1}g$, with a positive diagonal mass matrix $D$ (e.g. $D=hI$), is a *cooperative* system
  under the twist condition, hence monotone (Kamke's condition; Smith 1995) [K];
- under (B), solutions exist for all time and stay in a compact set [PS]. For real-analytic data (the Panda's forward
  kinematics, the coning circle, and polynomial or rational barriers are analytic), the Łojasiewicz–Simon inequality
  (Simon 1983) gives **convergence to a single harmonic lift** [PS];
- the limit is a critical point, generically a local minimiser. Global minimality is not guaranteed; started from the
  rest-posture lift it found the global minimiser on the Panda (same as Newton) [N].

In practice a linearly implicit Euler step $u^{n+1}=u^n-(\Delta t^{-1}D+H)^{-1}g$ interpolates between the heat flow
(small $\Delta t$) and Newton ($\Delta t\to\infty$). This is the recommended globalisation.

### 9.6 Continuation in $k$

Solve the pointwise lift (the paper's existing code), then follow the $\mathcal E_k$ branch from large $k$ down to
$k=0$ (Prop. E4), with a few Newton steps per $k$ value. On the Panda the branch is nondegenerate throughout (table in
7.1) [N]. This gives a robust initialisation that does not rely on the rest-posture lift being in the basin of
attraction of Newton.

---

## 10. Numerical stress test (Panda, coning circle, $\ell=0.1$ m, $K=I$)

All values are from the scripts in this folder ([N]). The rest-posture lift reproduces the design doc
($q^\star(0)=(0.01457,-0.13186,-0.02614,-2.21552,0.22458,2.24779,0.6471)$, margin 0.35918 rad at joint 4 near
$s=0.367$, closure $6\cdot10^{-15}$).

| quantity | rest-posture lift | harmonic lift (full energy) | vertically harmonic |
|---|---|---|---|
| $\mathcal E$ ($m=1000$; $m=500$) | 7.33978 | 6.66956; 6.66947 | 7.5621 ($m=250$) |
| $\tfrac12\oint U$ / $\mathcal E^v$ | 7.2552 / 0.0846 | 6.6005 / 0.0691 | 7.529 / 0.0328 |
| offset from rest lift along self-motion | 0 | +0.639 … +0.868 rad | −0.236 … −0.144 rad |
| $\max_s\lvert q-q_{\rm rest}\rvert$ | 0 | 0.8707 rad (at $s=0.092$) | 0.2364 |
| smallest joint margin | 0.3592 (joint 4) | 0.3596 (joint 4) | 0.3617 (joint 4) |
| smallest $\sigma_{\min}(J_b)$ | 0.2070 | 0.1589 | 0.2039 |
| self-motion speed $v=x\cdot q'$ (range) | — | $-0.035$ … $0.547$ | $0.268$ … $0.330$ (never 0: Cor. 4.4) |
| $\max\lvert q'\rvert$ | 5.28 | 5.02 | — |
| lowest Hill eigenvalue | — | 2.94 | 0.072 |
| Floquet multipliers | — | 3.813, 0.262 | 1.233, 0.811 |
| Newton iterations from rest lift | — | 4 | 4 |

Further checks:

- Pseudoinverse holonomy from $q^\star_{\rm rest}(0)$: $0.24230$ rad (2000 and 4000 RK4 steps agree to $10^{-11}$).
  On the 35 of 41 starts along the admissible arc over $\gamma(0)$ whose transport stays in the box, $T(t)-t\in[-0.65,-0.22]$; the other 6 leave the box. No admissible fixed point.
- Admissible self-motion component through the rest lift at $s\in\{0,0.1,\dots,0.7\}$: always a single arc (joint
  arclength 5.5–8.6 rad), ended at both sides by joint 2's limits, with $\sigma_{\min}(J_b)\ge0.10$ on it. This
  supports Setting C with $F=\mathbb R$; it was not checked for every $s$.
- Shooting: (4.2) with $K=I$ from the harmonic lift's data closes after one lap to $7.0\cdot10^{-6}$ (position) and
  $1.5\cdot10^{-5}$ (velocity); from the rest lift's data it misses by 1.24 rad. (The mid-lap values of about
  $10^{-3}$ printed in `run_shoot.log` compare against the grid node nearest to $s=(n+1)L/4000$. They are dominated by
  that offset in $s$, $|q'|\cdot L/4000\approx10^{-3}$. Only the end-of-lap closure is meaningful.)
- Critical points found along the single class (full energy, $m=250$): minimiser $\mathcal E=6.6691$ (index 0,
  tr 4.08); index-1 point $\mathcal E=7.4596$ (tr $-1.33$, elliptic); local minimiser $\mathcal E=6.7835$ (index 0, tr
  4.38, multipliers 4.14 / 0.242, margin 0.062 rad at joint 2). Strictly ordered by offset.
- Theorem A: see Section 3.2.
- Not computed: $K=M$. The repository has no Panda dynamic model; the identified model of the paper's robustness study
  is not in this repo.

---

## 11. What is new and what is known; open problems

**Known (cited).**

- Ehresmann connections and fibrations (Ehresmann 1951); Riemannian submersions, O'Neill's tensors, totally geodesic
  fibres (O'Neill 1966; Hermann 1960; Vilms 1970); complete metrics (Nomizu and Ozeki 1961).
- Harmonic maps and the heat flow (Eells and Sampson 1964); harmonic sections for the vertical energy (Wood 1987).
- The bundle and connection view of redundancy and non-repeatability (Klein and Huang 1983; Shamir and Yomdin 1988;
  Baker and Wampler 1988; Burdick 1989); the dynamically consistent inverse (Khatib 1987); global integral-cost
  redundancy resolution with fixed end points (Nakamura and Hanafusa 1987; Kazerounian and Wang 1988); task
  compatibility (Chiu 1988).
- Aubry–Mather theory: crossing lemma, ordering, minimax orbits (Aubry and Le Daeron 1983; Mather 1982, 1991; Katok
  1982; Moser 1986; Bangert 1988; Mañé 1996 for genericity of minimising measures); Hill's equation (Magnus and Winkler
  1966; Eastham 1973); residues and the second variation (MacKay and Meiss 1983); parabolic flows for twist maps
  (Angenent 1988); Łojasiewicz–Simon (Simon 1983); variational integrators (Marsden and West 2001).
- Rotation-minimising frames (Bishop 1975) and the uniform closing twist for closed curves (graphics: Wang et al. 2008;
  elastic rods: Bergou et al. 2008). Coulomb gauge on $S^1$ is textbook gauge theory.

**New in this note, to our knowledge.** We have not done an exhaustive literature search; treat priority claims as
provisional.

1. Defining the paper's closed lift as a minimising harmonic section of the self-motion bundle over the path loop,
   with closure built in, and the $r=1$ existence, ordering, genericity and hyperbolicity theory (Theorem C) obtained
   from Aubry–Mather / Tonelli theory.
2. Theorem A as stated: the paper's holonomy-cancelled frames are exactly the minimal-energy Coulomb gauges, unique
   modulo $O(5)$ iff $-1\notin\operatorname{spec}H$, and Lemma 3 likewise. The $SO(2)$ case for space curves is
   folklore.
3. The characterisations in Theorem B(e)–(g) in the *periodic* setting: harmonic lift = periodic orbit of
   acceleration-level resolution = passive dragged arm = zero dynamically consistent null-space torque. Combined with
   C5, this explains why the canonical lift is a hyperbolic (unstable) orbit and cannot be found by forward
   integration.
4. The full-versus-vertical distinction for non-Riemannian submersions (Prop. 4.3, Lemma 4.2): the difference is the
   operational-space inertia potential. This corrects the uniform-drift claim (Prop. D2) and gives the stretch-factor
   drift formula (D3).
5. Conformal barrier metrics leave the redundancy connection unchanged but make the full energy coercive (Prop. 8.1).
6. The discrete harmonic-lift problem is a twist (Frenkel–Kontorova) chain with a cyclic tridiagonal Newton system, and
   the Panda numerics of Section 10.

**Open problems [C].**

- (O1) For $r\ge2$: ordering and uniqueness theory, and whether nondegenerate minimisers are hyperbolic.
- (O2) Generic uniqueness with respect to the **path** $\gamma$ rather than an added potential.
- (O3) Existence on $\mathsf Q_{\rm adm}$ for $K=I$ or $K=M$ **without** a barrier, under checkable conditions.
- (O4) The $K=M$ harmonic lift of the Panda, and its effect on torques in the paper's robustness study.
- (O5) Whether the closed geodesic of the self-motion tube (4.5) is always a section.
- (O6) A joint harmonic section of the fibre product (lift, kernel frame) for $r\ge2$, where the minimal frame energy
  $\sum\theta_i(q^\star)^2/L$ depends on the lift and couples the two problems.

---

## Appendix A. Proofs

### A.1 Lemma 1.4

$\mathcal S^\circ$ is a 2-manifold. Choose a complete Riemannian metric $g$ on it (Nomizu and Ozeki 1961). The
vertical line bundle is oriented ($V$ is), so there is a smooth $g$-unit vertical field $\hat\nu$. *Its flow is
complete.* An integral curve $c$ on a maximal interval $(\alpha,\beta)$ with $\beta<\infty$ has unit speed, so it is
Cauchy as $t\to\beta$. It converges in the complete space $\mathcal S^\circ$ to a point of the same fibre, because
fibres are closed in $\mathcal S^\circ$ (preimages of points). So it extends, a contradiction. *Local triviality.*
Given $s_0$, take a local section $\sigma$ over an interval $U\ni s_0$ (e.g. by flowing any lift of $\partial_s$ for
short times from a point of $\mathcal S^\circ_{s_0}$). Define $\Psi(s,t)=\Phi^{\hat\nu}_t(\sigma(s))$ on
$U\times\mathbb R$. For fixed $s$, $t\mapsto\Psi(s,t)$ is the integral curve of a nonvanishing field on the connected
1-manifold $\mathcal S^\circ_s\cong\mathbb R$. Its image is open and closed, hence everything, and it is injective
because there are no closed orbits on $\mathbb R$. $d\Psi$ is invertible ($\partial_t\Psi$ vertical, $dp(\partial_s
\Psi)=1$). So $\Psi$ is a bijective local diffeomorphism, hence a diffeomorphism onto $p^{-1}(U)$. $\square$

### A.2 Distance in $SO(m)$ for $\langle A,B\rangle=\operatorname{tr}A^\top B$

Geodesics through $I$ are $t\mapsto e^{tY}$. By Hopf–Rinow (compactness), $d(I,R)=\min\{|Y|_F:e^Y=R,\ Y\in
\mathfrak{so}(m)\}$. $Y$ and $R=e^Y$ commute and are normal. On a common invariant plane $Y$ acts as rotation by
$\phi_j$, and $R$ by $\phi_j$ mod $2\pi$. So $|\phi_j|\ge\theta_j$, where $\theta_j\in[0,\pi]$ is the corresponding
angle of $R$, and $|Y|_F^2=2\sum\phi_j^2\ge2\sum\theta_j^2$. Equality holds iff every $|\phi_j|=\theta_j$: the principal
logarithms. These are unique iff no $\theta_j=\pi$. $\square$

### A.3 Proposition E2 (uniform convergence)

Let $m_f(s)=f(s,q_f(s))$. By minimality, $\mathcal E_K(q_k)+k\oint f(q_k)\le\mathcal E_K(q_f)+k\oint m_f$. Since
$f(q_k)\ge m_f$, this gives (i) $\mathcal E_K(q_k)\le\mathcal E_K(q_f)$ and (ii) $\oint(f(s,q_k)-m_f)\le\mathcal E_K
(q_f)/k$. By (i) and (E2), $(q_k)$ lies in a compact set and is bounded in $H^1$, hence equicontinuous. Any uniform
limit point $q_\infty$ is a continuous section with $\oint(f(q_\infty)-m_f)=0$ by (ii) and continuity. So
$f(s,q_\infty(s))=m_f(s)$ for all $s$, i.e. $q_\infty=q_f$ by uniqueness (E1). Hence the whole family converges. $\square$

### A.4 Crossing lemma (Theorem C2)

Let $u_1\ne u_2$ be minimisers of the same class; for $F=S^1$ use lifts with the same $k$. Then $\underline u=\min
(u_1,u_2)$ and $\bar u=\max(u_1,u_2)$ are $H^1$ sections of the same class. A.e. on $\{u_1<u_2\}$ and on
$\{u_1>u_2\}$ the pair of integrands is permuted. On $\{u_1=u_2\}$ the values agree, and the derivatives agree a.e.
(Stampacchia). So $\mathcal E(\underline u)+\mathcal E(\bar u)=\mathcal E(u_1)+\mathcal E(u_2)=2m$, both are $\ge m$,
and hence both are minimisers, smooth by Lemma 4.1. If $u_1(s_0)=u_2(s_0)$, then either $u_1'(s_0)=u_2'(s_0)$, and
uniqueness for the Euler–Lagrange ODE gives $u_1\equiv u_2$, or $u_1'(s_0)\ne u_2'(s_0)$, and $\underline u$ has a
corner at $s_0$, contradicting smoothness. So $u_1-u_2$ never vanishes. For $F=S^1$, apply this to $\tilde u_1$ and
$\tilde u_2+j\ell$ for every $j$. $\square$

### A.5 Hyperbolicity and parity (Theorem C5)

Along $u^\star$ the Jacobi operator is $\mathcal J\phi=-(\mathcal P\phi')'+\mathcal R\phi$ with $\mathcal P=\mathcal L_
{pp}>0$ and $\mathcal R=\mathcal L_{uu}-(\mathcal L_{up})'$, both $L$-periodic. The linearised Euler–Lagrange flow over
one period is the monodromy $M$ of $\mathcal J\phi=0$. It is unimodular (constant Wronskian) and, up to conjugation,
independent of the coordinates used. Hill's oscillation theorem (Magnus and Winkler 1966; Eastham 1973) [K]: the
periodic eigenvalues $\lambda_0<\lambda_3\le\lambda_4<\lambda_7\le\dots$ and the antiperiodic ones
$\lambda_1\le\lambda_2<\lambda_5\le\dots$ interlace as $\lambda_0<\lambda_1\le\lambda_2<\lambda_3\le\lambda_4<\dots$.
The discriminant $\Delta(\lambda)=\operatorname{tr}M(\lambda)$ of $\mathcal J\phi=\lambda\phi$ satisfies $\Delta>2$ on
$(-\infty,\lambda_0)$ and on $(\lambda_{4j+3},\lambda_{4j+4})$, $\Delta<-2$ on $(\lambda_{4j+1},\lambda_{4j+2})$, and
$|\Delta|<2$ on the remaining open gaps. If $0$ is not a periodic eigenvalue and the Morse index $\iota$ (the number of
negative periodic eigenvalues) is even, $\iota=2j$, then $0\in(\lambda_{4j-1},\lambda_{4j})$ (or $(-\infty,\lambda_0)$
for $j=0$), so $\Delta(0)>2$. If $\iota=2j+1$, then $0\in(\lambda_{4j},\lambda_{4j+3})$, where $\Delta<2$. Hence
$\operatorname{sign}(\Delta(0)-2)=(-1)^\iota$. With $\det M=1$ and $\operatorname{tr}M>2$, the multipliers are real,
positive and $\ne1$: hyperbolic. $\square$

### A.6 Proposition 8.2 (barrier potential)

Let $\mathcal E_K+\oint\mathsf B\le C$, $\delta=d(q(s_0))=\min d(q(s))$, and $A=\sqrt{2C/\underline k}$. By
Cauchy–Schwarz, $d(q(s))\le\delta+A|s-s_0|^{1/2}$. Then
$\oint\mathsf B\ge\int_{\delta^2/A^2\le|t|\le t_0}\frac{c\,dt}{(2A|t|^{1/2})^2}=\frac{c}{2A^2}\log\frac{t_0A^2}{\delta^2}$,
which exceeds $C$ unless $\delta\ge\delta(C):=A\sqrt{t_0}\,e^{-A^2C/c}$. The sections therefore stay in the compact set
$\{d\ge\delta(C)\}$, and (B) follows as in Prop. 8.1(c). $\square$

---

## Appendix B. Formula sheet for the implementation ($r=1$, $K=I$)

- Unit kernel: $x(q)=n/\|n\|$, $n_j=(-1)^j\det(J_b\text{ without column }j)$ (paper).
- Fibre curvature: $\kappa(q)=-J_b^+\big(dJ_b(q)[x]\,x\big)$, with $dJ_b[x]$ by a central difference along $x$.
- Euler–Lagrange: $x^\top q''=0$; ODE form $q''=J_b^+(\zeta_\gamma'-dJ_b[q']\,q')$.
- Hill potential: $\mathcal R=|x'|^2-\kappa\cdot q''$; discrete Hessian in 9.2; Floquet test in 9.3.
- Energy split: $U=\zeta_\gamma^\top(J_bJ_b^\top)^{-1}\zeta_\gamma$; $\mathcal E=\frac12\oint U+\frac12\oint(x\cdot q')^2$.
- Prop. E corrector: $\delta=\frac1k\,\frac{x\cdot q_f''}{1+(q_f-q_r)\cdot\kappa}$.
- Barrier: $\varphi$ of Prop. 8.1, or $\mathsf B=\sum_ic_i\big[(q_i-q_i^-)^{-2}+(q_i^+-q_i)^{-2}\big]+c_\Sigma
  \det(J_bJ_b^\top)^{-1}$. This is $\ge c/d^2$: near a face, the $i$-th term is $\ge c_i/d^2$; near $\Sigma$,
  $\det(J_bJ_b^\top)\le(C\ell_Jd)^2$ by the bound in the proof of 8.1(a).

---

## References

- Ambrosetti, A., Rabinowitz, P. H. (1973). Dual variational methods in critical point theory and applications.
  *J. Funct. Anal.* 14, 349–381.
- Angenent, S. (1988). The periodic orbits of an area preserving twist map. *Comm. Math. Phys.* 115, 353–374.
- Aubry, S., Le Daeron, P. Y. (1983). The discrete Frenkel–Kontorova model and its extensions I. *Physica D* 8, 381–422.
- Baker, D. R., Wampler, C. W. (1988). On the inverse kinematics of redundant manipulators. *Int. J. Robotics Research*
  7(2), 3–21.
- Bangert, V. (1988). Mather sets for twist maps and geodesics on tori. *Dynamics Reported* 1, 1–56.
- Bergou, M., Wardetzky, M., Robinson, S., Audoly, B., Grinspun, E. (2008). Discrete elastic rods. *ACM Trans. Graphics*
  27(3).
- Bishop, R. L. (1975). There is more than one way to frame a curve. *Amer. Math. Monthly* 82(3), 246–251.
- Burdick, J. W. (1989). On the inverse kinematics of redundant manipulators: characterization of the self-motion
  manifolds. *Proc. IEEE ICRA*, 264–270.
- Buttazzo, G., Giaquinta, M., Hildebrandt, S. (1998). *One-dimensional Variational Problems*. Oxford.
- Chiu, S. L. (1988). Task compatibility of manipulator postures. *Int. J. Robotics Research* 7(5), 13–21.
- Eastham, M. S. P. (1973). *The Spectral Theory of Periodic Differential Equations*. Scottish Academic Press.
- Eells, J., Sampson, J. H. (1964). Harmonic mappings of Riemannian manifolds. *Amer. J. Math.* 86, 109–160.
- Ehresmann, C. (1951). Les connexions infinitésimales dans un espace fibré différentiable. *Colloque de Topologie
  (Espaces Fibrés), Bruxelles 1950*, 29–55.
- Hermann, R. (1960). A sufficient condition that a mapping of Riemannian manifolds be a fibre bundle. *Proc. AMS* 11,
  236–242.
- Hsu, P., Hauser, J., Sastry, S. (1989). Dynamic control of redundant manipulators. *J. Robotic Systems* 6(2),
  133–148.
- Jiang, B. (1983). *Lectures on Nielsen Fixed Point Theory*. Contemp. Math. 14, AMS.
- Katok, A. (1982). Some remarks on Birkhoff and Mather twist map theorems. *Ergodic Theory Dynam. Systems* 2, 185–194.
- Kazerounian, K., Wang, Z. (1988). Global versus local optimization in redundancy resolution of robotic manipulators.
  *Int. J. Robotics Research* 7(5), 3–12.
- Khatib, O. (1987). A unified approach for motion and force control of robot manipulators: the operational space
  formulation. *IEEE J. Robotics and Automation* 3(1), 43–53.
- Klein, C. A., Huang, C.-H. (1983). Review of pseudoinverse control for use with kinematically redundant manipulators.
  *IEEE Trans. Systems, Man, Cybernetics* 13(3), 245–250.
- Kolář, I., Michor, P. W., Slovák, J. (1993). *Natural Operations in Differential Geometry*. Springer.
- MacKay, R. S., Meiss, J. D. (1983). Linear stability of periodic orbits in Lagrangian systems. *Phys. Lett. A* 98,
  92–94.
- Magnus, W., Winkler, S. (1966). *Hill's Equation*. Interscience/Wiley.
- Mañé, R. (1996). Generic properties and problems of minimizing measures of Lagrangian systems. *Nonlinearity* 9,
  273–310.
- Marsden, J. E., West, M. (2001). Discrete mechanics and variational integrators. *Acta Numerica* 10, 357–514.
- Mather, J. N. (1982). Existence of quasi-periodic orbits for twist homeomorphisms of the annulus. *Topology* 21,
  457–467.
- Mather, J. N. (1991). Action minimizing invariant measures for positive definite Lagrangian systems. *Math. Z.* 207,
  169–207.
- Moser, J. (1986). Monotone twist mappings and the calculus of variations. *Ergodic Theory Dynam. Systems* 6, 401–413.
- Nakamura, Y., Hanafusa, H. (1987). Optimal redundancy control of robot manipulators. *Int. J. Robotics Research*
  6(1), 32–42.
- Nomizu, K., Ozeki, H. (1961). The existence of complete Riemannian metrics. *Proc. AMS* 12, 889–891.
- O'Neill, B. (1966). The fundamental equations of a submersion. *Michigan Math. J.* 13, 459–469.
- Smith, H. L. (1995). *Monotone Dynamical Systems*. Math. Surveys Monogr. 41, AMS.
- Shamir, T., Yomdin, Y. (1988). Repeatability of redundant manipulators: mathematical solution of the problem.
  *IEEE Trans. Automatic Control* 33(11), 1004–1009.
- Simon, L. (1983). Asymptotics for a class of non-linear evolution equations, with applications to geometric
  problems. *Annals of Math.* 118, 525–571.
- Vilms, J. (1970). Totally geodesic maps. *J. Differential Geom.* 4, 73–79.
- Wang, W., Jüttler, B., Zheng, D., Liu, Y. (2008). Computation of rotation minimizing frames. *ACM Trans. Graphics*
  27(1).
- Wood, C. M. (1987). Harmonic sections and Yang–Mills fields. *Proc. London Math. Soc.* (3) 54, 544–558.
