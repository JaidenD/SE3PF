// Oracle for the Panda closed loops (the expected values of tests/test_panda_example.cpp and apps/panda_sim.cpp): the
// paper's two Panda runs computed with the reference library pf (~/Desktop/panda_path_following/include), never
// modified. Not part of the CMake build: it needs pf, which is not in this repository.
// pf has no rest-posture lift, so this program adds the paper's lift condition h(q) = x(q)^T (q - q_r) (RestPosture
// below). Settings as in se3pf's PandaExample: W = diag(0.01 I, I), the coning circle (p0 = (0.40, 0, 0.45), r =
// 0.12, alpha = 0.25, R0 = flange attitude at q_r), pf's default Gains (K_p = K_p,beta = 16, K_d = K_d,beta = 8,
// k_t = 5), nu_d = 0.08, no ramp on mu (ramp time 0), dt = 1 ms, and the held-acceleration update
// q += dt qd + dt^2/2 qdd*, qd += dt qdd*. pf_sim's panda demos differ (other gains, a ramp on mu, other lifts),
// so they do not reproduce these numbers.
// Prints: L and the lift's closure and margin; the paper run (start, max ||xi| - design|, max |beta|, final |xi|
// and eta_dot, joint margin, sigma_min(A), velocity and acceleration fractions, Newton updates); the self-motion
// run (max |xi|, max |beta - design|, max_i |q_i - q*_i(eta)| and its 2-norm at the last row); the step time.
// Build and run (a few seconds):
//   g++ -O2 -std=c++17 -I ~/Desktop/panda_path_following/include tools/oracle/panda_oracle.cpp -o /tmp/panda_oracle
//   /tmp/panda_oracle
#include <chrono>
#include <cstdio>
#include <vector>

#include "pf/demos.hpp"

using namespace pf;

class RestPosture final : public LiftCondition {
 public:
  explicit RestPosture(const Vec7& qr) : qr_(qr) {}
  double value(const Vec7& q) const override { return dot(panda::kernelVector(panda::bodyJacobian(q)), q - qr_); }
  Vec7 gradient(const Vec7& q) const override {
    const Vec7 x = panda::kernelVector(panda::bodyJacobian(q)), e = q - qr_;
    Vec7 g;
    constexpr double d = 1e-6;
    for (std::size_t j = 0; j < 7; ++j) {
      Vec7 u;
      u[j] = d;
      const Vec7 dx =
          (panda::kernelVector(panda::bodyJacobian(q + u)) - panda::kernelVector(panda::bodyJacobian(q - u))) / (2.0 * d);
      g[j] = dot(dx, e) + x[j];
    }
    return g;
  }

 private:
  Vec7 qr_;
};

int main() {
  const Metric metric = Metric::product(0.01, 1.0);
  const Mat3 R0 = rotation(panda::forwardKinematics(panda::readyPose()));
  const ConingCirclePath path(ConingCircle(R0, Vec3{0.40, 0.0, 0.45}, 0.12, 0.25), metric);
  const Tube tube(path, metric);
  const RestPosture h(panda::readyPose());
  const Lift lift(path, h, panda::readyPose());
  const KernelFrame frame(lift);
  const KernelFrameOutput beta(lift, frame);
  const double L = path.length();
  std::printf("L = %.8f, lift closure %.1e, margin %.4f, min sigma(Jb) %.4f\n", L, lift.closureGap(),
              lift.minJointMargin(), lift.minSigma());
  {
    const Vec7 q0 = lift.q(0.0);
    std::printf("q*(0) =");
    for (std::size_t i = 0; i < 7; ++i) std::printf(" %.5f", q0[i]);
    std::printf("\n");
  }

  TangentialReference ref;
  ref.velocity = [](double, double, double& v, double& dv, double& dve) {
    v = 0.08;
    dv = 0.0;
    dve = 0.0;
  };
  const double dt = 1e-3;
  const Vec5 xi0{0.02, -0.015, 0.01, -0.005, 0.015};

  // Paper run
  {
    ArmController ctrl(tube, Gains{}, 0.0, &beta);
    ctrl.setReference(ref);
    const Mat4 g0 = tube.Psi(0.1, xi0);
    Vec7 q0;
    if (!startConfiguration(beta, lift, g0, 0.1, 0.0, q0)) return 1;
    std::printf("paper start q0 =");
    for (std::size_t i = 0; i < 7; ++i) std::printf(" %.8f", q0[i]);
    std::printf("\n");
    const Vec6 e0 = productLog(invSE3(path.pose(0.1)) * panda::forwardKinematics(q0));
    std::printf("start: |xi0| = %.5f, %.3f deg, %.3f mm\n", xi0.norm(), head(e0).norm() * 180.0 / kPi,
                1000.0 * tail(e0).norm());
    ctrl.reset(0.1);
    Vec7 q = q0, dq;
    double dev = 0.0, devT = 0.0, bmax = 0.0, margin = 1e9, sig = 1e9, dqf = 0.0, ddqf = 0.0, devEarly = 0.0;
    int marginJoint = -1, dqJoint = -1, ddqJoint = -1, newton = 0;
    double marginT = 0, sigT = 0, dqT = 0, ddqT = 0, bmaxT = 0;
    ControlStatus last;
    const long steps = std::lround(20.0 / dt);
    for (long k = 0; k <= steps; ++k) {
      const double t = k * dt;
      Vec7 qdd;
      ControlStatus st;
      if (!ctrl.step(q, dq, t, dt, qdd, st)) {
        std::printf("FAILED at %g\n", t);
        return 2;
      }
      const double d = std::abs(st.xi.norm() - xi0.norm() * (1 + 4 * t) * std::exp(-4 * t));
      if (d > dev) { dev = d; devT = t; }
      if (t <= 2.0 + 1e-12 && d > devEarly) devEarly = d;
      if (std::abs(st.beta) > bmax) { bmax = std::abs(st.beta); bmaxT = t; }
      for (std::size_t i = 0; i < 7; ++i) {
        const double m = std::min(q[i] - panda::kQMin[i], panda::kQMax[i] - q[i]);
        if (m < margin) { margin = m; marginJoint = int(i); marginT = t; }
        if (std::abs(dq[i]) / panda::kDQMax[i] > dqf) { dqf = std::abs(dq[i]) / panda::kDQMax[i]; dqJoint = int(i); dqT = t; }
        if (std::abs(qdd[i]) / panda::kDDQMax[i] > ddqf) { ddqf = std::abs(qdd[i]) / panda::kDDQMax[i]; ddqJoint = int(i); ddqT = t; }
      }
      if (st.sigma_min < sig) { sig = st.sigma_min; sigT = t; }
      newton = std::max(newton, st.iterations);
      last = st;
      q += dt * dq + (0.5 * dt * dt) * qdd;
      dq += dt * qdd;
    }
    std::printf("paper run: max dev %.6e at t=%.3f (first 2 s: %.6e); max|beta| %.6e at t=%.3f; final |xi| %.6e; final eta_dot %.10f; final beta %.3e\n",
                dev, devT, devEarly, bmax, bmaxT, last.xi.norm(), last.eta_dot, last.beta);
    std::printf("  margin %.6f joint %d at t=%.3f (eta %.4f); min sigma(A) %.6f at t=%.3f; max dq frac %.6f joint %d at t=%.3f; max ddq frac %.6f joint %d at t=%.3f; Newton <= %d\n",
                margin, marginJoint + 1, marginT, 0.0, sig, sigT, dqf, dqJoint + 1, dqT, ddqf, ddqJoint + 1, ddqT, newton);
  }

  // Self-motion run
  {
    ArmController ctrl(tube, Gains{}, 0.0, &beta);
    ctrl.setReference(ref);
    const double s0 = 0.1, beta0 = 0.2, nu = 0.08;
    Vec7 q0;
    if (!startConfiguration(beta, lift, path.pose(s0), s0, beta0, q0)) return 3;
    OutputState<7> o;
    ctrl.outputs(q0, s0, o);
    Vec7 ydot;
    ydot[0] = nu;
    const Vec7 dq0 = LU<7>(o.A).solve(ydot);
    std::printf("self-motion start q0 =");
    for (std::size_t i = 0; i < 7; ++i) std::printf(" %.8f", q0[i]);
    std::printf("\n  dq0 =");
    for (std::size_t i = 0; i < 7; ++i) std::printf(" %.8f", dq0[i]);
    std::printf("\n");
    ctrl.reset(s0);
    Vec7 q = q0, dq = dq0;
    const double T = L / nu;
    double dxi = 0.0, dbeta = 0.0, dbetaT = 0.0, margin = 1e9, sig = 1e9, dqf = 0;
    int k = 0;
    ControlStatus last;
    Vec7 qLast;
    for (; k * dt < T; ++k) {
      const double t = k * dt;
      Vec7 qdd;
      ControlStatus st;
      if (!ctrl.step(q, dq, t, dt, qdd, st)) { std::printf("FAILED\n"); return 4; }
      dxi = std::max(dxi, st.xi.norm());
      const double db = std::abs(st.beta - beta0 * (1.0 + 4.0 * t) * std::exp(-4.0 * t));
      if (db > dbeta) { dbeta = db; dbetaT = t; }
      margin = std::min(margin, panda::jointMargin(q));
      sig = std::min(sig, st.sigma_min);
      for (std::size_t i = 0; i < 7; ++i) dqf = std::max(dqf, std::abs(dq[i]) / panda::kDQMax[i]);
      last = st;
      qLast = q;
      q += dt * dq + (0.5 * dt * dt) * qdd;
      dq += dt * qdd;
    }
    std::printf("self-motion: %d steps; max|xi| %.6e; max|beta - design| %.6e at t=%.3f; margin %.4f; sigma %.4f; dq frac %.4f\n",
                k, dxi, dbeta, dbetaT, margin, sig, dqf);
    // At the last logged row (q before the final update, eta from the step)
    {
      const Vec7 qs = lift.q(last.eta);
      const Vec7 e = qLast - qs;
      std::printf("  at last row: eta %.6f, max|q - q*| %.6e, 2-norm %.6e, beta %.3e, |xi| %.3e\n", last.eta, e.maxAbs(),
                  e.norm(), last.beta, last.xi.norm());
    }
    // After the final update
    {
      ctrl.outputs(q, s0, o);
      const Vec7 qs = lift.q(o.eta);
      const Vec7 e = q - qs;
      std::printf("  after final update: eta %.6f (eta - s0 = %.3e), max|q - q*| %.6e, 2-norm %.6e, beta %.3e, |xi| %.3e\n",
                  o.eta, wrapSymmetric(o.eta - s0, L), e.maxAbs(), e.norm(), o.beta, o.xi.norm());
    }
  }

  // Step time (paper start, moving)
  {
    ArmController ctrl(tube, Gains{}, 0.0, &beta);
    ctrl.setReference(ref);
    Vec7 q0;
    startConfiguration(beta, lift, tube.Psi(0.1, xi0), 0.1, 0.0, q0);
    Vec7 dq{0.1, -0.05, 0.08, 0.1, -0.1, 0.05, 0.12};
    Vec7 qdd;
    ControlStatus st;
    const int n = 20000;
    double sink = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < n; ++k) {
      ctrl.reset(0.1);
      ctrl.step(q0, dq, 0.0, 1e-3, qdd, st);
      sink += qdd[k % 7];
    }
    const double us = 1e6 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / n;
    std::printf("step time %.2f us (sink %g)\n", us, sink);
  }
  return 0;
}
