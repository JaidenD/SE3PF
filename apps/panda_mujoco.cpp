// The Panda example at the torque level in MuJoCo (docs/MUJOCO.md): the last check before hardware.
//
// The controller (Controller<PandaArm>, the paper's gains and the rest-posture lift) drives the MuJoCo Menagerie
// Panda with joint torques. Runs:
//   (0) kinematics: PandaArm::output against MuJoCo at 1000 random q, and at every step of run A;
//   (A) exact model: tau = M qdd* + C qdot + G at every MuJoCo step, so qacc = qdd*; compared with the same
//       controller at the joint-acceleration level (the design doc's numbers), and the difference's O(h) scaling;
//  (A2) the self-motion run with the exact model;
//   (B) zero-order hold: tau from the state at the start of each 1 ms period, held (what the robot does);
//   (C) hardware-like: hold, mu ramp T = 1 s, joint damping on and not compensated; every margin of the hardware
//       section; (C0) the same without the ramp, and C without the damping (to tell their effects apart, and to
//       check that the damping is felt);
//   (D) model mismatch: the plant carries the Franka Hand (0.76 kg), the torque uses the arm alone; (D-G) the same
//       with the hand's weight in G, as libfranka adds it when the hand is configured.
// Each run writes out/mujoco_<run>.csv (and out/mujoco_path.csv samples the path gamma, for the video, and
// out/panda_limits.csv holds the joint limits, for the plot); the summary at the end lists every check with its
// expected value. The hardware bounds are HardwareBounds (se3pf/panda.hpp).
//   Usage:  panda_mujoco [output directory, default out/] [Menagerie franka_emika_panda directory]
//   Plot:   python3 tools/plot_mujoco.py   Video: tools/render_mujoco.py   Live replay: tools/view_mujoco.py
#include "mujoco_panda.hpp"
#include "panda_limits.hpp"
#include "se3pf/panda.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace se3pf;
using Start = MujocoRunOptions::Start;
using Torque = MujocoRunOptions::Torque;

namespace {

    // ---------------------------------------------------------------------------------------------------------
    // The summary table: one row per check, with the expected value or bound in brackets

    enum class Status { Pass, Fail, Info };

    struct Check {
        std::string name;
        std::string value;
        std::string expected;
        Status status = Status::Info;
    };

    std::vector<Check> gChecks;

    std::string fmt(const char* format, double x)
    {
        char buf[64];
        std::snprintf(buf, sizeof(buf), format, x);
        return buf;
    }

    void check(const std::string& name, const std::string& value, const std::string& expected, bool pass)
    {
        gChecks.push_back({name, value, expected, pass ? Status::Pass : Status::Fail});
    }

    void info(const std::string& name, const std::string& value, const std::string& expected = "")
    {
        gChecks.push_back({name, value, expected, Status::Info});
    }

    void heading(const std::string& title) { gChecks.push_back({title, "", "", Status::Info}); }

    std::string pct(double fraction) { return fmt("%.1f %%", 100.0 * fraction); }

    // ---------------------------------------------------------------------------------------------------------
    // Logs

    bool writeCsv(const std::filesystem::path& file, const MujocoRun& run)
    {
        std::FILE* f = std::fopen(file.string().c_str(), "w");
        if (!f)
            return false;
        std::fprintf(f, "t,eta,eta_dot,xi1,xi2,xi3,xi4,xi5,xi_norm,xi_design,beta,beta_design");
        for (const char* name : {"q", "qd", "tau", "grav", "qdd"})
            for (int i = 1; i <= 7; ++i)
                std::fprintf(f, ",%s%d", name, i);
        std::fprintf(f, ",sigma_min,joint_margin,flange_speed,tau_rate_cmd,tau_rate,ncon,nefc,ramp,"
                        "xi_norm_kin,beta_kin,q_err_kin\n");
        for (const MujocoLogRow& r : run.log) {
            std::fprintf(f, "%.6g,%.12g,%.12g", r.t, r.eta, r.etaDot);
            for (int i = 0; i < 5; ++i)
                std::fprintf(f, ",%.9g", r.xi[i]);
            std::fprintf(f, ",%.9g,%.9g,%.9g,%.9g", r.xiNorm, r.xiDesign, r.beta, r.betaDesign);
            for (const Vec7* v : {&r.q, &r.qd, &r.tau, &r.G, &r.qdd})
                for (int i = 0; i < 7; ++i)
                    std::fprintf(f, ",%.12g", (*v)[i]);
            std::fprintf(f, ",%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%.9g,%.9g,%.9g,%.9g\n", r.sigmaMin, r.jointMargin,
                         r.flangeSpeed, r.tauRateCmd, r.tauRate, r.contacts, r.constraintRows, r.ramp, r.xiNormKin,
                         r.betaKin, r.qErrKin);
        }
        const bool ok = std::ferror(f) == 0;
        return std::fclose(f) == 0 && ok;
    }

    // The path gamma(s) at n samples s_k = k L / n, for the video and the live replay (tools/render_mujoco.py):
    // columns s, p (3), R row by row (9).
    bool writePathCsv(const std::filesystem::path& file, const PandaExample& ex, int n)
    {
        std::FILE* f = std::fopen(file.string().c_str(), "w");
        if (!f)
            return false;
        std::fprintf(f, "s,x,y,z,r11,r12,r13,r21,r22,r23,r31,r32,r33\n");
        const double L = ex.path().length();
        for (int k = 0; k < n; ++k) {
            const double s = k * L / n;
            const Pose g = ex.path().pose(s);
            std::fprintf(f, "%.12g,%.12g,%.12g,%.12g", s, g(0, 3), g(1, 3), g(2, 3));
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    std::fprintf(f, ",%.12g", g(r, c));
            std::fprintf(f, "\n");
        }
        const bool ok = std::ferror(f) == 0;
        return std::fclose(f) == 0 && ok;
    }

    // Runs, prints one line, and writes the log (if a file name is given).
    MujocoRun runAndLog(const PandaExample& ex, const char* name, const MujocoRunOptions& o,
                        const std::filesystem::path& dir, const char* csv)
    {
        MujocoRun run = simulateMujoco(ex, o);
        std::printf("  %-26s %6ld periods, %5.2f s%s\n", name, run.steps, run.wallTime,
                    run.ok ? "" : "  (FAILED)");
        if (!run.ok && run.failureTime >= 0.0)
            std::printf("    step failed at t = %.3f s\n", run.failureTime);
        if (csv) {
            const std::filesystem::path file = dir / csv;
            if (!writeCsv(file, run))
                throw std::runtime_error("cannot write " + file.string());
        }
        return run;
    }

    // ---------------------------------------------------------------------------------------------------------
    // The checks of each run

    // (0) PandaArm::output against MuJoCo at random configurations within the joint limits
    void kinematicsCheck(MujocoPanda& arm, int samples)
    {
        std::mt19937 gen(2026);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        const Vec7 lo = PandaArm::qMin(), hi = PandaArm::qMax();
        double poseErr = 0.0, jacErr = 0.0;
        for (int k = 0; k < samples; ++k) {
            Vec7 q = Vec7::Zero();
            for (int i = 0; i < 7; ++i)
                q[i] = lo[i] + (hi[i] - lo[i]) * u(gen);
            Mat67 Jb = Mat67::Zero(), JbMj = Mat67::Zero();
            const Pose g = PandaArm::output(q, &Jb);
            const Pose gMj = arm.flange(q, &JbMj);
            poseErr = std::max(poseErr, (g - gMj).cwiseAbs().maxCoeff());
            jacErr = std::max(jacErr, (Jb - JbMj).cwiseAbs().maxCoeff());
        }
        heading("(0) Kinematics: PandaArm::output against MuJoCo (flange = link7 + 0.107 m along z)");
        check("  FK, 1000 random q: max |g - g_mj| (entries)", fmt("%.1e", poseErr), "[~1e-15; < 1e-10]",
              poseErr < 1e-10);
        check("  J_b, 1000 random q: max |J_b - J_b,mj|", fmt("%.1e", jacErr), "[~5e-15; < 1e-10]", jacErr < 1e-10);
    }

    // " (joint j)": the joint with the largest entry of a per-joint peak vector, one-based as on the robot
    std::string peakJoint(const Vec7& perJoint)
    {
        int j = 0;
        perJoint.maxCoeff(&j);
        return " (joint " + std::to_string(j + 1) + ")";
    }

    // Per-joint peaks in percent, "j1 j2 ... j7"
    std::string perJointPct(const Vec7& perJoint)
    {
        std::string s;
        for (int i = 0; i < 7; ++i)
            s += fmt(i == 0 ? "%.1f" : " %.1f", 100.0 * perJoint[i]);
        return s + " %";
    }

    // The torque-rate bound, N m/s: 60 % of PandaArm::tauRateMax
    constexpr double kTauRateBound = HardwareBounds::tauRateFrac * PandaArm::tauRateMax;

    // Peaks every run reports, with the hardware margins as bounds (HardwareBounds; design doc, Hardware)
    void marginChecks(const MujocoRun& r, bool asChecks)
    {
        using B = HardwareBounds;
        auto row = [&](const std::string& name, const std::string& value, const std::string& bound, bool ok) {
            if (asChecks)
                check(name, value, bound, ok);
            else
                info(name, value + (ok ? "" : "  VIOLATED"), bound);
        };
        row("  smallest joint margin", fmt("%.4f rad", r.minJointMargin) + " (joint " +
                std::to_string(r.minMarginJoint + 1) + ")", fmt("[>= %.2f rad]", B::jointMargin),
            r.minJointMargin >= B::jointMargin);
        row("  largest |qd_i| / qd_max", pct(r.maxDqFrac) + peakJoint(r.peakDqFrac),
            fmt("[<= %.0f %%]", 100.0 * B::dqFrac), r.maxDqFrac <= B::dqFrac);
        row("  largest commanded |qdd*_i| / qdd_max", pct(r.maxDdqFrac) + peakJoint(r.peakDdqFrac) +
            fmt(" at t = %.3f s", r.maxDdqFracTime), fmt("[<= %.0f %%]", 100.0 * B::ddqFrac),
            r.maxDdqFrac <= B::ddqFrac);
        row("  largest delivered |tau_i| / tau_max", pct(r.maxTauFrac) + peakJoint(r.peakTauFrac),
            fmt("[<= %.0f %%]", 100.0 * B::tauFrac), r.maxTauFrac <= B::tauFrac);
        row("  torque rate |d(tau - G)/dt| (commanded)", fmt("%.1f N m/s", r.maxTauRateCmd),
            fmt("[<= %.0f N m/s]", kTauRateBound), r.maxTauRateCmd <= kTauRateBound);
        info("  torque rate |d tau/dt| (full, with G)", fmt("%.1f N m/s", r.maxTauRate), "[for comparison]");
        row("  flange speed", fmt("%.3f m/s", r.maxFlangeSpeed), fmt("[<= %.1f m/s]", B::flangeSpeed),
            r.maxFlangeSpeed <= B::flangeSpeed);
        row("  max |xi|", fmt("%.4f", r.maxXi), fmt("[<= %.2f]", B::xi), r.maxXi <= B::xi);
        row("  max |beta|", fmt("%.2e rad", r.maxBeta), fmt("[<= %.2f rad]", B::beta), r.maxBeta <= B::beta);
        static_assert(B::sigmaMin == 1e-3, "update the label below");
        row("  smallest sigma_min(A)", fmt("%.4f", r.minSigma), "[>= 1e-3]", r.minSigma >= B::sigmaMin);
        row("  contacts; torque clamped; MuJoCo warnings",
            std::to_string(r.maxContacts) + "; " + std::to_string(r.clampedSteps) + "; " +
                std::to_string(r.warnings), "[0; 0; 0]",
            r.maxContacts == 0 && r.clampedSteps == 0 && r.warnings == 0);
    }

    bool near(double value, double expected, double tol) { return std::abs(value - expected) <= tol; }

    // The baseline (joint-acceleration level) against the design doc's Panda rows (Testing, Closed loops)
    void baselinePaperChecks(const MujocoRun& r)
    {
        check("  baseline: max ||xi| - design|", fmt("%.3e", r.baselineMaxDeviation), "[2.43e-5; paper 2.4e-5]",
              near(r.baselineMaxDeviation, 2.43e-5, 0.01e-5));
        check("  baseline: max |beta|", fmt("%.2e", r.baselineMaxBeta), "[1.5e-5]",
              near(r.baselineMaxBeta, 1.5e-5, 0.05e-5));
        check("  baseline: final |xi|", fmt("%.2e", r.baselineFinalXi), "[7.1e-7]",
              near(r.baselineFinalXi, 7.1e-7, 0.05e-7));
        check("  baseline: final eta_dot", fmt("%.7f", r.baselineFinalEtaDot), "[0.0800025]",
              near(r.baselineFinalEtaDot, 0.0800025, 5e-8));
    }

    // ---------------------------------------------------------------------------------------------------------

    void printSummary()
    {
        std::printf("\nSummary (expected values and bounds in brackets)\n");
        std::printf("%-58s %-42s %-32s %s\n", "check", "value", "[expected]", "result");
        for (const Check& c : gChecks) {
            if (c.value.empty() && c.expected.empty()) {
                std::printf("\n%s\n", c.name.c_str());
                continue;
            }
            const char* status = c.status == Status::Pass ? "PASS" : c.status == Status::Fail ? "FAIL" : "info";
            std::printf("%-58s %-42s %-32s %s\n", c.name.c_str(), c.value.c_str(), c.expected.c_str(), status);
        }
    }

}  // namespace

int main(int argc, char** argv)
{
    const std::filesystem::path dir = argc > 1 ? argv[1] : "out/";
    const std::filesystem::path models = argc > 2 ? argv[2] : defaultPandaModelDir();
    try {
        std::filesystem::create_directories(dir);
        const std::string noHand = (models / "panda_nohand.xml").string();
        const std::string withHand = (models / "panda.xml").string();

        const PandaExample ex;  // path, tube, lift, kernel frame, starts (offline, well under a second)
        std::printf("Panda example in MuJoCo: coning circle, l = 0.1 m, L = %.8f, nu_d = %.2f, dt = %g s\n",
                    ex.tube().length(), PandaExample::nuD, PandaExample::dt);
        std::printf("models: %s\n\nRuns:\n", models.string().c_str());

        if (!writePathCsv(dir / "mujoco_path.csv", ex, 2000))
            throw std::runtime_error("cannot write " + (dir / "mujoco_path.csv").string());
        if (!writePandaLimits(dir / "panda_limits.csv"))  // tau_max for tools/plot_mujoco.py
            throw std::runtime_error("cannot write " + (dir / "panda_limits.csv").string());

        // (0) Kinematics at random q
        {
            MujocoPanda arm(noHand);
            kinematicsCheck(arm, 1000);
        }

        // (A) The exact model: panda_nohand.xml, damping off, exact computed torque at every MuJoCo step
        MujocoRunOptions a;
        a.plant = noHand;
        a.torque = Torque::Exact;
        a.baseline = true;
        a.crossCheck = true;
        const MujocoRun A = runAndLog(ex, "A  exact, h = 1 ms", a, dir, "mujoco_A.csv");
        MujocoRunOptions aFine = a;
        aFine.keepLog = false;
        aFine.crossCheck = false;
        aFine.substeps = 2;
        const MujocoRun A2h = runAndLog(ex, "A  exact, h = 0.5 ms", aFine, dir, nullptr);
        aFine.substeps = 4;
        const MujocoRun A4h = runAndLog(ex, "A  exact, h = 0.25 ms", aFine, dir, nullptr);

        // (A2) Self-motion, exact model
        MujocoRunOptions sm = a;
        sm.start = Start::SelfMotion;
        sm.crossCheck = false;
        const MujocoRun S = runAndLog(ex, "A2 self-motion, h = 1 ms", sm, dir, "mujoco_A2.csv");
        sm.keepLog = false;
        sm.substeps = 2;
        const MujocoRun S2h = runAndLog(ex, "A2 self-motion, h = 0.5 ms", sm, dir, nullptr);

        // (B) Zero-order hold, RK4, damping off
        MujocoRunOptions b;
        b.plant = noHand;
        b.torque = Torque::Hold;
        b.integrator = MujocoPanda::Integrator::RK4;
        b.baseline = true;
        const MujocoRun B = runAndLog(ex, "B  hold, dt = 1 ms", b, dir, "mujoco_B.csv");
        MujocoRunOptions bHalf = b;
        bHalf.dt = 0.5e-3;
        bHalf.keepLog = false;
        const MujocoRun Bh = runAndLog(ex, "B  hold, dt = 0.5 ms", bHalf, dir, nullptr);

        // (C) Hardware-like: hold, ramp T = 1 s, damping on and not compensated; (C0) without the ramp
        MujocoRunOptions c = b;
        c.damping = true;
        c.rampTime = 1.0;
        const MujocoRun C = runAndLog(ex, "C  hardware-like", c, dir, "mujoco_C.csv");
        MujocoRunOptions c0 = c;
        c0.rampTime = 0.0;
        const MujocoRun C0 = runAndLog(ex, "C0 hardware-like, no ramp", c0, dir, "mujoco_C0.csv");
        MujocoRunOptions cNoDamping = c;  // to tell the damping's effect from the ramp's and the hold's
        cNoDamping.damping = false;
        cNoDamping.keepLog = false;
        const MujocoRun Cnd = runAndLog(ex, "C  without damping", cNoDamping, dir, nullptr);

        // (D) Model mismatch: the plant carries the hand, the torque uses the arm alone; (D-G) G with the hand
        MujocoRunOptions d = c;
        d.plant = withHand;
        d.model = noHand;
        d.baseline = false;
        const MujocoRun D = runAndLog(ex, "D  hand not modeled", d, dir, "mujoco_D.csv");
        MujocoRunOptions dg = d;
        dg.plantGravity = true;
        const MujocoRun DG = runAndLog(ex, "D-G hand only in G", dg, dir, "mujoco_DG.csv");

        // ---------------------------------------------------------------------------------------------------
        // Checks

        check("  FK and J_b along run A (20001 q): max error", fmt("%.1e", std::max(A.maxPoseError,
              A.maxJacobianError)), "[~1e-15; < 1e-10]", A.maxPoseError < 1e-10 && A.maxJacobianError < 1e-10);

        const bool allOk = A.ok && A2h.ok && A4h.ok && S.ok && S2h.ok && B.ok && Bh.ok && C.ok && C0.ok && Cnd.ok &&
                           D.ok && DG.ok;
        heading("(A) Exact model: panda_nohand.xml, damping off, tau = M qdd* + C qdot + G at every step");
        check("  max |qacc - qdd*|", fmt("%.1e", A.maxQaccError), "[round-off; < 1e-10]", A.maxQaccError < 1e-10);
        baselinePaperChecks(A);
        // MuJoCo's semi-implicit Euler moves q by h qdot_{k+1} = h qdot_k + h^2 qdd*, h^2/2 qdd* more per step than
        // the exact update of the baseline: an O(h) difference, checked by halving the MuJoCo step below.
        check("  MuJoCo: max ||xi| - design|", fmt("%.3e", A.maxDeviation) + fmt(" at t = %.2f s",
              A.maxDeviationTime), "[2.43e-5 + O(h) Euler; < 5e-5]", A.maxDeviation < 5e-5);
        check("  MuJoCo: max |beta|", fmt("%.2e", A.maxBeta), "[1.5e-5 + O(h); < 5e-5]", A.maxBeta < 5e-5);
        check("  MuJoCo: final |xi|", fmt("%.2e", A.finalXi), "[7.1e-7 + O(h); < 2e-5]", A.finalXi < 2e-5);
        check("  MuJoCo: final eta_dot", fmt("%.7f", A.finalEtaDot), "[0.0800025; +-1e-5]",
              near(A.finalEtaDot, 0.0800025, 1e-5));
        check("  MuJoCo: smallest joint margin", fmt("%.4f rad", A.minJointMargin) + " (joint " +
              std::to_string(A.minMarginJoint + 1) + ")", "[0.3592, joint 4]",
              near(A.minJointMargin, 0.3592, 1e-3) && A.minMarginJoint == 3);
        check("  MuJoCo: smallest sigma_min(A)", fmt("%.4f", A.minSigma), "[0.0474]", near(A.minSigma, 0.0474, 1e-3));
        check("  MuJoCo: largest |qd_i| / qd_max", pct(A.maxDqFrac), "[26 %]", near(A.maxDqFrac, 0.26, 0.01));
        check("  MuJoCo: Newton updates per cycle, at most", std::to_string(A.maxIterations), "[1]",
              A.maxIterations == 1);
        check("  MuJoCo: contacts; constraint rows", std::to_string(A.maxContacts) + "; " +
              std::to_string(A.maxConstraintRows), "[0; 0]", A.maxContacts == 0 && A.maxConstraintRows == 0);
        {
            const double r1 = A.maxBaselineQError / A2h.maxBaselineQError;
            const double r2 = A2h.maxBaselineQError / A4h.maxBaselineQError;
            check("  max |q_mj - q_kin|, h = 1, 0.5, 0.25 ms", fmt("%.2e, ", A.maxBaselineQError) +
                  fmt("%.2e, ", A2h.maxBaselineQError) + fmt("%.2e rad", A4h.maxBaselineQError),
                  fmt("[O(h): ratios %.2f", r1) + fmt(", %.2f ~ 2]", r2),
                  r1 > 1.8 && r1 < 2.2 && r2 > 1.8 && r2 < 2.2);
            info("  max ||xi| - design| at h = 0.5, 0.25 ms", fmt("%.3e, ", A2h.maxDeviation) +
                 fmt("%.3e", A4h.maxDeviation), "[-> baseline 2.43e-5 as h -> 0]");
        }

        heading("(A2) Self-motion, exact model: gamma(0.1), beta0 = 0.2, A(q0) qd0 = (0.08, 0, ..., 0), one loop");
        check("  baseline: max |xi|", fmt("%.2e", S.baselineMaxXi), "[1.2e-6]", near(S.baselineMaxXi, 1.2e-6, 0.1e-6));
        check("  baseline: max |beta - beta0 (1 + 4t) e^-4t|", fmt("%.3e", S.baselineMaxBetaDeviation),
              "[1.51e-4; < 2e-4]", near(S.baselineMaxBetaDeviation, 1.51e-4, 0.01e-4));
        check("  baseline: final max_i |q_i - q*_i(eta)|; 2-norm", fmt("%.2e; ", S.baselineFinalLiftError) +
              fmt("%.2e", S.baselineFinalLiftError2), "[9.9e-6; 1.3e-5 (< 2e-5)]",
              near(S.baselineFinalLiftError, 9.9e-6, 0.1e-6) && S.baselineFinalLiftError2 < 2e-5);
        check("  MuJoCo: max |xi|", fmt("%.2e", S.maxXi), "[1.2e-6 + O(h); < 5e-5]", S.maxXi < 5e-5);
        check("  MuJoCo: max |beta - beta0 (1 + 4t) e^-4t|", fmt("%.3e", S.maxBetaDeviation),
              "[1.51e-4 + O(h); < 3e-4]", S.maxBetaDeviation < 3e-4);
        check("  MuJoCo: final max_i |q_i - q*_i(eta)|; 2-norm", fmt("%.2e; ", S.finalLiftError) +
              fmt("%.2e", S.finalLiftError2), "[9.9e-6 + O(h); < 1e-4]", S.finalLiftError < 1e-4);
        {
            const double r = S.maxBaselineQError / S2h.maxBaselineQError;
            check("  max |q_mj - q_kin|, h = 1, 0.5 ms", fmt("%.2e, ", S.maxBaselineQError) +
                  fmt("%.2e rad", S2h.maxBaselineQError), fmt("[O(h): ratio %.2f ~ 2]", r), r > 1.8 && r < 2.2);
            info("  at h = 0.5 ms: max |xi|; beta dev.; lift error", fmt("%.2e; ", S2h.maxXi) +
                 fmt("%.2e; ", S2h.maxBetaDeviation) + fmt("%.2e", S2h.finalLiftError), "[all closer to baseline]");
        }
        check("  MuJoCo: contacts; constraint rows", std::to_string(S.maxContacts) + "; " +
              std::to_string(S.maxConstraintRows), "[0; 0]", S.maxContacts == 0 && S.maxConstraintRows == 0);

        heading("(B) Zero-order hold: tau from the period's start, held 1 ms (RK4), damping off");
        check("  max ||xi| - design|", fmt("%.3e", B.maxDeviation), "[O(dt); < 1e-4]", B.maxDeviation < 1e-4);
        check("  max |beta|", fmt("%.2e", B.maxBeta), "[O(dt); < 5e-4]", B.maxBeta < 5e-4);
        check("  final |xi|", fmt("%.2e", B.finalXi), "[7.1e-7 + O(dt); < 1e-4]", B.finalXi < 1e-4);
        check("  final eta_dot", fmt("%.7f", B.finalEtaDot), "[0.08 + O(dt); +-5e-4]", near(B.finalEtaDot, 0.08, 5e-4));
        check("  smallest joint margin; sigma_min(A); largest |qd_i|",
              fmt("%.4f; ", B.minJointMargin) + fmt("%.4f; ", B.minSigma) + pct(B.maxDqFrac),
              "[0.3592; 0.0474; 26 %]", near(B.minJointMargin, 0.3592, 1e-3) && near(B.minSigma, 0.0474, 1e-3) &&
              near(B.maxDqFrac, 0.26, 0.01));
        {
            const double rBeta = B.maxBeta / Bh.maxBeta;
            const double rEta = (B.finalEtaDot - 0.08) / (Bh.finalEtaDot - 0.08);
            check("  max |beta| at dt = 1, 0.5 ms", fmt("%.2e, ", B.maxBeta) + fmt("%.2e", Bh.maxBeta),
                  fmt("[O(dt): ratio %.2f ~ 2]", rBeta), rBeta > 1.8 && rBeta < 2.2);
            check("  eta_dot - nu_d at dt = 1, 0.5 ms", fmt("%.2e, ", B.finalEtaDot - 0.08) +
                  fmt("%.2e", Bh.finalEtaDot - 0.08), fmt("[O(dt): ratio %.2f ~ 2]", rEta), rEta > 1.8 && rEta < 2.2);
        }
        check("  contacts; constraint rows; clamped; warnings", std::to_string(B.maxContacts) + "; " +
              std::to_string(B.maxConstraintRows) + "; " + std::to_string(B.clampedSteps) + "; " +
              std::to_string(B.warnings), "[0; 0; 0; 0]",
              B.maxContacts == 0 && B.maxConstraintRows == 0 && B.clampedSteps == 0 && B.warnings == 0);

        heading("(C) Hardware-like: hold 1 ms, mu ramp T = 1 s, joint damping 1 N m s/rad not compensated, 20 s");
        marginChecks(C, true);
        info("  per joint: peak |tau_i| / tau_max,i, joints 1-7", perJointPct(C.peakTauFrac),
             "[joints 5-7: of 12 N m]");
        // The damping is in the plant and not compensated: it drags eta_dot below nu_d (k_t = 5 has no integral
        // action) and leaves steady errors that the same run without damping does not have. Were the damping off
        // in the plant, or compensated in the torque, C would be the run without damping, and these checks fail.
        check("  damping felt: final eta_dot; without damping", fmt("%.4f; ", C.finalEtaDot) +
              fmt("%.4f", Cnd.finalEtaDot), "[< 0.07; 0.08 +- 5e-4]",
              C.finalEtaDot < 0.07 && near(Cnd.finalEtaDot, 0.08, 5e-4));
        check("  damping felt: last 5 s max |xi|; without damping", fmt("%.2e; ", C.steadyXi) +
              fmt("%.2e", Cnd.steadyXi), "[ratio > 10]", C.steadyXi > 10.0 * Cnd.steadyXi);
        check("  damping felt: last 5 s max |beta|; without damping", fmt("%.2e; ", C.steadyBeta) +
              fmt("%.2e", Cnd.steadyBeta), "[ratio > 10]", C.steadyBeta > 10.0 * Cnd.steadyBeta);
        check("  baseline (no damping) with the ramp: |qdd*|; |qd|", pct(C.baselineMaxDdqFrac) + "; " +
              pct(C.baselineMaxDqFrac), "[11.5 %; 22 % (design doc)]",
              near(C.baselineMaxDdqFrac, 0.115, 0.005) && near(C.baselineMaxDqFrac, 0.22, 0.01));

        heading("(C0) The same without the ramp: what the ramp is for");
        check("  first command: largest |qdd*_i| / qdd_max", pct(C0.maxDdqFrac) + fmt(" at t = %.3f s",
              C0.maxDdqFracTime), "[49 % > 30 %: violates]", near(C0.maxDdqFrac, 0.49, 0.01) &&
              C0.maxDdqFracTime == 0.0);
        check("  torque rate |d(tau - G)/dt| (jump at t = 0)", fmt("%.0f N m/s", C0.maxTauRateCmd),
              fmt("[> %.0f: violates]", kTauRateBound), C0.maxTauRateCmd > kTauRateBound);
        info("  with the ramp (run C)", pct(C.maxDdqFrac) + fmt(", %.1f N m/s", C.maxTauRateCmd), "[within 30 %, 600]");

        heading("(D) Model mismatch: plant panda.xml (hand 0.76 kg, fingers open), torque from the arm alone");
        info("  last 5 s: max |xi|; max |beta|", fmt("%.2e; ", D.steadyXi) + fmt("%.2e", D.steadyBeta),
             "[no integral action]");
        info("  final eta_dot", fmt("%.4f", D.finalEtaDot), "[nu_d = 0.08]");
        marginChecks(D, false);
        check("  contacts; constraint rows (finger coupling, limits)", std::to_string(D.maxContacts) + "; " +
              std::to_string(D.maxConstraintRows), "[0; >= 1]", D.maxContacts == 0 && D.maxConstraintRows >= 1);

        heading("(D-G) The same with the hand's weight in G (libfranka adds G with the configured end effector)");
        marginChecks(DG, true);
        info("  last 5 s: max |xi|; max |beta|", fmt("%.2e; ", DG.steadyXi) + fmt("%.2e", DG.steadyBeta),
             fmt("[run C: %.2e; ", C.steadyXi) + fmt("%.2e]", C.steadyBeta));
        info("  final eta_dot", fmt("%.4f", DG.finalEtaDot), fmt("[run C: %.4f]", C.finalEtaDot));

        heading("All runs (A, A2, B, C, C0, D, D-G and the variants above)");
        check("  completed without a failed step or a MuJoCo warning", allOk ? "yes" : "no", "[yes]", allOk);

        printSummary();
        const bool allPass = std::none_of(gChecks.begin(), gChecks.end(),
                                          [](const Check& c) { return c.status == Status::Fail; });
        std::printf("\n%s. Logs: %s; plot with tools/plot_mujoco.py\n",
                    allPass ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED",
                    (dir / "mujoco_{A,A2,B,C,C0,D,DG}.csv").string().c_str());
        return allPass ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
