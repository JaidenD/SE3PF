// The paper's Panda example, end to end (design doc, Panda simulation): the closed loop at the joint-acceleration
// level (eq:fb cancels the dynamics exactly, so qdd = qdd*), 1 ms period, nu_d = 0.08. Two runs:
//   paper:       at rest at Psi(0.1, xi0) with beta = 0, 20 s (the paper's Panda results);
//   self-motion: on the path at gamma(0.1) with beta = 0.2, moving at nu_d, one loop (checks the beta channel).
// Writes one CSV log per run and prints the numbers of the Closed loops table next to their expected values.
//   Usage:  panda_sim [output directory, default out/]
//   Output: <dir>/panda_paper.csv, <dir>/panda_self_motion.csv, and <dir>/panda_limits.csv (the joint limits the
//           summary uses, so the plot needs no copy of them; apps/panda_limits.hpp); plot them with
//           tools/plot_panda.py.
//   Exit status: 0, or 2 if a run failed or left the joint range (the simulation itself has no joint limits).
#include "panda_limits.hpp"
#include "se3pf/panda_example.hpp"

#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using namespace se3pf;

namespace {

    // Columns: t, eta, eta_dot, ref, e_tan, xi1..xi5, xi_norm, xi_design, beta, beta_design, q1..q7, qd1..qd7,
    // qdd1..qdd7, sigma_min, joint_margin, newton (joints numbered from 1, as on the robot).
    bool writeCsv(const std::filesystem::path& file, const PandaRun& run)
    {
        std::ofstream f(file);
        if (!f)
            return false;
        f << "t,eta,eta_dot,ref,e_tan,xi1,xi2,xi3,xi4,xi5,xi_norm,xi_design,beta,beta_design";
        for (const char* name : {"q", "qd", "qdd"})
            for (int i = 1; i <= 7; ++i)
                f << ',' << name << i;
        f << ",sigma_min,joint_margin,newton\n";
        f.precision(10);
        for (const PandaLogRow& r : run.log) {
            f << r.t << ',' << r.eta << ',' << r.etaDot << ',' << r.ref << ',' << r.eTan;
            for (int i = 0; i < 5; ++i)
                f << ',' << r.xi[i];
            f << ',' << r.xiNorm << ',' << r.xiDesign << ',' << r.beta << ',' << r.betaDesign;
            for (const Vec7* v : {&r.q, &r.qd, &r.qdd})
                for (int i = 0; i < 7; ++i)
                    f << ',' << (*v)[i];
            f << ',' << r.sigmaMin << ',' << r.jointMargin << ',' << r.iterations << '\n';
        }
        return static_cast<bool>(f);
    }

    // A run is a valid result for the arm if every step succeeded and the arm stayed inside its joint range,
    // which the simulation does not enforce (PandaRun::ok covers the controller only).
    bool valid(const PandaRun& run)
    {
        return run.ok && run.minJointMargin >= 0.0;
    }

    void printHeader(const char* name, const PandaRun& run)
    {
        std::printf("%s run (nu_d = %g, dt = %g s, %zu rows, last t = %.3f s)%s\n", name, PandaExample::nuD, run.dt,
                    run.log.size(), run.log.empty() ? 0.0 : run.log.back().t, run.ok ? "" : "  FAILED");
        if (!run.ok)
            std::printf("  step failed at t = %.3f s (the pose left the tube, or A became singular)\n",
                        run.failureTime);
    }

    // Ends a summary line with the expected value in brackets, if there is one.
    void endLine(const char* expected)
    {
        if (expected[0] != '\0')
            std::printf("   %s", expected);
        std::printf("\n");
    }

    // Expected values of the arm's numbers, empty where the design doc gives none
    struct ArmExpected {
        const char* margin = "";
        const char* sigma = "";
        const char* velocity = "";
        const char* acceleration = "";
    };

    // The arm's numbers, for both runs (joints numbered from 1, as on the robot)
    void printArm(const PandaRun& run, const ArmExpected& expected)
    {
        std::printf("  smallest joint margin       %.4f rad, joint %d, t = %.3f s%s", run.minJointMargin,
                    run.minMarginJoint + 1, run.minMarginTime,
                    run.minJointMargin >= 0.0 ? "" : "  OUTSIDE THE JOINT RANGE");
        endLine(expected.margin);
        std::printf("  smallest sigma_min(A)       %.4f at t = %.3f s", run.minSigma, run.minSigmaTime);
        endLine(expected.sigma);
        std::printf("  largest |qd_i| / qd_max     %.1f %%, joint %d, t = %.3f s", 100.0 * run.maxVelocityFraction,
                    run.maxVelocityJoint + 1, run.maxVelocityTime);
        endLine(expected.velocity);
        std::printf("  largest |qdd*_i| / qdd_max  %.1f %%, joint %d, t = %.3f s", 100.0 * run.maxAccelerationFraction,
                    run.maxAccelerationJoint + 1, run.maxAccelerationTime);
        endLine(expected.acceleration);
        std::printf("  Newton updates per cycle    at most %d", run.maxIterations);
        endLine("[1]");
    }

}  // namespace

int main(int argc, char** argv)
{
    const std::filesystem::path dir = argc > 1 ? argv[1] : "out/";
    try {
        std::filesystem::create_directories(dir);

        const PandaExample example;  // tube, lift, kernel frame and the two starts (offline)
        std::printf("Panda example: coning circle, l = 0.1 m, L = %.8f; rest-posture lift (closure %.1e), "
                    "beta = X^T (q - q*)\n",
                    example.tube().length(), example.lift().diagnostics().closureGap);

        const PandaRun paper = simulate(example, PandaExample::RunKind::Paper);
        const PandaRun self = simulate(example, PandaExample::RunKind::SelfMotion);

        // Expected values: the design doc (Testing, Closed loops), the paper's in brackets after "paper"
        printHeader("paper", paper);
        std::printf("  start: |xi0| = %.5f; %.2f deg and %.2f mm from gamma(0.1), at rest   "
                    "[0.03122; 15.07 deg, 16.82 mm; paper 15.1 deg, 16.8 mm]\n",
                    paper.startXiNorm, paper.startAttitudeDeg, 1000.0 * paper.startPosition);
        std::printf("  max ||xi| - design|         %.3e at t = %.3f s         [2.43e-5, paper 2.4e-5]\n",
                    paper.maxXiDeviation, paper.maxXiDeviationTime);
        std::printf("  max |beta|                  %.2e at t = %.3f s          [1.5e-5]\n", paper.maxBeta,
                    paper.maxBetaTime);
        std::printf("  final |xi|                  %.2e                        [7.1e-7]\n", paper.finalXiNorm);
        std::printf("  final eta_dot               %.7f                       [0.0800025]\n", paper.finalEtaDot);
        printArm(paper, {"[0.3592 rad, joint 4]", "[0.0474]", "[26 %]", "[49 % at t = 0, without the hardware ramp]"});

        printHeader("self-motion", self);
        std::printf("  start: on the path, |xi0| = %.1e, beta0 = %.4f, eta_dot = %.4f   [0, 0.2, 0.08]\n",
                    self.startXiNorm, self.startBeta, self.log.empty() ? 0.0 : self.log.front().etaDot);
        std::printf("  max |xi|                    %.2e                        [1.2e-6]\n", self.maxXiNorm);
        std::printf("  max |beta - design|         %.3e at t = %.3f s         [1.51e-4, accept < 2e-4]\n",
                    self.maxBetaDeviation, self.maxBetaDeviationTime);
        std::printf("  at the end: max_i |q_i - q*_i(eta)| %.2e, 2-norm %.2e   [9.9e-6, 1.3e-5; accept < 2e-5]\n",
                    self.finalLiftGap, self.finalLiftGapNorm);
        printArm(self, {});

        const std::filesystem::path paperFile = dir / "panda_paper.csv";
        const std::filesystem::path selfFile = dir / "panda_self_motion.csv";
        const std::filesystem::path limitsFile = dir / "panda_limits.csv";
        if (!writeCsv(paperFile, paper) || !writeCsv(selfFile, self) || !writePandaLimits(limitsFile)) {
            std::cerr << "error: cannot write the logs to " << dir << "\n";
            return 1;
        }
        std::printf("wrote %s, %s and %s\n", paperFile.string().c_str(), selfFile.string().c_str(),
                    limitsFile.string().c_str());
        return valid(paper) && valid(self) ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
