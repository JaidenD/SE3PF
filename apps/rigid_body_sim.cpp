// The paper's rigid-body example, end to end (design doc, Rigid-body simulation): runs the closed loop with the
// pose reference and with the velocity reference (20 s each, 2 ms sample-and-hold, 4 plant substeps), writes one
// CSV log per run and prints the numbers of the Closed loops table next to their expected values.
//   Usage:  rigid_body_sim [output directory, default out/]
//   Output: <dir>/rigid_body_pose.csv, <dir>/rigid_body_velocity.csv; plot them with tools/plot_rigid_body.py.
#include "se3pf/rigid_body_example.hpp"

#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using namespace se3pf;

namespace {

    // Columns: t, eta, eta_dot, ref, e_tan, xi1..xi5, xi_norm, xi_design, u1..u6 (u = (tau, f)).
    bool writeCsv(const std::filesystem::path& file, const RigidBodyRun& run)
    {
        std::ofstream f(file);
        if (!f)
            return false;
        f << "t,eta,eta_dot,ref,e_tan,xi1,xi2,xi3,xi4,xi5,xi_norm,xi_design,u1,u2,u3,u4,u5,u6\n";
        f.precision(10);
        for (const RigidBodyLogRow& r : run.log) {
            f << r.t << ',' << r.eta << ',' << r.etaDot << ',' << r.ref << ',' << r.eTan;
            for (int i = 0; i < 5; ++i)
                f << ',' << r.xi[i];
            f << ',' << r.xiNorm << ',' << r.xiDesign;
            for (int i = 0; i < 6; ++i)
                f << ',' << r.u[i];
            f << '\n';
        }
        return static_cast<bool>(f);
    }

    // One run's numbers, with the expected values of the design doc (Testing, Closed loops) in brackets.
    void printSummary(const char* name, const RigidBodyRun& run, double expDeviation, double expXi, double expTan,
                      double expInput)
    {
        std::printf("%s reference (dt = %g s, %d substeps, %zu rows)%s\n", name, run.dt, run.substeps, run.log.size(),
                    run.ok ? "" : "  FAILED");
        if (!run.ok)
            std::printf("  step failed at t = %.3f s (the pose left the tube)\n", run.failureTime);
        std::printf("  max ||xi| - design|  %.3e at t = %.3f s   [%.2e, paper: within 1e-3]\n", run.maxDeviation,
                    run.maxDeviationTime, expDeviation);
        std::printf("  final |xi|           %.2e                  [%.1e]\n", run.finalXiNorm, expXi);
        std::printf("  final e_tan          %.2e                 [%.1e]\n", run.finalTanError, expTan);
        std::printf("  max |u|_inf          %.2f at t = %.3f s      [%.1f at t = 0]\n", run.maxInput,
                    run.maxInputTime, expInput);
        std::printf("  Newton updates per projection: at most %d\n", run.maxIterations);
    }

}  // namespace

int main(int argc, char** argv)
{
    const std::filesystem::path dir = argc > 1 ? argv[1] : "out/";
    try {
        std::filesystem::create_directories(dir);

        const RigidBodyExample example;  // builds the tube (offline)
        std::printf("Rigid-body example: tumbling circle, r = 1 m, W = I, L = %.8f\n", example.tube().length());

        const RigidBodyRun pose = simulate(example, RigidBodyExample::Reference::Position);
        const RigidBodyRun velocity = simulate(example, RigidBodyExample::Reference::Velocity);

        std::printf("start: |xi0| = %.4f; %.2f deg and %.4f m from gamma(0.1)   [0.6245; 25.45 deg, 0.4389 m]\n",
                    RigidBodyExample::xi0().norm(), pose.startAttitudeDeg, pose.startPosition);
        printSummary("pose", pose, 9.69e-4, 1.8e-6, -1.7e-6, 29.1);
        printSummary("velocity", velocity, 9.51e-4, 2.8e-6, -3.2e-6, 14.5);

        const std::filesystem::path poseFile = dir / "rigid_body_pose.csv";
        const std::filesystem::path velocityFile = dir / "rigid_body_velocity.csv";
        if (!writeCsv(poseFile, pose) || !writeCsv(velocityFile, velocity)) {
            std::cerr << "error: cannot write the logs to " << dir << "\n";
            return 1;
        }
        std::printf("wrote %s and %s\n", poseFile.string().c_str(), velocityFile.string().c_str());
        return pose.ok && velocity.ok ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
