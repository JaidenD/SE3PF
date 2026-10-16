// The hardware experiments (docs/HARDWARE.md) on the MuJoCo Menagerie Panda: the best prediction of the robot.
//
// For every experiment, in the order the lab runs them:
//   preflight   kinematic only, with MuJoCo's arm model, and with the arm + Franka Hand model (MujocoDynamics);
//   closed loop at the torque level, run as the robot driver runs it (ExperimentLoop, safety stop, end conditions),
//     with the command held for each 1 ms period, libfranka's low-pass filter and rate limiter, the robot's G added,
//     and the plant's joint damping (1 N m s/rad, the robot's friction) uncompensated, as run C of panda_mujoco:
//       arm   plant panda_nohand.xml;
//       hand  plant panda.xml (the Franka Hand as a load), M and C qd of the arm alone, G with the hand (the hand
//             configured in Desk, as run D-G).
//   safety stop (SafetyStop, the damped hold): the hand missing in Desk (plant panda.xml, G without the hand: run
//     D) must be stopped before |xi| leaves its bound and held within 0.03 rad of where it stopped, at rest; and the
//     stop alone from 50 % of every velocity limit must travel less than the 0.10 rad joint margin and settle.
// Writes out/mujoco_<experiment>_<arm|hand>.csv in the robot log's format (apps/robot/robot_log.hpp), to overlay on
// a robot log: tools/plot_robot_log.py <robot log> --overlay out/mujoco_paper_hand.csv
//   Usage:  experiments_mujoco [output directory, default out/] [Menagerie franka_emika_panda directory]
//   Exit status: 0 if every preflight passes, every closed loop ends at rest with no safety stop, and every safety
//   stop above holds the arm as required; 2 otherwise.
#include "mujoco_experiments.hpp"
#include "robot_log.hpp"
#include "se3pf/experiments.hpp"
#include "se3pf/panda_example.hpp"
#include "se3pf/preflight.hpp"

#include <algorithm>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <string>
#include <vector>

using namespace se3pf;

int main(int argc, char** argv)
try {
    const std::filesystem::path out = argc > 1 ? argv[1] : "out";
    const std::string modelDir = argc > 2 ? argv[2] : defaultPandaModelDir();
    std::filesystem::create_directories(out);
    const std::string noHand = modelDir + "/panda_nohand.xml", withHand = modelDir + "/panda.xml";

    const PandaExample example;
    const MujocoDynamics armModel(noHand), handModel(withHand);
    bool allOk = true;

    std::printf("The hardware experiments on the MuJoCo Panda (bounds in brackets; docs/HARDWARE.md)\n\n");
    for (const ExperimentSpec& spec : experiments()) {
        const ExperimentSetup setup(example, spec);
        const SafetyMonitor monitor = setup.monitor();
        std::printf("=== %s: %s (%.0f s)\n", spec.name.c_str(), spec.title.c_str(), spec.duration);

        struct Pre {
            const char* name;
            const DynamicsModel* dyn;
        };
        for (const Pre& p : {Pre{"kinematic", nullptr}, Pre{"MuJoCo arm", &armModel}, Pre{"MuJoCo arm + hand", &handModel}}) {
            const PreflightReport rep = preflight(setup, p.dyn, monitor);
            std::printf("  preflight, %s: %s\n", p.name, rep.summary().c_str());
            allOk = allOk && rep.ok;
        }

        struct Loop {
            const char* name;
            MujocoExperimentOptions options;
        };
        MujocoExperimentOptions arm;
        arm.plant = noHand;
        MujocoExperimentOptions hand;
        hand.plant = withHand;
        hand.model = noHand;
        hand.plantGravity = true;
        for (const Loop& l : {Loop{"arm", arm}, Loop{"hand", hand}}) {
            const MujocoExperimentRun run = simulateExperimentMujoco(setup, monitor, l.options);
            const bool ok = run.finished && run.violation == Violation::None && run.warnings == 0;
            allOk = allOk && ok;
            std::printf("  closed loop, %s (hold 1 ms, libfranka filter, damping 1 N m s/rad): %s at t = %.3f s",
                        l.name, ok ? "finished at rest" : "FAILED", run.tEnd);
            if (run.violation != Violation::None)
                std::printf(" (safety stop: %s at t = %.3f s)", toString(run.violation), run.tStop);
            std::printf("\n    %s\n", run.peaks.summary(monitor.limits()).c_str());
            std::printf("    last 5 s of following: |xi| <= %.2e, |beta| <= %.2e rad", run.steadyXi, run.steadyBeta);
            if (spec.mode == TangentialReference::Mode::Velocity)
                std::printf("; mean eta_dot while cruising %.4f (nu_d = %.3f)", run.cruiseEtaDot, spec.speed.nuMax);
            std::printf("; contacts %d, clamped periods %ld, MuJoCo warnings %d\n", run.maxContacts,
                        run.clampedPeriods, run.warnings);
            const std::filesystem::path file = out / ("mujoco_" + spec.name + "_" + l.name + ".csv");
            if (!robot::writeLogCsv(file.string(), run.log, run.log.size())) {
                std::printf("    could not write %s\n", file.string().c_str());
                allOk = false;
            }
        }

        // The hand missing in Desk: the monitor must stop the sagging arm, and the stop must hold it
        {
            MujocoExperimentOptions missing;
            missing.plant = withHand;
            missing.model = noHand;
            missing.plantGravity = false;
            const MujocoExperimentRun run = simulateExperimentMujoco(setup, monitor, missing);
            double maxXi = 0.0;
            for (const ExperimentLogRow& row : run.log)
                if (row.phase == 0)
                    maxXi = std::max(maxXi, row.xiNorm);
            const bool ok = run.violation != Violation::None && maxXi <= monitor.limits().xi && run.stopSettled
                            && run.maxStopDrift < 0.03 && run.endSpeed < RunConstants::stopSpeed && run.warnings == 0;
            allOk = allOk && ok;
            std::printf("  safety stop, hand missing in Desk: %s at t = %.3f s (|xi| <= %.4f before it); held for "
                        "%.3f s, largest joint drift %.4f rad [< 0.03], final |dq| %.4f rad/s: %s\n",
                        toString(run.violation), run.tStop, maxXi, run.tEnd - run.tStop, run.maxStopDrift,
                        run.endSpeed, ok ? "OK" : "FAILED");
        }
        // The stop alone from 50 % of every velocity limit (alternating signs), at the start, with the hand
        {
            Vec7 qd0 = 0.5 * PandaArm::dqMax();
            for (int i = 1; i < 7; i += 2)
                qd0[i] = -qd0[i];
            const MujocoStopRun r = simulateStopMujoco(withHand, "", setup.start(), qd0);
            const bool ok = r.drift.maxCoeff() < HardwareBounds::jointMargin && r.settled && r.warnings == 0;
            allOk = allOk && ok;
            std::printf("  safety stop from 50 %% of every velocity limit: travel %.3f rad [< %.2f], at rest after "
                        "%.3f s, peak stop torque %.1f N m: %s\n",
                        r.drift.maxCoeff(), HardwareBounds::jointMargin, r.duration, r.maxTorque, ok ? "OK" : "FAILED");
        }
        std::printf("\n");
    }
    std::printf("%s\n", allOk ? "ALL PASSED" : "SOME FAILED");
    return allOk ? 0 : 2;
} catch (const std::exception& e) {
    std::fprintf(stderr, "experiments_mujoco: %s\n", e.what());
    return 2;
}
