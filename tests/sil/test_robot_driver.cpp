// Software in the loop: the robot driver (apps/robot/robot_driver.cpp) built against mock/franka, run end to end as
// a separate process (robot_driver_mock), as the operator runs robot_driver on the lab PC: every experiment with
// --yes (the mock-only option that answers the prompts), the operator's answers through standard input, faults
// injected into the simulated robot through the mock's environment variables (mock/franka/README.md), every kind of
// violation through the driver's safety stop (--inject, mock only), the friction feedforward against a plant with
// Coulomb friction (CTRL-1), and the driver's refusals. Checks the exit codes, the output and the CSV logs. A
// separate test program (tests_sil); the runs take about 45 s in all in Release.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "se3pf/experiments.hpp"
#include "se3pf/panda.hpp"
#include "test_helpers.hpp"

#include <sys/wait.h>
#include <unistd.h>  // getpid

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <utility>
#include <sstream>
#include <string>
#include <vector>

#ifndef SE3PF_ROBOT_DRIVER_MOCK
#error "SE3PF_ROBOT_DRIVER_MOCK must be the path of robot_driver_mock (CMakeLists.txt)"
#endif

using namespace se3pf;
namespace fs = std::filesystem;

namespace {

    // The runs' logs (tens of MB) go to a temporary directory, removed at exit unless SE3PF_SIL_KEEP is set. Every run
    // that can reach the preflight gets a --log there (the default log path is in the working directory).
    struct TempDir {
        fs::path path;
        TempDir() : path(fs::temp_directory_path() / ("se3pf_sil_" + std::to_string(::getpid())))
        {
            fs::create_directories(path);
        }
        ~TempDir()
        {
            std::error_code ec;
            if (std::getenv("SE3PF_SIL_KEEP") == nullptr)
                fs::remove_all(path, ec);
        }
    };

    fs::path workDir()
    {
        static const TempDir dir;
        return dir.path;
    }

    struct Run {
        int code = -1;
        std::string out;
        bool contains(const std::string& s) const { return out.find(s) != std::string::npos; }
    };

    // Runs robot_driver_mock with the arguments, environment assignments (VAR=value ...) and standard input.
    Run run(const std::string& args, const std::string& env = "", const std::string& input = "")
    {
        const fs::path outFile = workDir() / "out.txt";
        const fs::path inFile = workDir() / "in.txt";
        {
            std::ofstream in(inFile);
            in << input;
        }
        const std::string cmd = "env " + env + " '" + std::string(SE3PF_ROBOT_DRIVER_MOCK) + "' " + args + " < '"
                                + inFile.string() + "' > '" + outFile.string() + "' 2>&1";
        const int status = std::system(cmd.c_str());
        Run r;
        r.code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
        std::ifstream f(outFile);
        std::stringstream ss;
        ss << f.rdbuf();
        r.out = ss.str();
        return r;
    }

    // A CSV log of the driver (apps/robot/robot_log.hpp), by column name
    struct Csv {
        std::map<std::string, std::size_t> column;
        std::vector<std::vector<double>> rows;
        double at(std::size_t row, const std::string& name) const { return rows.at(row).at(column.at(name)); }
        Vec7 vec(std::size_t row, const std::string& prefix) const
        {
            Vec7 v = Vec7::Zero();
            for (int i = 0; i < 7; ++i)
                v[i] = at(row, prefix + std::to_string(i + 1));
            return v;
        }
    };

    Csv readCsv(const fs::path& path)
    {
        Csv c;
        std::ifstream f(path);
        std::string line;
        if (!std::getline(f, line))
            return c;
        std::stringstream header(line);
        std::string name;
        for (std::size_t k = 0; std::getline(header, name, ','); ++k)
            c.column[name] = k;
        while (std::getline(f, line)) {
            std::vector<double> row;
            std::stringstream ss(line);
            std::string cell;
            while (std::getline(ss, cell, ','))
                row.push_back(std::strtod(cell.c_str(), nullptr));
            c.rows.push_back(row);
        }
        return c;
    }

    std::string logArg(const std::string& name) { return "--log=" + (workDir() / (name + ".csv")).string(); }
    fs::path logPath(const std::string& name) { return workDir() / (name + ".csv"); }

    // FRANKA_MOCK_COULOMB=c1,...,c7: the plant's Coulomb friction, e.g. what an experiment's feedforward compensates
    std::string coulombEnv(const Vec7& c)
    {
        std::string env = "FRANKA_MOCK_COULOMB=";
        for (int i = 0; i < 7; ++i)
            env += (i ? "," : "") + std::to_string(c[i]);
        return env;
    }

    // Mean eta_dot over the cruise of a velocity reference (rows with ref at its maximum), and the mean |dq_i|
    // there; NaN without such rows
    struct Cruise {
        double etaDot = std::nan("");
        Vec7 meanAbsDq = Vec7::Constant(std::nan(""));
        double maxXi = 0.0, maxBeta = 0.0;  // over every following row
    };
    Cruise cruiseOf(const Csv& log)
    {
        Cruise c;
        double nuMax = 0.0;
        for (std::size_t k = 0; k < log.rows.size(); ++k)
            nuMax = std::max(nuMax, log.at(k, "ref"));
        double sum = 0.0;
        Vec7 dqSum = Vec7::Zero();
        long n = 0;
        for (std::size_t k = 0; k < log.rows.size(); ++k) {
            if (log.at(k, "phase") != 0.0)
                continue;
            c.maxXi = std::max(c.maxXi, log.at(k, "xi_norm"));
            c.maxBeta = std::max(c.maxBeta, std::abs(log.at(k, "beta")));
            if (nuMax > 0.0 && log.at(k, "ref") >= 0.999 * nuMax) {
                sum += log.at(k, "eta_dot");
                dqSum += log.vec(k, "dq").cwiseAbs();
                ++n;
            }
        }
        if (n > 0) {
            c.etaDot = sum / n;
            c.meanAbsDq = dqSum / n;
        }
        return c;
    }

    // The first row of the safety stop (phase 1), or rows.size()
    std::size_t firstStopRow(const Csv& log)
    {
        for (std::size_t k = 0; k < log.rows.size(); ++k)
            if (log.at(k, "phase") == 1.0)
                return k;
        return log.rows.size();
    }

    // The safety stop in a log: from the first stop row on, every row is a stop row whose torque is the damped hold
    // at the configuration of the first one (to the CSV's 10 digits: K_s (q - q_stop) to about 3e-7 N m); the largest
    // joint drift from it; the last row at rest.
    struct StopCheck {
        std::size_t first = 0;
        double tStop = -1.0, duration = 0.0, maxDrift = 0.0, maxLawError = 0.0, endSpeed = 0.0;
        bool allStopRows = true;
    };
    StopCheck checkStop(const Csv& log)
    {
        StopCheck c;
        c.first = firstStopRow(log);
        if (c.first >= log.rows.size())
            return c;
        const Vec7 qStop = log.vec(c.first, "q");
        c.tStop = log.at(c.first, "t");
        for (std::size_t k = c.first; k < log.rows.size(); ++k) {
            c.allStopRows = c.allStopRows && log.at(k, "phase") == 1.0;
            const Vec7 q = log.vec(k, "q"), dq = log.vec(k, "dq");
            c.maxLawError = std::max(c.maxLawError, testing::maxDiff(log.vec(k, "tau"), SafetyStop::torque(q, dq, qStop)));
            c.maxDrift = std::max(c.maxDrift, (q - qStop).cwiseAbs().maxCoeff());
        }
        c.duration = log.at(log.rows.size() - 1, "t") - c.tStop;
        c.endSpeed = log.vec(log.rows.size() - 1, "dq").norm();
        return c;
    }

}  // namespace

TEST_CASE("robot_driver_mock: list, usage errors, unknown experiments")
{
    const Run list = run("list");
    CHECK(list.code == 0);
    for (const ExperimentSpec& e : experiments())
        CHECK(list.contains(e.name));
    CHECK(run("").code == 1);
    CHECK(run("0.0.0.0").code == 1);
    CHECK(run("0.0.0.0 hold --bogus").code == 1);
    CHECK(run("0.0.0.0 hold --log=").code == 1);
    CHECK(run("0.0.0.0 hold --yes --inject=bogus@1").code == 1);
    CHECK(run("0.0.0.0 hold --yes --inject=tube").code == 1);
    CHECK(run("0.0.0.0 hold --yes --inject=tube@x").code == 1);
    CHECK(run("0.0.0.0 hold --yes --friction-scale=").code == 1);
    CHECK(run("0.0.0.0 hold --yes --friction-scale=x").code == 1);
    CHECK(run("0.0.0.0 hold --yes --friction-scale=-1").code == 1);
    CHECK(run("0.0.0.0 hold --yes --friction-scale=0.5x").code == 1);
    CHECK(run("0.0.0.0 no-such-experiment --yes").code == 2);
    CHECK(run("0.0.0.0 preflight no-such-experiment").code == 2);
}

TEST_CASE("robot_driver_mock: check and preflight on the simulated robot; no motion")
{
    const Run check = run("0.0.0.0 check");
    INFO(check.out);
    CHECK(check.code == 0);
    CHECK(check.contains("Kinematics check"));
    CHECK(check.contains("OK"));
    CHECK(check.contains("realtime kernel"));
    CHECK(check.contains("without swapping libfranka's (v, omega) rows"));
    CHECK(check.contains("m_total = 0.730 kg"));  // the mock's configured end effector: a Franka Hand

    const fs::path dir = workDir() / "predictions";
    fs::create_directories(dir);
    const Run pre = run("0.0.0.0 preflight all --log-dir=" + dir.string());
    INFO(pre.out);
    CHECK(pre.code == 0);
    for (const ExperimentSpec& e : experiments()) {
        CAPTURE(e.name);
        CHECK(pre.contains("Preflight [" + e.name + "] with the robot's model: OK"));
        const Csv c = readCsv(dir / ("preflight_" + e.name + ".csv"));
        CHECK(c.rows.size() > static_cast<std::size_t>(e.duration * 1000.0));
    }
    const Run one = run("0.0.0.0 preflight paper");
    CHECK(one.code == 0);
    CHECK(one.contains("Preflight [paper]"));
    CHECK_FALSE(one.contains("Preflight [hold]"));
}

TEST_CASE("robot_driver_mock: every experiment end to end: approach, torque loop, at rest, log", SLOW_IN_DEBUG)
{
    for (const ExperimentSpec& e : experiments()) {
        CAPTURE(e.name);
        // The simulated robot has the Coulomb friction the experiment's feedforward compensates (loops: Gaz et
        // al. 2019 with joint 6 raised; the others: none)
        const bool friction = e.friction.active();
        const Run r = run("0.0.0.0 " + e.name + " --yes " + logArg(e.name), coulombEnv(e.friction.coulomb));
        INFO(r.out);
        REQUIRE(r.code == 0);
        CHECK(r.contains("Approach done"));
        CHECK(r.contains("At the start: |xi| = "));
        CHECK(r.contains("Preflight [" + e.name + "] from the measured start: OK"));
        CHECK(r.contains("Finished at rest."));
        CHECK_FALSE(r.contains("SAFETY STOP"));
        CHECK(r.contains("callback allocations: 0"));  // the 1 kHz callback allocates nothing

        const Csv log = readCsv(logPath(e.name));
        const Csv pred = readCsv(fs::path(logPath(e.name).string() + ".preflight.csv"));
        REQUIRE(log.rows.size() > static_cast<std::size_t>(e.duration * 1000.0));
        REQUIRE(pred.rows.size() > 0);
        const std::size_t last = log.rows.size() - 1;
        CHECK(log.at(last, "t") >= e.duration);
        CHECK(log.vec(last, "dq").norm() < RunConstants::restSpeed);
        double maxXiDiff = 0.0, maxXi = 0.0, maxBeta = 0.0;
        for (std::size_t k = 0; k < log.rows.size(); ++k) {
            CHECK(log.at(k, "phase") == 0.0);
            maxXi = std::max(maxXi, log.at(k, "xi_norm"));
            maxBeta = std::max(maxBeta, std::abs(log.at(k, "beta")));
            if (k < pred.rows.size())
                maxXiDiff = std::max(maxXiDiff, std::abs(log.at(k, "xi_norm") - pred.at(k, "xi_norm")));
        }
        // The simulated robot follows the prediction (from the measured start): the differences are the low-pass
        // filter of the torque command and the torque held over each period. With Coulomb friction in the plant
        // (loops) the feedforward cancels it while a joint moves, not while it sticks (the start, each reversal of
        // a joint of the lift), which the frictionless prediction does not see: |xi| 2e-3, |beta| 6e-3, well inside
        // the experiment's pass criteria (2e-2, 3e-2), and a small error stays at rest (static friction holds it,
        // as the robot's did in 'hold').
        CAPTURE(maxXiDiff, maxXi, maxBeta);
        CHECK(maxXiDiff < (friction ? 5e-3 : 1e-4));  // paper: 2.9e-5, the O(dt) of MuJoCo run B (3.4e-5)
        CHECK(maxXi <= HardwareBounds::xi);
        CHECK(maxBeta <= (friction ? 1e-2 : std::max(1e-3, std::abs(e.beta0) + 1e-4)));
        CHECK(log.at(last, "xi_norm") < 1e-3);
        CHECK(std::abs(log.at(last, "beta")) < (friction ? 5e-3 : 1e-3));
        // Two loops: the arm ends where it started (with friction, a little short along the path: eta_dot lags
        // nu_d by about 1 %, so the end is a few hundredths of a radian before the start, on the lift)
        if (e.name == "loops") {
            if (friction)
                CHECK(std::abs(testing::arcDiff(log.at(last, "eta"), 0.1, 0.77008765)) < 0.05);
            else
                CHECK(testing::maxDiff(log.vec(last, "q"), log.vec(0, "q")) < 2e-3);
        }
        // The friction feedforward in the log: zero at the first cycle, never above Fc + Fv |v|, part of tau
        if (friction) {
            const Vec7 fc = e.friction.coulomb;
            CHECK(log.vec(0, "tauF").isZero());
            double maxTauF = 0.0;
            for (std::size_t k = 0; k < log.rows.size(); ++k) {
                const Vec7 tauF = log.vec(k, "tauF");
                CHECK((tauF.cwiseAbs().array() <= fc.array() + 1e-9).all());
                maxTauF = std::max(maxTauF, tauF.cwiseAbs().maxCoeff());
            }
            CHECK(maxTauF > 0.9 * fc.maxCoeff());
            CHECK(r.contains("friction feedforward tau_f"));
        } else {
            CHECK(r.contains("friction feedforward: off"));
        }
    }
}

TEST_CASE("robot_driver_mock: every kind of violation goes through the safety stop, which holds the arm at rest")
{
    // --inject (mock build only) makes the driver treat that violation as found 1 s into the torque loop of 'hold';
    // the checks that find each one are tested in tests/test_safety.cpp and tests/test_preflight.cpp.
    const Violation all[] = {Violation::Projection, Violation::JointLimit, Violation::JointVelocity,
                             Violation::Acceleration, Violation::Torque, Violation::TorqueRate,
                             Violation::FlangeSpeed, Violation::Tube, Violation::Redundancy,
                             Violation::Conditioning, Violation::NonFinite, Violation::Timeout};
    const char* names[] = {"projection", "joint-limit", "joint-velocity", "acceleration", "torque", "torque-rate",
                           "flange-speed", "tube", "redundancy", "conditioning", "non-finite", "timeout"};
    for (int k = 0; k < 12; ++k) {
        CAPTURE(names[k]);
        const std::string name = std::string("inject_") + names[k];
        const Run r = run("0.0.0.0 hold --yes --inject=" + std::string(names[k]) + "@1.0 " + logArg(name));
        INFO(r.out);
        CHECK(r.code == 5);
        CHECK(r.contains(std::string("SAFETY STOP: ") + toString(all[k]) + " at t = 1.0"));
        CHECK(r.contains("(at rest)"));
        CHECK(r.contains("callback allocations: 0"));
        const Csv log = readCsv(logPath(name));
        const StopCheck c = checkStop(log);
        REQUIRE(c.first < log.rows.size());
        CHECK(c.tStop >= 1.0);
        CHECK(c.tStop < 1.0 + 2e-3);
        CHECK(c.allStopRows);
        CHECK(c.maxLawError < 1e-6);
        CHECK(c.maxDrift < 1e-3);  // at rest when it stopped: it stays
        CHECK(c.endSpeed < RunConstants::stopSpeed);
        CHECK(c.duration >= RunConstants::stopSettle);
        CHECK(c.duration < 0.2);
    }
    // While moving at the paper's top speed: the hold stops the arm within a few centimetres of joint travel
    const Run r = run("0.0.0.0 paper --yes --inject=tube@6 " + logArg("inject_moving"));
    INFO(r.out);
    CHECK(r.code == 5);
    const Csv log = readCsv(logPath("inject_moving"));
    const StopCheck c = checkStop(log);
    REQUIRE(c.first < log.rows.size());
    CAPTURE(c.maxDrift, c.duration, c.endSpeed);
    CHECK(log.vec(c.first, "dq").norm() > 0.1);  // it was moving
    CHECK(c.maxLawError < 1e-6);
    CHECK(c.maxDrift < 0.05);
    CHECK(c.endSpeed < RunConstants::stopSpeed);
    CHECK(c.duration < 1.0);
}

TEST_CASE("robot_driver_mock: a push on joint 4 during hold: the safety stop engages and holds the arm")
{
    // 12 N m on joint 4 for 0.3 s from t = 3 s (below the robot's collision threshold of 18 N m), as the reference
    // implementation's software-in-the-loop test
    const Run r = run("0.0.0.0 hold --yes " + logArg("push"), "FRANKA_MOCK_PUSH=4,12,3,0.3");
    INFO(r.out);
    CHECK(r.code == 5);
    CHECK(r.contains("SAFETY STOP"));
    CHECK(r.contains("(at rest)"));
    CHECK_FALSE(r.contains("REFLEX"));
    CHECK(r.contains("callback allocations: 0"));
    const Csv log = readCsv(logPath("push"));
    REQUIRE(log.rows.size() > 3000);
    const StopCheck c = checkStop(log);
    REQUIRE(c.first < log.rows.size());
    CAPTURE(c.tStop, c.maxDrift, c.duration, c.endSpeed);
    CHECK(c.tStop > 3.0);
    CHECK(c.tStop < 3.1);  // within 0.1 s of the push
    // From then on: the damped hold at the configuration of the violation
    CHECK(c.allStopRows);
    CHECK(c.maxLawError < 1e-6);
    // The push deflects the hold while it lasts (12 N m against 300 N m/rad), then the arm comes back and rests
    CHECK(c.maxDrift < 0.1);
    CHECK(testing::maxDiff(log.vec(log.rows.size() - 1, "q"), log.vec(c.first, "q")) < 5e-3);
    CHECK(c.endSpeed < RunConstants::stopSpeed);
    CHECK(c.duration < RunConstants::stopTime);
    // The arm moved, but stayed far from its joint limits
    double margin = 1e9;
    for (std::size_t k = 0; k < log.rows.size(); ++k)
        margin = std::min(margin, PandaArm::jointMargin(log.vec(k, "q")));
    CHECK(margin > 0.3);
}

TEST_CASE("robot_driver_mock: an end effector missing in Desk: the monitor stops the sagging arm, the hold keeps it")
{
    // A mass at the flange in the plant that the configured load (and so G) does not have: MuJoCo run D's failure.
    // The first violation comes within 0.2 s; the damped hold then keeps every joint within 0.03 rad of where it
    // stopped for the hand's 0.73 kg (0.028; damping alone let joint 2 drift 0.5 rad here), within 0.04 rad for 1 kg
    // (0.035: the toy arm's large joint inertias; MuJoCo's Panda with its 0.76 kg hand: 0.022), and the arm comes to
    // rest.
    for (const auto& [load, driftBound] : {std::pair<const char*, double>{"0.73", 0.03}, {"1.0", 0.04}}) {
        CAPTURE(load);
        const std::string name = std::string("sag_") + load;
        const Run r = run("0.0.0.0 paper --yes " + logArg(name), std::string("FRANKA_MOCK_PLANT_LOAD=") + load);
        INFO(r.out);
        CHECK(r.code == 5);
        CHECK(r.contains("SAFETY STOP"));
        CHECK(r.contains("(at rest)"));
        const Csv log = readCsv(logPath(name));
        REQUIRE(!log.rows.empty());
        const StopCheck c = checkStop(log);
        REQUIRE(c.first < log.rows.size());
        CAPTURE(c.tStop, c.maxDrift, c.duration, c.endSpeed);
        CHECK(c.tStop < 0.2);
        double maxXi = 0.0;
        for (std::size_t k = 0; k < c.first; ++k)
            maxXi = std::max(maxXi, log.at(k, "xi_norm"));
        CHECK(maxXi <= HardwareBounds::xi);
        CHECK(c.maxLawError < 1e-6);
        CHECK(c.maxDrift < driftBound);
        CHECK(c.endSpeed < RunConstants::stopSpeed);
        CHECK(c.duration < RunConstants::stopTime);
    }
}

TEST_CASE("robot_driver_mock: a NaN joint velocity from the robot: NonFinite, and the hold still holds")
{
    const Run r = run("0.0.0.0 hold --yes " + logArg("nan"), "FRANKA_MOCK_NAN=4,2.0");
    INFO(r.out);
    CHECK(r.code == 5);
    CHECK(r.contains(std::string("SAFETY STOP: ") + toString(Violation::NonFinite)));
    CHECK(r.contains("(at rest)"));
    CHECK_FALSE(r.contains("ERROR"));  // the command stayed finite
    const Csv log = readCsv(logPath("nan"));
    const StopCheck c = checkStop(log);
    REQUIRE(c.first < log.rows.size());
    CHECK(c.tStop >= 2.0);
    CHECK(c.tStop < 2.0 + 2e-3);
    CHECK(c.maxDrift < 1e-3);
    CHECK(c.endSpeed < RunConstants::stopSpeed);
}

TEST_CASE("robot_driver_mock: joint friction (1 N m s/rad, as MuJoCo run C) leaves errors but keeps the margins",
          SLOW_IN_DEBUG)
{
    const Run r = run("0.0.0.0 paper --yes " + logArg("friction"), "FRANKA_MOCK_FRICTION=1");
    INFO(r.out);
    CHECK(r.code == 0);
    CHECK_FALSE(r.contains("SAFETY STOP"));
    const Csv log = readCsv(logPath("friction"));
    REQUIRE(!log.rows.empty());
    double maxBeta = 0.0, minEtaDot = 1.0;
    for (std::size_t k = 0; k < log.rows.size(); ++k) {
        maxBeta = std::max(maxBeta, std::abs(log.at(k, "beta")));
        if (log.at(k, "t") > 6.0 && log.at(k, "t") < 15.0)
            minEtaDot = std::min(minEtaDot, log.at(k, "eta_dot"));
    }
    CAPTURE(maxBeta, minEtaDot);
    CHECK(maxBeta > 1e-2);   // the friction is felt ...
    CHECK(maxBeta < HardwareBounds::beta);  // ... within the margin
    CHECK(minEtaDot < 0.079);  // the path speed lags nu_d = 0.08 (no integral action: CTRL-1)
}

TEST_CASE("robot_driver_mock: Coulomb friction (the stall of 2026-10-19, CTRL-1): loops stalls without the friction "
          "feedforward and completes with it",
          SLOW_IN_DEBUG)
{
    const ExperimentSpec& loops = *findExperiment("loops");
    REQUIRE(loops.friction.active());
    const Vec7 fc = loops.friction.coulomb;  // Gaz et al. 2019, joint 6 raised
    const double nuD = loops.speed.nuMax;

    // 1. The stall of loops_1: the plant with the Coulomb levels of Gaz et al. 2019 and no feedforward
    // (--friction-scale=0, mock only: the spec as it ran on the robot). As on the
    // robot (logs/loops_1.csv: eta_dot 1 % of nu_d while cruising, |xi| 0.018, beta 0.030, the wrist never moving
    // under a commanded torque below its friction), the arm stays put: no reflex, no stop, "finished at rest" with
    // the errors frozen.
    {
        const Vec7 plant = FrictionFeedforward::gaz2019().coulomb;
        const Run r = run("0.0.0.0 loops --yes --friction-scale=0 " + logArg("stall"), coulombEnv(plant));
        INFO(r.out);
        CHECK(r.code == 0);
        CHECK(r.contains("friction feedforward: off"));
        CHECK_FALSE(r.contains("SAFETY STOP"));
        CHECK_FALSE(r.contains("REFLEX"));
        const Csv log = readCsv(logPath("stall"));
        REQUIRE(log.rows.size() > 30000);
        const Cruise c = cruiseOf(log);
        CAPTURE(c.etaDot, c.maxXi, c.maxBeta, c.meanAbsDq.transpose());
        CHECK(c.etaDot < 0.15 * nuD);      // 7 % here, 1 % on the robot
        CHECK(c.maxBeta > 0.03);           // 0.049 here, 0.030 on the robot
        CHECK(c.maxBeta <= HardwareBounds::beta);
        CHECK(c.maxXi > 0.01);             // 0.019 here, 0.018 on the robot
        for (int i : {0, 4, 6})            // joints 1, 5 and 7 never move
            CHECK(c.meanAbsDq[i] < 1e-3);
        // The commanded torque on the stuck wrist joints stays below their Coulomb level while qdd* asks for more
        // than 10 % of their acceleration limit: the controller cannot break them free
        double maxTau5 = 0.0, maxTau7 = 0.0, maxDdq = 0.0;
        for (std::size_t k = 0; k < log.rows.size(); ++k) {
            const Vec7 tau = log.vec(k, "tau").cwiseAbs();
            maxTau5 = std::max(maxTau5, tau[4]);
            maxTau7 = std::max(maxTau7, tau[6]);
            maxDdq = std::max(maxDdq, log.vec(k, "qdd").cwiseAbs().cwiseQuotient(PandaArm::ddqMax()).maxCoeff());
        }
        CAPTURE(maxTau5, maxTau7, maxDdq);
        CHECK(maxTau5 < plant[4]);
        CHECK(maxTau7 < plant[6]);
        CHECK(maxDdq > 0.10);
    }
    // 2. The plant with exactly the friction the spec compensates, without and with the feedforward
    {
        const Run off = run("0.0.0.0 loops --yes --friction-scale=0 " + logArg("coulomb_off"), coulombEnv(fc));
        INFO(off.out);
        CHECK(off.code == 0);
        CHECK_FALSE(off.contains("REFLEX"));
        const Cruise c = cruiseOf(readCsv(logPath("coulomb_off")));
        CAPTURE(c.etaDot, c.maxXi, c.maxBeta);
        CHECK(c.etaDot < 0.5 * nuD);   // 28 %: the arm lags far behind the reference ...
        CHECK(c.maxBeta > 0.02);       // ... and |beta| 0.040 reaches what the robot showed (0.030)
        CHECK(c.maxXi > 5e-3);

        const Run on = run("0.0.0.0 loops --yes " + logArg("coulomb_on"), coulombEnv(fc));
        INFO(on.out);
        CHECK(on.code == 0);
        CHECK(on.contains("Finished at rest."));
        CHECK_FALSE(on.contains("SAFETY STOP"));
        CHECK(on.contains("callback allocations: 0"));
        CHECK(on.contains("tau_f <= 0.64 N m (joint 4)"));
        const Csv log = readCsv(logPath("coulomb_on"));
        const Cruise d = cruiseOf(log);
        CAPTURE(d.etaDot, d.maxXi, d.maxBeta);
        CHECK(d.etaDot > 0.95 * nuD);  // 99 %
        CHECK(d.etaDot < 1.05 * nuD);
        CHECK(d.maxXi < 5e-3);         // 2e-3 (pass criterion 2e-2)
        CHECK(d.maxBeta < 1e-2);       // 6e-3 (3e-2)
        for (int i = 0; i < 7; ++i)    // every joint moves along the lift
            CHECK(d.meanAbsDq[i] > 0.01);
        // tau_f in the log: within Fc, reaching it on every joint, and smooth (the commanded torque rate stays far
        // below the 600 N m/s margin)
        Vec7 maxTauF = Vec7::Zero();
        double maxRate = 0.0;
        for (std::size_t k = 0; k < log.rows.size(); ++k) {
            const Vec7 tauF = log.vec(k, "tauF");
            CHECK((tauF.cwiseAbs().array() <= fc.array() + 1e-9).all());
            maxTauF = maxTauF.cwiseMax(tauF.cwiseAbs());
            if (k > 0)
                maxRate = std::max(maxRate, (log.vec(k, "tau") - log.vec(k - 1, "tau")).cwiseAbs().maxCoeff() / 1e-3);
        }
        CAPTURE(maxTauF.transpose(), maxRate);
        CHECK((maxTauF.array() > 0.9 * fc.array()).all());
        CHECK(maxRate < 20.0);  // 3 N m/s
        // The prediction carries the same feedforward
        const Csv pred = readCsv(fs::path(logPath("coulomb_on").string() + ".preflight.csv"));
        REQUIRE(pred.rows.size() > 30000);
        double predMax = 0.0;
        for (std::size_t k = 0; k < pred.rows.size(); ++k)
            predMax = std::max(predMax, pred.vec(k, "tauF").cwiseAbs().maxCoeff());
        CHECK_THAT(predMax, Catch::Matchers::WithinAbs(fc.maxCoeff(), 1e-3));
    }
}

TEST_CASE("robot_driver_mock: late periods (lost packets): the run keeps time and its margins", SLOW_IN_DEBUG)
{
    // Every 50th period lasts 3 ms (the robot keeps the last command; the next callback gets a 3 ms period)
    const Run r = run("0.0.0.0 paper --yes " + logArg("late"), "FRANKA_MOCK_LATE=50,3");
    INFO(r.out);
    CHECK(r.code == 0);
    CHECK(r.contains("Finished at rest."));
    CHECK_FALSE(r.contains("0 periods longer than 1 ms"));
    CHECK(r.contains("lowest command success rate 0.96"));
    const Csv log = readCsv(logPath("late"));
    REQUIRE(!log.rows.empty());
    long late = 0;
    double maxXi = 0.0;
    for (std::size_t k = 1; k < log.rows.size(); ++k) {
        const double dt = log.at(k, "t") - log.at(k - 1, "t");
        CHECK_THAT(dt, Catch::Matchers::WithinAbs(1e-3 * log.at(k, "period_ms"), 1e-9));  // t follows libfranka's periods
        late += log.at(k, "period_ms") > 1.0 ? 1 : 0;
        maxXi = std::max(maxXi, log.at(k, "xi_norm"));
    }
    CHECK(late > 300);
    CHECK(maxXi <= HardwareBounds::xi);
    CHECK(log.at(log.rows.size() - 1, "t") >= 20.0);
}

TEST_CASE("robot_driver_mock: the driver refuses to move: kinematics, robot mode, distance, operator, dry run")
{
    SECTION("a model whose flange is 0.1 mm off fails the kinematics check (exit 3), nothing moves")
    {
        const Run r = run("0.0.0.0 hold --yes " + logArg("kin"), "FRANKA_MOCK_POSE_ERROR=1e-4");
        CHECK(r.code == 3);
        CHECK(r.contains("MISMATCH"));
        CHECK(r.contains("position 1.0e-04 m"));  // the diagnosis names the kind of error
        CHECK_FALSE(r.contains("Approach"));
        CHECK_FALSE(fs::exists(logPath("kin")));
        CHECK(run("0.0.0.0 check", "FRANKA_MOCK_POSE_ERROR=1e-4").code == 3);
        CHECK(run("0.0.0.0 preflight", "FRANKA_MOCK_POSE_ERROR=1e-4").code == 3);
    }
    SECTION("a robot in reflex mode or with the user stop pressed (exit 7); a dry run works in any mode")
    {
        for (const char* mode : {"reflex", "user-stopped"}) {
            const Run r = run("0.0.0.0 hold --yes " + logArg("mode"), std::string("FRANKA_MOCK_MODE=") + mode);
            INFO(r.out);
            CHECK(r.code == 7);
            CHECK(r.contains("Not ready"));
            CHECK_FALSE(r.contains("The robot will MOVE"));
            const Run dry = run("0.0.0.0 hold --yes --dry-run " + logArg("mode_dry"),
                                std::string("FRANKA_MOCK_MODE=") + mode);
            CHECK(dry.code == 0);
            CHECK(dry.contains("(needs Idle)"));
        }
        // error-recovery clears a reflex (the experiment named 'recover' is something else)
        const Run rec = run("0.0.0.0 error-recovery --yes", "FRANKA_MOCK_MODE=reflex");
        CHECK(rec.code == 0);
        CHECK(rec.contains("Recovered. Robot mode: Idle"));
    }
    SECTION("an arm far from the start (exit 7)")
    {
        const Run r = run("0.0.0.0 hold --yes " + logArg("far"), "FRANKA_MOCK_Q0=1.5,-0.785,0,-2.356,0,1.571,0.785");
        CHECK(r.code == 7);
        CHECK(r.contains("Guide the arm"));
    }
    SECTION("the robot's model with a load that fails the preflight (exit 4): nothing moves")
    {
        // 5 kg configured in Desk with its centre of mass 0.15 m from the flange: joint 6 needs 71 % of its limit
        const Run r = run("0.0.0.0 paper --yes " + logArg("heavy"), "FRANKA_MOCK_CONFIGURED_LOAD=5,0,0,0.15");
        INFO(r.out);
        CHECK(r.code == 4);
        CHECK(r.contains("Preflight [paper] with the robot's model: FAILED: joint torque margin"));
        CHECK(r.contains("Nothing moved"));
        CHECK_FALSE(r.contains("The robot will MOVE"));
        CHECK_FALSE(fs::exists(logPath("heavy")));
        CHECK(run("0.0.0.0 preflight paper", "FRANKA_MOCK_CONFIGURED_LOAD=5,0,0,0.15").code == 4);
    }
    SECTION("the approach ends off the start (exit 7): no torque control")
    {
        // 0.02 rad on every joint (the joint check), and 3e-3 rad (inside the joint tolerance; xi and beta are off)
        for (const char* e : {"0.02", "0.003"}) {
            CAPTURE(e);
            const Run r = run("0.0.0.0 self-motion --yes " + logArg("off"), std::string("FRANKA_MOCK_APPROACH_ERROR=") + e);
            INFO(r.out);
            CHECK(r.code == 7);
            CHECK(r.contains("Approach done"));
            CHECK(r.contains("Not ready: not at rest at the start"));
            CHECK_FALSE(r.contains("TORQUE control) will start"));
            CHECK_FALSE(fs::exists(logPath("off")));
        }
        // 3e-4 rad: accepted, and the preflight is run again from there
        const Run ok = run("0.0.0.0 self-motion --yes " + logArg("off_small"), "FRANKA_MOCK_APPROACH_ERROR=3e-4");
        CHECK(ok.code == 0);
        CHECK(ok.contains("Preflight [self-motion] from the measured start: OK"));
    }
    SECTION("the approach is planned from q_d, not from the measured q (1 mm of deflection)")
    {
        // Planned from q, the first command would differ from q_d: a start-pose reflex (exit 8)
        const Run r = run("0.0.0.0 hold --yes " + logArg("qd"), "FRANKA_MOCK_HOLD_ERROR=1e-3");
        INFO(r.out);
        CHECK(r.code == 0);
        CHECK_FALSE(r.contains("REFLEX"));
    }
    SECTION("the state is read again after each Enter: a user stop pressed while waiting at a prompt")
    {
        // The run reads the state at the kinematics check (1), before the approach prompt (2), after its Enter (3),
        // after the approach (4) and after the torque prompt's Enter (5). Pressed after read 4: the arm stays at the
        // start, no torque loop (without read 5 libfranka would refuse the motion instead: exit 8). Pressed after
        // read 2: libfranka rejects the behavior setters that follow it (exit 8), and nothing moves.
        const Run a = run("0.0.0.0 hold --yes " + logArg("stop2"), "FRANKA_MOCK_USER_STOP_AFTER_READS=2");
        INFO(a.out);
        CHECK(a.code == 8);
        CHECK(a.contains("command rejected: command not possible in the current mode (\"User stopped\")"));
        CHECK(a.contains("The robot rejected a command"));
        CHECK_FALSE(a.contains("Approach done"));
        const Run b = run("0.0.0.0 hold --yes " + logArg("stop4"), "FRANKA_MOCK_USER_STOP_AFTER_READS=4");
        INFO(b.out);
        CHECK(b.code == 7);
        CHECK(b.contains("Approach done"));
        CHECK(b.contains("robot mode UserStopped"));
        CHECK_FALSE(b.contains("Run summary"));
        CHECK_FALSE(fs::exists(logPath("stop4")));
    }
    SECTION("no realtime kernel: libfranka refuses the first motion (exit 8), nothing moves")
    {
        const Run r = run("0.0.0.0 hold --yes " + logArg("nort"), "FRANKA_MOCK_NO_RT=1");
        CHECK(r.code == 8);
        CHECK(r.contains("realtime"));
        CHECK_FALSE(r.contains("Approach done"));
    }
    SECTION("the operator: only Enter alone moves the robot (exit 6 otherwise)")
    {
        const Run noLoad = run("0.0.0.0 hold " + logArg("op1"), "", "n\n");
        CHECK(noLoad.code == 6);
        CHECK_FALSE(noLoad.contains("Approach done"));
        const Run eof = run("0.0.0.0 hold " + logArg("op1"), "", "");  // no terminal: never moves
        CHECK(eof.code == 6);
        CHECK_FALSE(eof.contains("Approach done"));
        // At the approach prompt: any answer but Enter, and end of input, abort; nothing moves
        for (const char* input : {"y\nq\n", "y\nn\n", "y\nno\n", "y\nabort\n", "y\nstop\n", "y\n \n", "y\n q\n", "y\n"}) {
            CAPTURE(input);
            const Run r = run("0.0.0.0 hold " + logArg("op2"), "", input);
            CHECK(r.code == 6);
            CHECK(r.contains("Aborted"));
            CHECK_FALSE(r.contains("Approach done"));
            CHECK_FALSE(fs::exists(logPath("op2")));
        }
        // At the torque prompt: the same; the arm stays at the start, no torque loop, no log
        for (const char* input : {"y\n\nq\n", "y\n\nno\n", "y\n\nn\n", "y\n\n"}) {
            CAPTURE(input);
            const Run r = run("0.0.0.0 hold " + logArg("op3"), "", input);
            CHECK(r.code == 6);
            CHECK(r.contains("Approach done"));
            CHECK_FALSE(r.contains("Run summary"));
            CHECK_FALSE(fs::exists(logPath("op3")));
        }
        // Enter at every prompt runs it (the operator's normal session)
        const Run ok = run("0.0.0.0 self-motion " + logArg("op4"), "", "y\n\n\n");
        CHECK(ok.code == 0);
        CHECK(ok.contains("Finished at rest."));
    }
    SECTION("the log path is checked before anything moves (exit 9), and a log is never overwritten unasked")
    {
        const std::string missing = (workDir() / "no_such_dir" / "x.csv").string();
        const Run r = run("0.0.0.0 self-motion --yes --log=" + missing);
        INFO(r.out);
        CHECK(r.code == 9);
        CHECK(r.contains("Cannot write"));
        CHECK(r.contains("Nothing moved"));
        CHECK_FALSE(r.contains("Kinematics check"));
        CHECK_FALSE(r.contains("The robot will MOVE"));
        // An existing log
        {
            std::ofstream f(logPath("exists"));
            f << "keep me\n";
        }
        const Run e = run("0.0.0.0 self-motion --yes " + logArg("exists"));
        CHECK(e.code == 9);
        CHECK(e.contains("exists already"));
        CHECK_FALSE(e.contains("The robot will MOVE"));
        std::ifstream kept(logPath("exists"));
        std::string line;
        std::getline(kept, line);
        CHECK(line == "keep me");
        // ... unless --overwrite (a dry run here)
        CHECK(run("0.0.0.0 self-motion --yes --dry-run --overwrite " + logArg("exists")).code == 0);
        // The check leaves no file behind when the run is aborted
        const Run a = run("0.0.0.0 hold " + logArg("aborted"), "", "n\n");
        CHECK(a.code == 6);
        CHECK_FALSE(fs::exists(logPath("aborted")));
    }
    SECTION("a dry run stops after the preflight: no motion, a prediction, no log")
    {
        const Run r = run("0.0.0.0 paper --dry-run --yes " + logArg("dry"));
        CHECK(r.code == 0);
        CHECK(r.contains("Dry run: no motion."));
        CHECK_FALSE(r.contains("Approach"));
        CHECK_FALSE(fs::exists(logPath("dry")));
        CHECK(fs::exists(fs::path(logPath("dry").string() + ".preflight.csv")));
    }
}
