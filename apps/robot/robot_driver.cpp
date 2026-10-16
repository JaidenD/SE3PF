// robot_driver: the hardware experiments on a Franka Emika Panda (FER) with libfranka 0.9.2 (docs/HARDWARE.md).
//
// The only code of the project that includes libfranka. The same source builds two programs:
//   robot_driver       against the real libfranka 0.9.2 (cmake -DSE3PF_WITH_LIBFRANKA=ON -DFranka_DIR=...);
//   robot_driver_mock  against mock/franka (always built with the tests; software in the loop), with the test-only
//                      options --yes (answers every prompt) and --inject (forces a violation), compiled only into
//                      the mock build (SE3PF_MOCK_FRANKA).
//
// Usage:
//   robot_driver list                                          the experiments, in the order to run them
//   robot_driver <ip> check [--no-rt]                          connect, read the state, print q, the configured
//                                                              load and the realtime setup, kinematics check; no
//                                                              motion
//   robot_driver <ip> preflight [experiment|all] [--log-dir=DIR] [--no-rt]
//                                                              kinematics check and the preflight with the robot's
//                                                              own model (load included); no motion
//   robot_driver <ip> error-recovery                           automaticErrorRecovery() after a reflex; no motion
//                                                              (not to be confused with the experiment 'recover')
//   robot_driver <ip> <experiment> [--dry-run] [--log=file.csv] [--overwrite]
//                                                              one experiment (below); --no-rt only with --dry-run
//
// A run, in this order; each step must pass, or nothing (more) moves:
//   1. the log: --log must not exist (unless --overwrite), and it and <log>.preflight.csv must be writable;
//   2. the kinematics check (PandaArm against franka::Model, flange frame, 1e-5);
//   3. the configured end-effector load, which the operator confirms (y);
//   4. the preflight with franka::Model and that load (the prediction goes to <log>.preflight.csv); --dry-run
//      stops here;
//   5. a PREEMPT_RT kernel (real build); the robot Idle, the arm at rest and within 1 rad of the start; then the
//      collision behavior and impedances of the libfranka examples (libfranka accepts them only when Idle);
//   6. Enter (an empty line; any other answer, and end of input, aborts); the state is read again and step 5's checks
//      repeated (and the arm must not have moved since the prompt); a minimum-jerk joint motion to the start
//      configuration (franka::JointPositions, 20 % / 10 % / 5 % of the velocity / acceleration / jerk limits, at
//      least 3 s);
//   7. the arm at the start: Idle, at rest, joints within 0.01 rad of q0, xi and beta within 2e-3 of the preflight's
//      start (printed next to their bounds); and the preflight again, from the measured start configuration;
//   8. Enter; the state is read again and step 7 repeated; the 1 kHz torque loop: ExperimentLoop::follow (the
//      preflight's cycle: Controller::step, then the safety monitor's checks), tau = M(q) qdd* + C(q, qd) qd from
//      franka::Model (the robot adds G), sent with libfranka's rate limiter and low-pass filter.
// On any violation, or a failed step, the loop switches to the safety stop (SafetyStop, se3pf/experiments.hpp): a
// damped hold at the configuration of the violation, tau = -K_s (q - q_stop) - D qd, until the arm has been at rest
// (|qd| < 0.01 rad/s) for 0.1 s, or for at most 2 s; then the loop ends and the robot's own controller takes the arm.
// The log (one row per callback, buffer allocated before the loop) is written after the loop, also after a reflex.
//
// Exit codes:
//   0  done (or check / preflight / dry run passed)
//   1  usage error (also --no-rt for a run that moves the robot)
//   2  unknown experiment
//   3  kinematics check failed (nothing moved)
//   4  preflight failed: before the approach (nothing moved), or from the measured start after it (the arm stays at
//      the start; no torque control)
//   5  safety stop: the monitor stopped the torque loop and the damped hold brought the arm to rest
//   6  aborted by the operator (any answer but Enter alone, end of input, or the load not confirmed)
//   7  not ready: no realtime kernel, robot mode not Idle, the arm moving, too far from the start, moved while
//      waiting, or not at the start after the approach
//   8  libfranka error (connection, reflex, user stop, realtime, ...): see the message
//   9  a log file: --log exists or cannot be written (checked first: nothing moved), or the run's log could not be
//      written after the run
#include <franka/control_types.h>
#include <franka/duration.h>
#include <franka/exception.h>
#include <franka/log.h>
#include <franka/lowpass_filter.h>
#include <franka/model.h>
#include <franka/robot.h>
#include <franka/robot_state.h>

#include "franka_convert.hpp"
#include "robot_log.hpp"
#include "se3pf/experiments.hpp"
#include "se3pf/panda.hpp"
#include "se3pf/panda_example.hpp"
#include "se3pf/preflight.hpp"
#include "se3pf/safety.hpp"

#include <sys/resource.h>  // getrlimit (RLIMIT_RTPRIO)
#include <sys/utsname.h>   // uname
#include <termios.h>       // tcflush
#include <unistd.h>        // isatty, STDIN_FILENO

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#ifdef SE3PF_MOCK_FRANKA
#include <atomic>
#include <cstdlib>
#include <new>
#endif

using namespace se3pf;
using namespace se3pf::robot;

#ifdef SE3PF_MOCK_FRANKA
// Mock build only: count the allocations made inside the 1 kHz callback (Rules for the 1 kHz loop), by replacing
// the global operator new of this program. The summary prints the count and the software-in-the-loop test requires 0.
namespace {
    std::atomic<long> gCallbackAllocations{0};
    thread_local bool gInCallback = false;
}  // namespace

void* operator new(std::size_t size)
{
    if (gInCallback)
        gCallbackAllocations.fetch_add(1, std::memory_order_relaxed);
    if (void* p = std::malloc(size == 0 ? 1 : size))
        return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
#define SE3PF_CALLBACK_BEGIN() (gInCallback = true)
#define SE3PF_CALLBACK_END() (gInCallback = false)
#else
#define SE3PF_CALLBACK_BEGIN() ((void)0)
#define SE3PF_CALLBACK_END() ((void)0)
#endif

namespace {

    // ---------------------------------------------------------------------------------------------------------
    // Settings

    enum Exit : int {
        kOk = 0,
        kUsage = 1,
        kUnknownExperiment = 2,
        kKinematics = 3,
        kPreflight = 4,
        kSafetyStop = 5,
        kAborted = 6,
        kNotReady = 7,
        kRobotError = 8,
        kLogError = 9
    };

    constexpr double kKinematicsTolerance = 1e-5;  // PandaArm against franka::Model, pose and J_b entries
    constexpr double kMaxApproachDistance = 1.0;   // rad, any joint: farther, guide the arm closer in Desk first
    constexpr double kStartTolerance = 0.01;       // rad, any joint, after the approach
    // |xi - xi0| and |beta - beta0| after the approach (5 % of the 0.04 tube bound, 4 % of the 0.05 beta bound; the
    // mock's approach ends 1e-5 off, the robot's should be about 1e-4). A sanity check: the margins themselves are
    // checked by a second preflight from the measured start.
    constexpr double kStartOutputTolerance = 2e-3;
    constexpr double kRestSpeed = 0.02;            // rad/s: the arm must be at rest before a motion

    struct Options {
        std::string ip;
        std::string command;       // check, preflight, error-recovery, or an experiment name
        std::string which = "all"; // preflight: one experiment or all
        std::string logPath;       // run: the log file (default robot_<name>_<time>.csv)
        std::string logDir;        // preflight: write the predictions there
        bool dryRun = false;
        bool noRt = false;
        bool overwrite = false;    // run: an existing --log may be replaced
        bool yes = false;          // mock build only: answer every prompt
        Violation inject = Violation::None;  // mock build only: force this violation from injectTime on ...
        double injectTime = 0.0;             // ... s after the torque loop started
    };

    void usage(const char* prog)
    {
        std::cerr << "usage: " << prog << " list\n"
                  << "       " << prog << " <robot-ip> check [--no-rt]\n"
                  << "       " << prog << " <robot-ip> preflight [experiment|all] [--log-dir=DIR] [--no-rt]\n"
                  << "       " << prog << " <robot-ip> error-recovery\n"
                  << "       " << prog << " <robot-ip> <experiment> [--dry-run [--no-rt]] [--log=file.csv] "
                                          "[--overwrite]\n"
#ifdef SE3PF_MOCK_FRANKA
                  << "  (mock build) --yes answers every prompt; --inject=<violation>@<t> forces that violation t s "
                     "into the torque loop (violations: projection, joint-limit, joint-velocity, acceleration, "
                     "torque, torque-rate, flange-speed, tube, redundancy, conditioning, non-finite, timeout)\n"
#endif
                  << "exit codes: 0 ok, 1 usage, 2 unknown experiment, 3 kinematics check, 4 preflight, 5 safety "
                     "stop, 6 aborted, 7 not ready / not at the start, 8 libfranka error, 9 log file\n";
    }

    // ---------------------------------------------------------------------------------------------------------
    // The operator

    // Discards whatever was typed before a prompt, so that an Enter pressed during the approach cannot start the
    // torque loop. Nothing to do when standard input is not a terminal (the software-in-the-loop tests).
    void discardTypedInput()
    {
        if (::isatty(STDIN_FILENO))
            ::tcflush(STDIN_FILENO, TCIFLUSH);
    }

    // Only Enter alone (an empty line) continues; any other answer ('n', 'no', 'q', 'stop', ' ') and end of input
    // abort.
    bool confirm(const Options& opt, const std::string& what)
    {
        std::cout << what << "\n  Press Enter to continue; anything else, then Enter, aborts: " << std::flush;
#ifdef SE3PF_MOCK_FRANKA
        if (opt.yes) {
            std::cout << "(--yes)\n";
            return true;
        }
#endif
        (void)opt;
        discardTypedInput();
        std::string line;
        if (!std::getline(std::cin, line)) {
            std::cout << "\n  Aborted (end of input).\n";
            return false;
        }
        if (!line.empty()) {
            std::cout << "  Aborted (the answer was not Enter alone).\n";
            return false;
        }
        return true;
    }

    // Only "y" or "yes" confirms.
    bool askYes(const Options& opt, const std::string& question)
    {
        std::cout << question << " [y/N]: " << std::flush;
#ifdef SE3PF_MOCK_FRANKA
        if (opt.yes) {
            std::cout << "(--yes)\n";
            return true;
        }
#endif
        (void)opt;
        discardTypedInput();
        std::string line;
        if (!std::getline(std::cin, line))
            return false;
        return line == "y" || line == "Y" || line == "yes";
    }

    std::string defaultLogPath(const std::string& name)
    {
        char stamp[32];
        const std::time_t now = std::time(nullptr);
        std::strftime(stamp, sizeof(stamp), "%Y%m%d-%H%M%S", std::localtime(&now));
        return "robot_" + name + "_" + stamp + ".csv";
    }

    std::string vecString(const Vec7& v, int precision = 4)
    {
        std::string s = "(";
        char buf[32];
        for (int i = 0; i < 7; ++i) {
            std::snprintf(buf, sizeof(buf), "%s%.*f", i ? ", " : "", precision, v[i]);
            s += buf;
        }
        return s + ")";
    }

    const char* modeName(franka::RobotMode m)
    {
        switch (m) {
            case franka::RobotMode::kOther: return "Other (joints locked?)";
            case franka::RobotMode::kIdle: return "Idle";
            case franka::RobotMode::kMove: return "Move";
            case franka::RobotMode::kGuiding: return "Guiding";
            case franka::RobotMode::kReflex: return "Reflex";
            case franka::RobotMode::kUserStopped: return "UserStopped (user stop pressed)";
            case franka::RobotMode::kAutomaticErrorRecovery: return "AutomaticErrorRecovery";
        }
        return "unknown";
    }

#ifdef SE3PF_MOCK_FRANKA
    // --inject's names of the violations
    Violation violationFromName(const std::string& name)
    {
        static const std::pair<const char*, Violation> names[] = {
            {"projection", Violation::Projection},     {"joint-limit", Violation::JointLimit},
            {"joint-velocity", Violation::JointVelocity}, {"acceleration", Violation::Acceleration},
            {"torque", Violation::Torque},             {"torque-rate", Violation::TorqueRate},
            {"flange-speed", Violation::FlangeSpeed},  {"tube", Violation::Tube},
            {"redundancy", Violation::Redundancy},     {"conditioning", Violation::Conditioning},
            {"non-finite", Violation::NonFinite},      {"timeout", Violation::Timeout}};
        for (const auto& n : names)
            if (name == n.first)
                return n.second;
        return Violation::None;
    }
#endif

    // ---------------------------------------------------------------------------------------------------------
    // The lab PC

    // A PREEMPT_RT kernel: /sys/kernel/realtime reads 1 (the RT patch), or uname -v names PREEMPT_RT (or "PREEMPT
    // RT"). libfranka 0.9.2 reads only the file, and reads an uninitialized value when it is missing
    // (src/control_tools.cpp), so the driver checks this itself before any motion. detail says what was found,
    // with the realtime priority limit (ulimit -r: libfranka needs 99, or it refuses to start a motion).
    bool realtimeKernel(std::string& detail)
    {
        int flag = 0;
        std::ifstream f("/sys/kernel/realtime");
        const bool haveFile = static_cast<bool>(f >> flag);
        struct utsname u {};
        const bool named = ::uname(&u) == 0
                           && (std::strstr(u.version, "PREEMPT_RT") != nullptr
                               || std::strstr(u.version, "PREEMPT RT") != nullptr);
        struct rlimit rl {};
        char rtprio[32] = "unknown";
        if (::getrlimit(RLIMIT_RTPRIO, &rl) == 0) {
            if (rl.rlim_cur == RLIM_INFINITY)
                std::snprintf(rtprio, sizeof(rtprio), "unlimited");
            else
                std::snprintf(rtprio, sizeof(rtprio), "%lu", static_cast<unsigned long>(rl.rlim_cur));
        }
        char buf[200];
        std::snprintf(buf, sizeof(buf), "/sys/kernel/realtime %s; uname -v %s PREEMPT_RT; ulimit -r %s (needs 99)",
                      haveFile ? (flag == 1 ? "= 1" : "is not 1") : "missing", named ? "names" : "does not name",
                      rtprio);
        detail = buf;
        return (haveFile && flag == 1) || named;
    }

    // The run's log and its prediction, before anything moves: an existing log is refused unless --overwrite, and
    // both files must be writable (each is opened for appending, and removed again if it did not exist).
    bool checkLogPath(const std::string& path, bool overwrite)
    {
        std::error_code ec;
        if (!overwrite && std::filesystem::exists(path, ec)) {
            std::printf("The log %s exists already: choose another --log, or add --overwrite. Nothing moved.\n",
                        path.c_str());
            return false;
        }
        for (const std::string& p : {path, path + ".preflight.csv"}) {
            const bool existed = std::filesystem::exists(p, ec);
            bool writable = false;
            {
                std::ofstream f(p, std::ios::app);
                writable = static_cast<bool>(f);
            }
            if (!writable) {
                std::printf("Cannot write %s (does its directory exist?). Nothing moved.\n", p.c_str());
                return false;
            }
            if (!existed)
                std::filesystem::remove(p, ec);
        }
        return true;
    }

    // ---------------------------------------------------------------------------------------------------------
    // The robot's model

    // franka::Model as a DynamicsModel for the preflight, with the end-effector load of the robot state (m_total,
    // F_x_Ctotal, I_total: the end effector configured in Desk plus any load set by setLoad), as the RobotState
    // overloads of mass, coriolis and gravity use it in the torque loop.
    class FrankaDynamics final : public DynamicsModel {
        public:
            FrankaDynamics(const franka::Model& model, const franka::RobotState& state)
                : model_(model), I_total_(state.I_total), m_total_(state.m_total), F_x_Ctotal_(state.F_x_Ctotal)
            {
            }
            Mat<7, 7> mass(const Vec7& q) const override
            {
                return mat7FromColumnMajor(model_.mass(array7(q), I_total_, m_total_, F_x_Ctotal_));
            }
            Vec7 coriolis(const Vec7& q, const Vec7& dq) const override
            {
                return vec7(model_.coriolis(array7(q), array7(dq), I_total_, m_total_, F_x_Ctotal_));
            }
            Vec7 gravity(const Vec7& q) const override { return vec7(model_.gravity(array7(q), m_total_, F_x_Ctotal_)); }

        private:
            const franka::Model& model_;
            std::array<double, 9> I_total_;
            double m_total_;
            std::array<double, 3> F_x_Ctotal_;
    };

    // Collision thresholds of libfranka's torque-control example (joint_impedance_control.cpp) and the joint and
    // Cartesian impedances of examples_common.cpp (they act during the joint position approach). libfranka rejects
    // them unless the robot is Idle (CommandException), so they are set after the Idle check.
    void setDefaultBehavior(franka::Robot& robot)
    {
        robot.setCollisionBehavior({{20.0, 20.0, 18.0, 18.0, 16.0, 14.0, 12.0}},
                                   {{20.0, 20.0, 18.0, 18.0, 16.0, 14.0, 12.0}},
                                   {{20.0, 20.0, 18.0, 18.0, 16.0, 14.0, 12.0}},
                                   {{20.0, 20.0, 18.0, 18.0, 16.0, 14.0, 12.0}},
                                   {{20.0, 20.0, 20.0, 25.0, 25.0, 25.0}}, {{20.0, 20.0, 20.0, 25.0, 25.0, 25.0}},
                                   {{20.0, 20.0, 20.0, 25.0, 25.0, 25.0}}, {{20.0, 20.0, 20.0, 25.0, 25.0, 25.0}});
        robot.setJointImpedance({{3000, 3000, 3000, 2500, 2500, 2000, 2000}});
        robot.setCartesianImpedance({{3000, 3000, 3000, 300, 300, 300}});
    }

    struct NamedConfiguration {
        std::string name;
        Vec7 q;
    };

    // The configurations of the kinematics check (besides the robot's own q): the rest posture, every experiment's
    // start and eight points of the lift q*(s).
    std::vector<NamedConfiguration> checkConfigurations(const PandaExample& example)
    {
        std::vector<NamedConfiguration> qs;
        qs.push_back({"the rest posture", PandaArm::restPosture()});
        for (const ExperimentSpec& e : experiments())
            qs.push_back({"the start of " + e.name, ExperimentSetup(example, e).start()});
        for (int k = 0; k < 8; ++k) {
            Vec7 q = Vec7::Zero();
            const double s = k * example.tube().length() / 8.0;
            example.lift().eval(s, q);
            char name[64];
            std::snprintf(name, sizeof(name), "the lift at s = %.3f", s);
            qs.push_back({name, q});
        }
        return qs;
    }

    // The largest of one kind of error over the configurations, and where. A NaN sticks (and fails the check).
    struct WorstError {
        double value = 0.0;
        std::size_t where = 0;
        void add(double e, std::size_t k)
        {
            if (std::isnan(value))
                return;
            if (!(e <= value)) {
                value = e;
                where = k;
            }
        }
    };

    // PandaArm::output against model.pose(Frame::kFlange, ...) and model.bodyJacobian(Frame::kFlange, ...): the
    // flange (Franka frame F), not kEndEffector (which includes the end-effector offset set in Desk), and the body
    // Jacobian with libfranka's (v, omega) rows swapped to (omega, v). First with the RobotState overloads at the
    // robot's q, then with the q overloads at the check configurations. Fails closed (NaN fails). Prints the position
    // and rotation errors of the pose and the omega and v blocks of J_b separately, each with its worst
    // configuration, and J_b without the row swap, so that a convention error can be told from a small numeric
    // difference.
    bool checkKinematics(const franka::Model& model, const franka::RobotState& state,
                         const std::vector<NamedConfiguration>& configurations)
    {
        std::vector<std::string> names{"the robot's q"};
        for (const NamedConfiguration& c : configurations)
            names.push_back(c.name);
        WorstError pose, position, rotation, omega, v, unswapped;
        auto compare = [&](std::size_t k, const Vec7& q, const std::array<double, 16>& poseArray,
                           const std::array<double, 42>& jac) {
            Mat67 Jb = Mat67::Zero();
            const Pose g = PandaArm::output(q, &Jb);
            const Pose gf = poseFromColumnMajor(poseArray);
            const Mat67 Jf = bodyJacobianFromFranka(jac);
            Mat67 Jraw = Mat67::Zero();  // libfranka's rows as they come, (v, omega)
            for (int c = 0; c < 7; ++c)
                for (int r = 0; r < 6; ++r)
                    Jraw(r, c) = jac[static_cast<std::size_t>(c * 6 + r)];
            pose.add((g - gf).cwiseAbs().maxCoeff(), k);
            position.add((g.block<3, 1>(0, 3) - gf.block<3, 1>(0, 3)).cwiseAbs().maxCoeff(), k);
            rotation.add((g.block<3, 3>(0, 0) - gf.block<3, 3>(0, 0)).cwiseAbs().maxCoeff(), k);
            omega.add((Jb.topRows<3>() - Jf.topRows<3>()).cwiseAbs().maxCoeff(), k);
            v.add((Jb.bottomRows<3>() - Jf.bottomRows<3>()).cwiseAbs().maxCoeff(), k);
            unswapped.add((Jb - Jraw).cwiseAbs().maxCoeff(), k);
        };
        compare(0, vec7(state.q), model.pose(franka::Frame::kFlange, state),
                model.bodyJacobian(franka::Frame::kFlange, state));
        for (std::size_t k = 0; k < configurations.size(); ++k) {
            const Vec7& q = configurations[k].q;
            compare(k + 1, q, model.pose(franka::Frame::kFlange, array7(q), state.F_T_EE, state.EE_T_K),
                    model.bodyJacobian(franka::Frame::kFlange, array7(q), state.F_T_EE, state.EE_T_K));
        }
        const double jacobian = std::isnan(omega.value) || std::isnan(v.value) ? std::nan("")
                                                                                 : std::max(omega.value, v.value);
        const bool ok = pose.value <= kKinematicsTolerance && jacobian <= kKinematicsTolerance;
        std::printf("Kinematics check (PandaArm vs franka::Model at the flange, %zu configurations, tolerance %.0e): "
                    "pose %.1e, J_b %.1e  %s\n",
                    names.size(), kKinematicsTolerance, pose.value, jacobian,
                    ok ? "OK" : "MISMATCH: do not run anything.");
        std::printf("  pose: position %.1e m (worst at %s), rotation %.1e (worst at %s)\n", position.value,
                    names[position.where].c_str(), rotation.value, names[rotation.where].c_str());
        std::printf("  J_b, rows (omega, v): omega %.1e (worst at %s), v %.1e (worst at %s); without swapping "
                    "libfranka's (v, omega) rows: %.1e\n",
                    omega.value, names[omega.where].c_str(), v.value, names[v.where].c_str(), unswapped.value);
        if (!ok)
            std::printf("  A difference of order 1 (or J_b fine only without the row swap) is a convention error: the "
                        "robot (FER, not FR3), libfranka 0.9.2, or a change of PandaArm (src/panda.cpp). A small one "
                        "(1e-5 to 1e-3) is numeric: bring this output back, do not loosen the tolerance.\n");
        return ok;
    }

    void printState(const franka::RobotState& s)
    {
        std::printf("Robot mode: %s\n", modeName(s.robot_mode));
        std::printf("q = %s rad\n", vecString(vec7(s.q)).c_str());
        std::printf("|dq| = %.4f rad/s, joint margin %.3f rad\n", vec7(s.dq).norm(), PandaArm::jointMargin(vec7(s.q)));
    }

    // The end-effector load the robot's model (and so G, M and C) includes, as configured in Desk
    void printLoad(const franka::RobotState& s)
    {
        std::printf("Configured load (Desk end effector + setLoad), in the flange frame:\n");
        std::printf("  end effector: m_ee = %.3f kg, centre of mass F_x_Cee = (%.4f, %.4f, %.4f) m\n", s.m_ee,
                    s.F_x_Cee[0], s.F_x_Cee[1], s.F_x_Cee[2]);
        std::printf("  extra load:   m_load = %.3f kg, F_x_Cload = (%.4f, %.4f, %.4f) m\n", s.m_load, s.F_x_Cload[0],
                    s.F_x_Cload[1], s.F_x_Cload[2]);
        const Matrix3d I = mat3FromColumnMajor(s.I_total);
        std::printf("  total:        m_total = %.3f kg, F_x_Ctotal = (%.4f, %.4f, %.4f) m, I_total diag = (%.4f, "
                    "%.4f, %.4f) kg m^2\n",
                    s.m_total, s.F_x_Ctotal[0], s.F_x_Ctotal[1], s.F_x_Ctotal[2], I(0, 0), I(1, 1), I(2, 2));
        std::printf("  (Franka Hand: about 0.73 kg; no end effector: 0 kg. An end effector missing here makes the arm "
                    "sag: MuJoCo run D.)\n");
    }

    // Runs the preflight of one experiment with the robot's model; prints the report. log (optional) receives the
    // prediction.
    PreflightReport runPreflight(const ExperimentSetup& setup, const FrankaDynamics& dyn,
                                 std::vector<ExperimentLogRow>* log)
    {
        const PreflightReport rep = preflight(setup, &dyn, setup.monitor(), log);
        std::printf("Preflight [%s] with the robot's model: %s\n", setup.spec().name.c_str(), rep.summary().c_str());
        return rep;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Readiness

    // Before the approach: the robot Idle, the arm at rest and within kMaxApproachDistance of the start. Prints why
    // not.
    bool readyToApproach(const franka::RobotState& s, const Vec7& start)
    {
        if (s.robot_mode != franka::RobotMode::kIdle) {
            std::printf("Not ready: the robot is in mode %s, not Idle. Release the user stop, unlock the joints and "
                        "activate FCI in Desk, or recover from a reflex (robot_driver <ip> error-recovery).\n",
                        modeName(s.robot_mode));
            return false;
        }
        const double speed = vec7(s.dq).norm();
        if (!(speed <= kRestSpeed)) {
            std::printf("Not ready: the arm is moving (|dq| = %.3f rad/s).\n", speed);
            return false;
        }
        const double distance = (start - vec7(s.q)).cwiseAbs().maxCoeff();
        if (!(distance <= kMaxApproachDistance)) {
            std::printf("Not ready: joint distance to the start %.2f rad > %.1f rad. Guide the arm by hand close to "
                        "the start (near the rest posture %s) and run again.\n",
                        distance, kMaxApproachDistance, vecString(PandaArm::restPosture(), 2).c_str());
            return false;
        }
        return true;
    }

    // After the approach: the robot Idle, the arm at rest, the joints within kStartTolerance of q0 and the outputs
    // within kStartOutputTolerance of the preflight's start. verbose prints the measured start next to the bounds.
    bool atStart(const franka::RobotState& s, const ExperimentSetup& setup, bool verbose)
    {
        const ExperimentSpec& spec = setup.spec();
        const SafetyLimits& b = setup.monitor().limits();
        const ExperimentLoop probe(setup, setup.monitor());
        ExperimentLoop::Status st;
        const Vec7 q = vec7(s.q), dq = vec7(s.dq);
        const double qError = (q - setup.start()).cwiseAbs().maxCoeff();
        const bool outputsOk = probe.observe(q, dq, st);
        const double xiError = outputsOk ? (st.xi - spec.xi0).norm() : std::nan("");
        const double betaError = outputsOk ? std::abs(st.beta - spec.beta0) : std::nan("");
        const double speed = dq.norm();
        if (verbose) {
            std::printf("At the start: |xi| = %.4f [bound %.3f], |beta| = %.4f rad [bound %.3f]\n",
                        outputsOk ? st.xi.norm() : std::nan(""), b.xi, outputsOk ? std::abs(st.beta) : std::nan(""),
                        b.beta);
            std::printf("  |xi - xi0| = %.1e, |beta - beta0| = %.1e [each <= %.0e]; max |q - q0| = %.1e rad [<= "
                        "%.2f]; |dq| = %.1e rad/s [<= %.2f]; mode %s\n",
                        xiError, betaError, kStartOutputTolerance, qError, kStartTolerance, speed, kRestSpeed,
                        modeName(s.robot_mode));
        }
        std::string why;
        if (!(qError <= kStartTolerance))
            why += " joints off by more than the tolerance;";
        if (!(xiError <= kStartOutputTolerance))
            why += " xi off;";
        if (!(betaError <= kStartOutputTolerance))
            why += " beta off;";
        if (!(speed <= kRestSpeed))
            why += " the arm is moving;";
        if (s.robot_mode != franka::RobotMode::kIdle)
            why += std::string(" robot mode ") + modeName(s.robot_mode) + ", not Idle;";
        if (!why.empty())
            std::printf("Not ready: not at rest at the start:%s no torque control. Nothing moves; run again (the "
                        "approach starts from where the arm is).\n",
                        why.c_str());
        return why.empty();
    }

    // ---------------------------------------------------------------------------------------------------------
    // Motions

    // Minimum-jerk joint motion from the robot's q_d at the first callback to goal (ApproachMotion). False (and no
    // motion) if the arm is farther than kMaxApproachDistance from the goal at the first callback.
    bool approach(franka::Robot& robot, const Vec7& goal, double& duration)
    {
        std::optional<ApproachMotion> motion;
        bool refused = false;
        double t = 0.0;
        std::function<franka::JointPositions(const franka::RobotState&, franka::Duration)> callback =
            [&](const franka::RobotState& s, franka::Duration period) -> franka::JointPositions {
            if (!motion) {
                // Plan from q_d, the robot's last commanded position: the first command must equal it.
                const Vec7 qd = vec7(s.q_d);
                if (!((goal - qd).cwiseAbs().maxCoeff() <= kMaxApproachDistance + 0.05)) {
                    refused = true;
                    return franka::MotionFinished(franka::JointPositions(s.q_d));
                }
                motion.emplace(qd, goal);
                t = 0.0;
            } else {
                t += period.toSec();
            }
            const franka::JointPositions cmd(array7(motion->position(t)));
            return t >= motion->duration() ? franka::MotionFinished(cmd) : cmd;
        };
        robot.control(callback, franka::ControllerMode::kJointImpedance, true, franka::kDefaultCutoffFrequency);
        duration = motion ? motion->duration() : 0.0;
        return !refused;
    }

    struct FollowResult {
        Violation violation = Violation::None;
        double tStop = -1.0;        // time of the violation (the safety stop starts), s
        double tEnd = 0.0;          // time of the last callback, s
        bool finished = false;      // the run ended normally (at rest after the reference)
        bool stopSettled = false;   // the safety stop ended with the arm at rest
        double maxStopDrift = 0.0;  // max_i |q_i - q_stop,i| during the safety stop, rad
        std::size_t rows = 0;       // log rows written
    };

    // The 1 kHz torque loop. Allocation free: everything it uses is built before robot.control, and the log is a
    // buffer of rows allocated before the loop. res is updated at every callback, so it also tells what happened
    // before an exception (a reflex) ended robot.control.
    void follow(franka::Robot& robot, const franka::Model& model, ExperimentLoop& loop, const Options& opt,
                std::vector<ExperimentLogRow>& log, FollowResult& res)
    {
        res = FollowResult();
        double t = 0.0;
        bool stopping = false;
        SafetyStop stop;
        ExperimentLoop::Cycle cycle;
        ModelTerms terms;
        (void)opt;  // the real build has no test options
        std::function<franka::Torques(const franka::RobotState&, franka::Duration)> callback =
            [&](const franka::RobotState& s, franka::Duration period) -> franka::Torques {
            SE3PF_CALLBACK_BEGIN();
            const auto start = std::chrono::steady_clock::now();
            t += period.toSec();
            const Vec7 q = vec7(s.q), dq = vec7(s.dq);
            Vec7 tau = Vec7::Zero(), tauG = Vec7::Zero(), qdd = Vec7::Zero();
            ExperimentLoop::Status st;
            bool done = false;
            bool violatedNow = false;

            if (!stopping) {
                // eq:fb with the robot's model at the measured state; the RobotState overloads include the load
                terms.M = mat7FromColumnMajor(model.mass(s));
                terms.c = vec7(model.coriolis(s));
                terms.G = vec7(model.gravity(s));
                loop.follow(t, q, dq, &terms, cycle);
#ifdef SE3PF_MOCK_FRANKA
                if (opt.inject != Violation::None && t >= opt.injectTime && cycle.violation == Violation::None)
                    cycle.violation = opt.inject;  // test only: the stop path for every kind of violation
#endif
                st = cycle.status;
                qdd = cycle.qdd;
                if (cycle.violation != Violation::None) {
                    stopping = true;
                    violatedNow = true;
                    res.violation = cycle.violation;
                    res.tStop = t;
                    stop.begin(t, q);
                } else {
                    tau = cycle.tau;
                    tauG = cycle.tauG;
                    done = cycle.finished;
                    res.finished = done;
                }
            }
            if (stopping) {
                // The safety stop: the damped hold at q_stop (SafetyStop); the robot adds G
                tau = stop.step(t, q, dq, done);
                tauG = tau + vec7(model.gravity(s));
                const double drift = (q - stop.holdConfiguration()).cwiseAbs().maxCoeff();
                if (!std::isnan(res.maxStopDrift) && !(drift <= res.maxStopDrift))
                    res.maxStopDrift = drift;
                res.stopSettled = stop.settled();
                if (!violatedNow && !loop.observe(q, dq, st)) {
                    st.xi = Vec5::Constant(std::nan(""));
                    st.beta = std::nan("");
                }
                st.t = t;
            }
            res.tEnd = t;
            if (res.rows < log.size()) {
                ExperimentLogRow& row = log[res.rows++];
                row = makeLogRow(t, stopping ? 1 : 0, st, q, dq, qdd, tau, tauG);
                row.periodMs = static_cast<double>(period.toMSec());
                row.successRate = s.control_command_success_rate;
                row.computeUs = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - start)
                                    .count();
            }
            const franka::Torques command(array7(tau));
            SE3PF_CALLBACK_END();
            return done ? franka::MotionFinished(command) : command;
        };
        robot.control(callback, true, franka::kDefaultCutoffFrequency);
    }

    // The summary of a run's log
    void printRunSummary(const std::vector<ExperimentLogRow>& log, std::size_t rows, const FollowResult& res,
                         const ExperimentSetup& setup, bool reflex)
    {
        RunPeaks peaks;
        peaks.torques = true;
        double maxCompute = 0.0, minSuccess = 1.0;
        long late = 0;
        for (std::size_t k = 0; k < rows; ++k) {
            peaks.add(log[k]);
            maxCompute = std::max(maxCompute, log[k].computeUs);
            if (k > 0 && log[k].periodMs > 1.0)
                ++late;
            minSuccess = std::min(minSuccess, log[k].successRate);
        }
        std::printf("Run summary [%s]: %zu callbacks over %.3f s\n", setup.spec().name.c_str(), rows, res.tEnd);
        std::printf("  last following cycle: t = %.3f s, |xi| = %.2e, beta = %.2e rad, tangential error = %.2e, "
                    "eta_dot = %.4f\n",
                    peaks.finalT, peaks.finalXi, peaks.finalBeta, peaks.finalETan, peaks.finalEtaDot);
        std::printf("  peaks: %s\n", peaks.summary(setup.monitor().limits()).c_str());
        std::printf("  timing: callback at most %.0f us, %ld periods longer than 1 ms, lowest command success rate "
                    "%.3f\n",
                    maxCompute, late, minSuccess);
#ifdef SE3PF_MOCK_FRANKA
        std::printf("  callback allocations: %ld\n", gCallbackAllocations.load());
#endif
        if (reflex) {
            std::printf("ENDED BY A ROBOT REFLEX at t = %.3f s (see above)%s\n", res.tEnd,
                        res.violation != Violation::None ? ", after a safety stop" : "");
            if (res.violation != Violation::None)
                std::printf("  safety stop: %s at t = %.3f s\n", toString(res.violation), res.tStop);
        } else if (res.violation != Violation::None) {
            const double endSpeed = rows > 0 ? log[rows - 1].dq.norm() : 0.0;
            std::printf("SAFETY STOP: %s at t = %.3f s; held at the stop configuration until t = %.3f s: largest "
                        "joint drift %.4f rad, final |dq| = %.4f rad/s%s\n",
                        toString(res.violation), res.tStop, res.tEnd, res.maxStopDrift, endSpeed,
                        res.stopSettled ? " (at rest)" : "");
            if (!res.stopSettled)
                std::printf("WARNING: the arm was not at rest when the stop ended after %.0f s; the robot's own "
                            "controller has it now. Check the end effector in Desk and the log before anything else "
                            "(docs/HARDWARE.md, section 7).\n",
                            RunConstants::stopTime);
        } else {
            std::printf("Finished at rest.\n");
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    // Commands

    int cmdList()
    {
        int k = 1;
        for (const ExperimentSpec& e : experiments()) {
            std::printf("%d. %-12s %s (%.0f s + time to rest)\n   %s\n   Watch: %s\n", k++, e.name.c_str(),
                        e.title.c_str(), e.duration, e.purpose.c_str(), e.watch.c_str());
        }
        return kOk;
    }

    int cmdCheck(franka::Robot& robot, const PandaExample& example)
    {
        std::printf("Robot server version %u (libfranka 0.9.x speaks version 5)\n",
                    static_cast<unsigned>(robot.serverVersion()));
        std::string rt;
        const bool realtime = realtimeKernel(rt);
        std::printf("This PC: %s realtime kernel (%s)\n", realtime ? "a" : "NO", rt.c_str());
        const franka::RobotState state = robot.readOnce();
        printState(state);
        printLoad(state);
        franka::Model model = robot.loadModel();
        return checkKinematics(model, state, checkConfigurations(example)) ? kOk : kKinematics;
    }

    int cmdPreflight(franka::Robot& robot, const PandaExample& example, const Options& opt)
    {
        const franka::RobotState state = robot.readOnce();
        printLoad(state);
        franka::Model model = robot.loadModel();
        const bool kinematicsOk = checkKinematics(model, state, checkConfigurations(example));
        const FrankaDynamics dyn(model, state);
        bool ok = true;
        for (const ExperimentSpec& e : experiments()) {
            if (opt.which != "all" && opt.which != e.name)
                continue;
            const ExperimentSetup setup(example, e);
            std::vector<ExperimentLogRow> prediction;
            ok = runPreflight(setup, dyn, &prediction).ok && ok;
            if (!opt.logDir.empty()) {
                const std::string path = opt.logDir + "/preflight_" + e.name + ".csv";
                if (!writeLogCsv(path, prediction, prediction.size()))
                    std::printf("  could not write %s\n", path.c_str());
            }
        }
        if (!kinematicsOk)
            return kKinematics;
        return ok ? kOk : kPreflight;
    }

    int cmdRecover(franka::Robot& robot, const Options& opt)
    {
        const franka::RobotState state = robot.readOnce();
        printState(state);
        if (!confirm(opt, "automaticErrorRecovery() clears a reflex. The arm does not move."))
            return kAborted;
        robot.automaticErrorRecovery();
        std::printf("Recovered. Robot mode: %s\n", modeName(robot.readOnce().robot_mode));
        return kOk;
    }

    int cmdRun(franka::Robot& robot, const PandaExample& example, const ExperimentSpec& spec, const Options& opt)
    {
        const ExperimentSetup setup(example, spec);
        const std::string logPath = opt.logPath.empty() ? defaultLogPath(spec.name) : opt.logPath;
        std::printf("Experiment '%s': %s\n  %s\n  Watch: %s\n", spec.name.c_str(), spec.title.c_str(),
                    spec.purpose.c_str(), spec.watch.c_str());
        std::printf("  start q0 = %s rad (|xi0| = %.4f, beta0 = %.3f rad), %.1f s + time to rest\n",
                    vecString(setup.start()).c_str(), spec.xi0.norm(), spec.beta0, spec.duration);

        // 1. The log, before anything else: a run whose log cannot be written must not move the robot
        if (!checkLogPath(logPath, opt.overwrite))
            return kLogError;

        // 2. The kinematics
        franka::Model model = robot.loadModel();
        franka::RobotState state = robot.readOnce();
        printState(state);
        if (!checkKinematics(model, state, checkConfigurations(example)))
            return kKinematics;

        // 3. The load
        printLoad(state);
        if (!askYes(opt, "Does this load match the end effector mounted on the robot right now?")) {
            std::printf("Aborted: set the end effector in Desk (Settings > End Effector), then run again.\n");
            return kAborted;
        }

        // 4. The preflight with the robot's own model and that load
        const FrankaDynamics dyn(model, state);
        std::vector<ExperimentLogRow> prediction;
        const PreflightReport rep = runPreflight(setup, dyn, &prediction);
        const std::string predictionPath = logPath + ".preflight.csv";
        const bool predictionWritten = writeLogCsv(predictionPath, prediction, prediction.size());
        std::printf("  prediction: %s%s\n", predictionPath.c_str(), predictionWritten ? "" : "  COULD NOT BE WRITTEN");
        if (!rep.ok) {
            std::printf("Aborted: the predicted motion violates a margin. Nothing moved. Do not run it: bring the "
                        "printed peaks back (docs/HARDWARE.md, section 7).\n");
            return kPreflight;
        }
        if (!predictionWritten) {
            std::printf("Aborted: the prediction could not be written. Nothing moved.\n");
            return kLogError;
        }

        // 5. The realtime kernel and the robot's state (read again: the preflight and the question can take minutes)
        std::string rtDetail;
        const bool realtime = realtimeKernel(rtDetail);
        state = robot.readOnce();
        if (opt.dryRun) {
            std::printf("Dry run: no motion. For a run: robot mode %s (needs Idle); %s realtime kernel (%s).\n",
                        modeName(state.robot_mode), realtime ? "a" : "NO", rtDetail.c_str());
            return kOk;
        }
#ifndef SE3PF_MOCK_FRANKA
        if (!realtime) {
            std::printf("Not ready: this PC has no realtime kernel (%s). The torque loop needs PREEMPT_RT "
                        "(docs/HARDWARE.md, section 2). Nothing moved.\n",
                        rtDetail.c_str());
            return kNotReady;
        }
#endif
        if (!readyToApproach(state, setup.start()))
            return kNotReady;
        setDefaultBehavior(robot);

        // 6. The approach, after Enter, from a state read again after Enter
        const Vec7 qPrompt = vec7(state.q);
        const ApproachMotion plan(qPrompt, setup.start());
        char what[300];
        std::snprintf(what, sizeof(what),
                      "The robot will MOVE to the start configuration: joint motion of up to %.3f rad, %.1f s "
                      "(20 %% of the velocity limits). Keep the user stop in hand.",
                      plan.maxJointDistance(), plan.duration());
        if (!confirm(opt, what))
            return kAborted;
        state = robot.readOnce();
        if (!readyToApproach(state, setup.start()))
            return kNotReady;
        const double moved = (vec7(state.q) - qPrompt).cwiseAbs().maxCoeff();
        if (!(moved <= kStartTolerance)) {
            std::printf("Not ready: the arm moved by %.3f rad while waiting for Enter. Nothing moved; run again.\n",
                        moved);
            return kNotReady;
        }
        double approachTime = 0.0;
        if (!approach(robot, setup.start(), approachTime)) {
            std::printf("Not ready: the arm moved away from where it was planned from; nothing moved.\n");
            return kNotReady;
        }
        std::printf("Approach done (%.1f s).\n", approachTime);

        // 7. At the start? Then the preflight again, from where the arm actually is
        state = robot.readOnce();
        if (!atStart(state, setup, true))
            return kNotReady;
        {
            const Vec7 qMeasured = vec7(state.q);
            std::vector<ExperimentLogRow> again;
            const PreflightReport from = preflight(setup, &dyn, setup.monitor(), &again, &qMeasured);
            std::printf("Preflight [%s] from the measured start: %s\n", spec.name.c_str(), from.summary().c_str());
            if (!from.ok) {
                std::printf("Aborted: from where the arm stopped, the motion would violate a margin. No torque "
                            "control (the arm stays at the start).\n");
                return kPreflight;
            }
            if (writeLogCsv(predictionPath, again, again.size()))
                std::printf("  prediction (from the measured start): %s\n", predictionPath.c_str());
        }

        // 8. The torque loop, after Enter, from a state read again after Enter
        std::snprintf(what, sizeof(what),
                      "Path following (TORQUE control) will start: '%s', %.0f s. Keep the user stop in hand.",
                      spec.name.c_str(), spec.duration);
        if (!confirm(opt, what))
            return kAborted;
        state = robot.readOnce();
        if (!atStart(state, setup, false))
            return kNotReady;

        // The log: one row per callback, allocated here (following until the timeout, then at most stopTime of the
        // safety stop, and a margin).
        const std::size_t capacity =
            static_cast<std::size_t>(std::ceil((setup.timeout() + RunConstants::stopTime + 1.0) / RunConstants::dt));
        std::vector<ExperimentLogRow> log(capacity);
        ExperimentLoop loop(setup, setup.monitor());
        FollowResult res;
        int code = kOk;
        try {
            follow(robot, model, loop, opt, log, res);
            code = res.violation != Violation::None ? kSafetyStop : kOk;
        } catch (const franka::ControlException& e) {
            std::printf("ROBOT REFLEX / CONTROL ERROR: %s\n", e.what());
            const std::string reflexPath = logPath + ".reflex.csv";
            std::ofstream f(reflexPath);
            if (f)
                f << franka::logToCSV(e.log);
            std::printf("  libfranka's last %zu states and commands: %s\n", e.log.size(), reflexPath.c_str());
            std::printf("  The robot stopped itself. Recover with Desk (unlock / acknowledge) or '<ip> error-recovery' "
                        "(automaticErrorRecovery), and keep the log for the analysis.\n");
            code = kRobotError;
        } catch (const std::exception& e) {
            // A network error, a realtime error (RealtimeException), or a non-finite command
            std::printf("ERROR in the torque loop: %s\n", e.what());
            std::printf("  The robot stopped. Check Desk, recover if needed, and keep the log.\n");
            code = kRobotError;
        }
        // The log, also after a reflex or an error (the rows before it)
        const std::size_t rows = res.rows;
        const bool written = writeLogCsv(logPath, log, rows);
        std::printf("Log: %s (%zu rows)%s\n", logPath.c_str(), rows, written ? "" : "  COULD NOT BE WRITTEN");
        printRunSummary(log, rows, res, setup, code == kRobotError);
        if (!written && code == kOk)
            code = kLogError;
        return code;
    }

    bool parse(int argc, char** argv, Options& opt)
    {
        if (argc < 3)
            return false;
        opt.ip = argv[1];
        opt.command = argv[2];
        for (int i = 3; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--dry-run")
                opt.dryRun = true;
            else if (a == "--no-rt")
                opt.noRt = true;
            else if (a == "--overwrite")
                opt.overwrite = true;
            else if (a.rfind("--log=", 0) == 0 && a.size() > 6)
                opt.logPath = a.substr(6);
            else if (a.rfind("--log-dir=", 0) == 0 && a.size() > 10)
                opt.logDir = a.substr(10);
#ifdef SE3PF_MOCK_FRANKA
            else if (a == "--yes")
                opt.yes = true;
            else if (a.rfind("--inject=", 0) == 0) {
                const std::string v = a.substr(9);
                const std::size_t at = v.find('@');
                if (at == std::string::npos)
                    return false;
                opt.inject = violationFromName(v.substr(0, at));
                char* end = nullptr;
                opt.injectTime = std::strtod(v.c_str() + at + 1, &end);
                if (opt.inject == Violation::None || end == v.c_str() + at + 1 || *end != '\0'
                    || !(opt.injectTime >= 0.0))
                    return false;
            }
#endif
            else if (opt.command == "preflight" && a.rfind("--", 0) != 0 && opt.which == "all")
                opt.which = a;
            else
                return false;
        }
        return true;
    }

}  // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "list")
        return cmdList();
    Options opt;
    if (!parse(argc, argv, opt)) {
        usage(argv[0]);
        return kUsage;
    }
    const ExperimentSpec* spec = nullptr;
    const bool known = opt.command == "check" || opt.command == "preflight" || opt.command == "error-recovery";
    if (!known) {
        spec = findExperiment(opt.command);
        if (spec == nullptr) {
            std::fprintf(stderr, "Unknown experiment or command '%s' (run '%s list').\n", opt.command.c_str(),
                         argv[0]);
            return kUnknownExperiment;
        }
    }
    if (opt.command == "preflight" && opt.which != "all" && findExperiment(opt.which) == nullptr) {
        std::fprintf(stderr, "Unknown experiment '%s' (run '%s list').\n", opt.which.c_str(), argv[0]);
        return kUnknownExperiment;
    }
#ifndef SE3PF_MOCK_FRANKA
    // libfranka checks realtime only when a motion starts, so check, preflight, error-recovery and --dry-run never
    // need --no-rt; its only effect on a run would be a robot moving without realtime scheduling.
    if (opt.noRt && spec != nullptr && !opt.dryRun) {
        std::fprintf(stderr, "--no-rt is refused for a run that moves the robot: the torque loop needs a PREEMPT_RT "
                             "kernel (docs/HARDWARE.md, section 2).\n");
        return kUsage;
    }
#endif
    if (opt.noRt)
        std::printf("WARNING: --no-rt: no realtime scheduling (no motion is allowed with it on the robot).\n");
    std::fflush(stdout);  // stdout before anything goes to stderr

    try {
        const PandaExample example;  // the paper's Panda example: path, tube, lift, kernel frame (offline, < 1 s)
        franka::Robot robot(opt.ip, opt.noRt ? franka::RealtimeConfig::kIgnore : franka::RealtimeConfig::kEnforce);
        if (opt.command == "check")
            return cmdCheck(robot, example);
        if (opt.command == "preflight")
            return cmdPreflight(robot, example, opt);
        if (opt.command == "error-recovery")
            return cmdRecover(robot, opt);
        return cmdRun(robot, example, *spec, opt);
    } catch (const franka::IncompatibleVersionException& e) {
        std::fflush(stdout);
        std::fprintf(stderr, "%s\n", e.what());
        std::fprintf(stderr, "The robot's server version is %u; this libfranka speaks %u. Server version 5: "
                             "libfranka 0.9.2 (system >= 4.2.1); 4: libfranka 0.8.0 (system 4.0 to 4.2.0); 3: "
                             "libfranka 0.7.1 (system 3.x). See docs/HARDWARE.md, section 2.\n",
                     static_cast<unsigned>(e.server_version), static_cast<unsigned>(e.library_version));
        return kRobotError;
    } catch (const franka::CommandException& e) {
        std::fflush(stdout);
        std::fprintf(stderr, "%s\n", e.what());
        std::fprintf(stderr, "The robot rejected a command: check its mode in Desk (joints unlocked, user stop "
                             "released, FCI active, no reflex: '%s <ip> error-recovery').\n",
                     argv[0]);
        return kRobotError;
    } catch (const franka::RealtimeException& e) {
        std::fflush(stdout);
        std::fprintf(stderr, "%s\n", e.what());
        std::fprintf(stderr, "No motion started: the realtime setup of this PC (PREEMPT_RT kernel, ulimit -r 99: "
                             "docs/HARDWARE.md, section 2).\n");
        return kRobotError;
    } catch (const franka::Exception& e) {
        std::fflush(stdout);
        std::fprintf(stderr, "%s\n", e.what());  // libfranka's messages start with "libfranka: "
        std::fprintf(stderr, "Check the network (ping the robot; the firewall must let the robot's UDP packets in), "
                             "that FCI is active in Desk, and the robot system version (libfranka 0.9.2 needs Panda "
                             "system 4.2.1 or newer). After a reflex: Desk, or '%s <ip> error-recovery'.\n",
                     argv[0]);
        return kRobotError;
    } catch (const std::exception& e) {
        std::fflush(stdout);
        std::fprintf(stderr, "error: %s\n", e.what());
        return kRobotError;
    }
}
