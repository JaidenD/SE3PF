// Mock of <franka/robot.h> (libfranka 0.9.2): a simulated Panda behind the same interface (the overloads of control
// the driver uses, readOnce, loadModel, the behavior setters, automaticErrorRecovery, stop, serverVersion). The
// control loops run at 1 kHz of simulated time as fast as the computer allows. See mock/franka/README.md for what
// the simulated robot does and the environment variables that configure it.
#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "franka/control_types.h"
#include "franka/duration.h"
#include "franka/exception.h"
#include "franka/log.h"
#include "franka/lowpass_filter.h"
#include "franka/mock_arm.h"
#include "franka/model.h"
#include "franka/rate_limiting.h"
#include "franka/robot_state.h"

namespace franka {

    class Robot {
        public:
            using ServerVersion = uint16_t;

            explicit Robot(const std::string& franka_address, RealtimeConfig realtime_config = RealtimeConfig::kEnforce,
                           size_t log_size = 50)
                : config_(mock::Config::fromEnvironment()), logSize_(log_size), realtime_(realtime_config)
            {
                (void)franka_address;
                configuredLoad_ = config_.configuredLoad ? config_.load : mock::frankaHand();
                plantLoad_ = configuredLoad_;
                if (config_.plantLoad > 0.0) {
                    // an extra mass at the flange that Desk does not know: combined centre of mass
                    const double m = plantLoad_.mass + config_.plantLoad;
                    plantLoad_.com = (plantLoad_.mass * plantLoad_.com + config_.plantLoad * Eigen::Vector3d(0, 0, 0.05)) / m;
                    plantLoad_.mass = m;
                }
                // The robot holds q_d; the measured q is off by a small deflection (the joints' elasticity)
                qD_ = config_.q0;
                q_ = qD_ + mock::Vec7::Constant(config_.holdError);
                dq_ = mock::Vec7::Zero();
                tauCmd_ = mock::Vec7::Zero();
                mode_ = config_.mode == "reflex"         ? RobotMode::kReflex
                        : config_.mode == "user-stopped" ? RobotMode::kUserStopped
                                                         : RobotMode::kIdle;
                upperTorqueThresholds_ = {{20.0, 20.0, 18.0, 18.0, 16.0, 14.0, 12.0}};
                updateState();
            }
            Robot(Robot&&) noexcept = default;
            Robot& operator=(Robot&&) noexcept = default;
            ~Robot() noexcept = default;
            Robot(const Robot&) = delete;
            Robot& operator=(const Robot&) = delete;

            // Torque control: the callback's torque is filtered and rate limited as libfranka does, the robot adds
            // G of its model, and the plant integrates for 1 ms. With FRANKA_MOCK_LATE some periods last longer: the
            // robot keeps the last command for the whole period and the next callback gets the longer period.
            void control(std::function<Torques(const RobotState&, franka::Duration)> control_callback,
                         bool limit_rate = true, double cutoff_frequency = kDefaultCutoffFrequency)
            {
                if (!control_callback)
                    throw std::invalid_argument("libfranka: Invalid control callback given.");
                requireRealtime();
                startMotion();
                const double dt = kDeltaT;
                const long maxSteps = static_cast<long>(config_.maxSeconds / dt);
                const mock::ToyArm plant(1.0), model(1.0);
                const Eigen::Vector3d g(0.0, 0.0, -9.81);
                Duration period(0);
                double tLoop = 0.0;
                bool nanDone = false;
                for (long k = 0; k < maxSteps; ++k) {
                    // A sensor glitch: one period with a NaN joint velocity
                    RobotState s = state_;
                    s.q_d = s.q;  // no motion generator in torque control
                    if (config_.nanJoint >= 0 && !nanDone && tLoop >= config_.nanTime) {
                        s.dq[static_cast<std::size_t>(config_.nanJoint)] = std::numeric_limits<double>::quiet_NaN();
                        nanDone = true;
                    }
                    const Torques cmd = control_callback(s, period);
                    std::array<double, 7> tau = cmd.tau_J;
                    if (cutoff_frequency < kMaxCutoffFrequency)
                        for (std::size_t i = 0; i < 7; ++i)
                            tau[i] = lowpassFilter(kDeltaT, tau[i], state_.tau_J_d[i], cutoff_frequency);
                    if (limit_rate)
                        tau = limitRate(kMaxTorqueRate, tau, state_.tau_J_d);
                    for (double x : tau)
                        if (!std::isfinite(x))
                            throw std::invalid_argument("Commanding value is infinite or NaN.");
                    tauCmd_ = mock::fromArray(tau);

                    const int periodMs =
                        config_.lateEvery > 0 && (k + 1) % config_.lateEvery == 0 ? config_.lateMs : 1;
                    for (int m = 0; m < periodMs; ++m) {
                        // The robot adds G of its model with the configured load; the plant may carry another load
                        const mock::Vec7 motor = tauCmd_ + model.gravity(q_, configuredLoad_, g);
                        const mock::Vec7 tauMax = se3pf::PandaArm::tauMax();
                        for (int i = 0; i < 7; ++i)
                            if (!(std::abs(motor[i]) <= tauMax[i]))
                                reflex("tau_J_range_violation (mock: joint torque limit)", cmd);
                        const mock::Vec7 ext = push(tLoop);
                        for (int i = 0; i < 7; ++i)
                            if (std::abs(ext[i]) > upperTorqueThresholds_[static_cast<std::size_t>(i)])
                                reflex("joint_reflex (mock: external torque above the collision threshold)", cmd);
                        tauExt_ = ext;

                        for (int sub = 0; sub < 2; ++sub) {  // semi-implicit Euler, 0.5 ms
                            const double h = 0.5 * dt;
                            const mock::Vec7 rhs = motor + ext - config_.friction * dq_
                                                   - plant.coriolis(q_, dq_, plantLoad_) - plant.gravity(q_, plantLoad_, g);
                            const mock::Vec7 ddq = plant.mass(q_, plantLoad_).ldlt().solve(rhs);
                            dq_ += h * ddq;
                            q_ += h * dq_;
                        }
                        const mock::Vec7 qMin = se3pf::PandaArm::qMin(), qMax = se3pf::PandaArm::qMax();
                        const mock::Vec7 dqMax = se3pf::PandaArm::dqMax();
                        for (int i = 0; i < 7; ++i) {
                            if (!(q_[i] >= qMin[i] && q_[i] <= qMax[i]))
                                reflex("joint_position_limits_violation (mock)", cmd);
                            if (!(std::abs(dq_[i]) <= dqMax[i]))
                                reflex("joint_velocity_violation (mock)", cmd);
                        }
                        tLoop += dt;
                    }
                    record(cmd);
                    ms_ += static_cast<uint64_t>(periodMs);
                    period = Duration(static_cast<uint64_t>(periodMs));
                    successRate_ = config_.lateEvery > 0 ? 1.0 - (config_.lateMs - 1.0) / config_.lateEvery : 1.0;
                    updateState();
                    if (cmd.motion_finished) {
                        finishMotion();
                        return;
                    }
                }
                reflex("mock: the control loop did not finish within FRANKA_MOCK_MAX_SECONDS", Torques(tau0()));
            }

            // Joint position motion: the arm follows the command (the measured q off by the hold deflection). As on
            // the robot, the first command must be q_d (the last commanded or held position, not the measured q),
            // and the commanded motion must stay within the Panda's velocity, acceleration and jerk limits and end at
            // rest. With FRANKA_MOCK_APPROACH_ERROR the arm settles short of the goal.
            void control(std::function<JointPositions(const RobotState&, franka::Duration)> motion_generator_callback,
                         ControllerMode controller_mode = ControllerMode::kJointImpedance, bool limit_rate = true,
                         double cutoff_frequency = kDefaultCutoffFrequency)
            {
                (void)controller_mode;
                (void)limit_rate;
                if (!motion_generator_callback)
                    throw std::invalid_argument("libfranka: Invalid motion callback given.");
                requireRealtime();
                startMotion();
                const double dt = kDeltaT;
                const long maxSteps = static_cast<long>(config_.maxSeconds / dt);
                Duration period(0);
                const mock::Vec7 qFrom = qD_;
                mock::Vec7 vPrev = mock::Vec7::Zero(), aPrev = mock::Vec7::Zero();
                const mock::Vec7 qMin = se3pf::PandaArm::qMin(), qMax = se3pf::PandaArm::qMax();
                const mock::Vec7 dqMax = se3pf::PandaArm::dqMax(), ddqMax = se3pf::PandaArm::ddqMax();
                const mock::Vec7 jerkMax = se3pf::PandaArm::jerkMax();
                const mock::Vec7 deflection = mock::Vec7::Constant(config_.holdError);
                for (long k = 0; k < maxSteps; ++k) {
                    const JointPositions cmd = motion_generator_callback(state_, period);
                    std::array<double, 7> qc = cmd.q;
                    for (double x : qc)
                        if (!std::isfinite(x))
                            throw std::invalid_argument("Commanding value is infinite or NaN.");
                    if (cutoff_frequency < kMaxCutoffFrequency)
                        for (std::size_t i = 0; i < 7; ++i)
                            qc[i] = lowpassFilter(kDeltaT, qc[i], state_.q_d[i], cutoff_frequency);
                    const mock::Vec7 qn = mock::fromArray(qc);
                    if (k == 0 && (mock::fromArray(cmd.q) - mock::fromArray(state_.q_d)).cwiseAbs().maxCoeff() > 1e-9)
                        reflex("joint_motion_generator_start_pose_invalid (mock: the first command must be q_d)",
                               Torques(tau0()));
                    // The commanded signal's derivatives, as the robot checks them
                    const mock::Vec7 v = (qn - qD_) / dt, a = (v - vPrev) / dt, j = (a - aPrev) / dt;
                    for (int i = 0; i < 7; ++i) {
                        if (!(qn[i] >= qMin[i] && qn[i] <= qMax[i]))
                            reflex("joint_motion_generator_position_limits_violation (mock)", Torques(tau0()));
                        if (std::abs(v[i]) > dqMax[i])
                            reflex("joint_motion_generator_velocity_limits_violation (mock)", Torques(tau0()));
                        if (std::abs(a[i]) > ddqMax[i])
                            reflex("joint_motion_generator_acceleration_discontinuity (mock)", Torques(tau0()));
                        if (std::abs(j[i]) > jerkMax[i])
                            reflex("joint_motion_generator_jerk_discontinuity (mock)", Torques(tau0()));
                    }
                    vPrev = v;
                    aPrev = a;
                    qD_ = qn;
                    q_ = qD_ + deflection;
                    dq_ = v;
                    ++ms_;
                    period = Duration(1);
                    updateState();
                    if (cmd.motion_finished) {
                        if (v.cwiseAbs().maxCoeff() > 1e-2)
                            reflex("mock: the joint motion finished while moving", Torques(tau0()));
                        // The arm settles approachError short of the goal (q_d stays the last command)
                        for (int i = 0; i < 7; ++i) {
                            const double dir = qD_[i] > qFrom[i] ? 1.0 : (qD_[i] < qFrom[i] ? -1.0 : 0.0);
                            q_[i] -= config_.approachError * dir;
                        }
                        finishMotion(false);
                        return;
                    }
                }
                reflex("mock: the motion did not finish within FRANKA_MOCK_MAX_SECONDS", Torques(tau0()));
            }

            // FRANKA_MOCK_USER_STOP_AFTER_READS=n: someone presses the user stop right after the n-th call (while the
            // operator waits at a prompt, say); the robot is UserStopped from then on
            RobotState readOnce()
            {
                const RobotState s = state_;
                if (config_.userStopAfterReads > 0 && ++reads_ == config_.userStopAfterReads) {
                    mode_ = RobotMode::kUserStopped;
                    updateState();
                }
                return s;
            }

            // The setters, as libfranka: rejected (CommandException) unless the robot is Idle
            void setCollisionBehavior(const std::array<double, 7>& lower_torque_thresholds_acceleration,
                                      const std::array<double, 7>& upper_torque_thresholds_acceleration,
                                      const std::array<double, 7>& lower_torque_thresholds_nominal,
                                      const std::array<double, 7>& upper_torque_thresholds_nominal,
                                      const std::array<double, 6>& lower_force_thresholds_acceleration,
                                      const std::array<double, 6>& upper_force_thresholds_acceleration,
                                      const std::array<double, 6>& lower_force_thresholds_nominal,
                                      const std::array<double, 6>& upper_force_thresholds_nominal)
            {
                (void)lower_torque_thresholds_acceleration;
                (void)upper_torque_thresholds_acceleration;
                (void)lower_torque_thresholds_nominal;
                (void)lower_force_thresholds_acceleration;
                (void)upper_force_thresholds_acceleration;
                (void)lower_force_thresholds_nominal;
                (void)upper_force_thresholds_nominal;
                requireIdle("Set Collision Behavior");
                upperTorqueThresholds_ = upper_torque_thresholds_nominal;
            }
            void setCollisionBehavior(const std::array<double, 7>& lower_torque_thresholds,
                                      const std::array<double, 7>& upper_torque_thresholds,
                                      const std::array<double, 6>& lower_force_thresholds,
                                      const std::array<double, 6>& upper_force_thresholds)
            {
                (void)lower_torque_thresholds;
                (void)lower_force_thresholds;
                (void)upper_force_thresholds;
                requireIdle("Set Collision Behavior");
                upperTorqueThresholds_ = upper_torque_thresholds;
            }
            void setJointImpedance(const std::array<double, 7>& K_theta)  // NOLINT
            {
                (void)K_theta;
                requireIdle("Set Joint Impedance");
            }
            void setCartesianImpedance(const std::array<double, 6>& K_x)  // NOLINT
            {
                (void)K_x;
                requireIdle("Set Cartesian Impedance");
            }

            void automaticErrorRecovery()
            {
                mode_ = RobotMode::kIdle;
                updateState();
            }
            void stop() {}
            Model loadModel() { return Model(mock::ToyArm(1.0), config_.poseError); }
            ServerVersion serverVersion() const noexcept { return 5; }  // the robot server version of libfranka 0.9.x

        private:
            static std::array<double, 7> tau0() { return {}; }

            static const char* modeString(RobotMode m)
            {
                switch (m) {
                    case RobotMode::kUserStopped: return "\"User stopped\"";
                    case RobotMode::kIdle: return "\"Idle\"";
                    case RobotMode::kMove: return "\"Move\"";
                    case RobotMode::kGuiding: return "\"Guiding\"";
                    case RobotMode::kReflex: return "\"Reflex\"";
                    case RobotMode::kAutomaticErrorRecovery: return "\"Automatic error recovery\"";
                    case RobotMode::kOther: return "\"Other\"";
                }
                return "\"Other\"";
            }
            // libfranka's message for a command the robot rejects in its current mode (robot_impl.h)
            void requireIdle(const char* command) const
            {
                if (mode_ != RobotMode::kIdle)
                    throw CommandException(std::string("libfranka: ") + command
                                           + " command rejected: command not possible in the current mode ("
                                           + modeString(mode_) + ")!"
                                           + (mode_ == RobotMode::kOther ? " Did you open the brakes?" : ""));
            }
            // libfranka 0.9.2 checks the kernel when a motion starts (ControlLoop), not when it connects
            void requireRealtime() const
            {
                if (realtime_ == RealtimeConfig::kEnforce && config_.noRealtime)
                    throw RealtimeException("libfranka: Running kernel does not have realtime capabilities.");
            }

            // As the robot: no motion unless it is idle (user stop released, no reflex)
            void startMotion()
            {
                requireIdle("Move");
                mode_ = RobotMode::kMove;
                records_.clear();
                updateState();
            }
            // afterTorque: the robot holds the arm where the torque loop left it
            void finishMotion(bool afterTorque = true)
            {
                dq_ = mock::Vec7::Zero();  // the robot holds the arm after a motion
                if (afterTorque)
                    qD_ = q_;
                tauCmd_ = mock::Vec7::Zero();
                tauExt_ = mock::Vec7::Zero();
                successRate_ = 1.0;
                mode_ = RobotMode::kIdle;
                updateState();
            }
            [[noreturn]] void reflex(const std::string& what, const Torques& cmd)
            {
                record(cmd);
                mode_ = RobotMode::kReflex;
                dq_ = mock::Vec7::Zero();  // the brakes stop the arm
                qD_ = q_;
                updateState();
                throw ControlException("libfranka: Move command aborted: motion aborted by reflex! [" + what + "]",
                                       std::vector<Record>(records_.begin(), records_.end()));
            }
            void record(const Torques& cmd)
            {
                Record r;
                r.state = state_;
                r.command.torques = cmd;
                records_.push_back(r);
                while (records_.size() > logSize_)
                    records_.pop_front();
            }
            mock::Vec7 push(double t) const
            {
                mock::Vec7 ext = mock::Vec7::Zero();
                if (config_.pushJoint >= 0 && t >= config_.pushStart && t < config_.pushStart + config_.pushDuration)
                    ext[config_.pushJoint] = config_.pushTorque;
                return ext;
            }
            void updateState()
            {
                const mock::ToyArm model(1.0);
                state_.q = mock::toArray(q_);
                state_.q_d = mock::toArray(qD_);
                state_.dq = mock::toArray(dq_);
                state_.dq_d = state_.dq;
                state_.tau_J_d = mock::toArray(tauCmd_);
                state_.tau_J = mock::toArray(
                    tauCmd_ + model.gravity(q_, configuredLoad_, Eigen::Vector3d(0, 0, -9.81)) + tauExt_);
                state_.tau_ext_hat_filtered = mock::toArray(tauExt_);
                const se3pf::Pose T = se3pf::PandaArm::output(q_, nullptr);
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c)
                        state_.O_T_EE[static_cast<std::size_t>(c * 4 + r)] = T(r, c);
                state_.F_T_EE = {{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}};
                state_.EE_T_K = state_.F_T_EE;
                state_.m_ee = configuredLoad_.mass;
                state_.F_x_Cee = {{configuredLoad_.com[0], configuredLoad_.com[1], configuredLoad_.com[2]}};
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c)
                        state_.I_ee[static_cast<std::size_t>(c * 3 + r)] = configuredLoad_.inertia(r, c);
                state_.m_load = 0.0;
                state_.m_total = state_.m_ee;
                state_.F_x_Ctotal = state_.F_x_Cee;
                state_.I_total = state_.I_ee;
                state_.control_command_success_rate = successRate_;
                state_.robot_mode = mode_;
                state_.time = Duration(ms_);
            }

            mock::Config config_;
            std::size_t logSize_;
            RealtimeConfig realtime_;
            mock::RigidBody configuredLoad_, plantLoad_;
            mock::Vec7 q_, qD_, dq_, tauCmd_, tauExt_ = mock::Vec7::Zero();
            RobotMode mode_ = RobotMode::kIdle;
            std::array<double, 7> upperTorqueThresholds_{};
            uint64_t ms_ = 0;
            int reads_ = 0;
            double successRate_ = 1.0;
            std::deque<Record> records_;
            RobotState state_;
    };

}  // namespace franka
