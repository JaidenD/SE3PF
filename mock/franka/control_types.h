// Mock of <franka/control_types.h> (libfranka 0.9.2): the subset the driver uses. See mock/franka/README.md.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace franka {

    enum class ControllerMode { kJointImpedance, kCartesianImpedance };
    enum class RealtimeConfig { kEnforce, kIgnore };

    struct Finishable {
        bool motion_finished = false;
    };

    class Torques : public Finishable {
        public:
            Torques(const std::array<double, 7>& torques) noexcept : tau_J(torques) {}
            Torques(std::initializer_list<double> torques)
            {
                if (torques.size() != 7)
                    throw std::invalid_argument("Invalid size of torques");
                std::copy(torques.begin(), torques.end(), tau_J.begin());
            }
            std::array<double, 7> tau_J{};  // NOLINT(readability-identifier-naming)
    };

    class JointPositions : public Finishable {
        public:
            JointPositions(const std::array<double, 7>& joint_positions) noexcept : q(joint_positions) {}
            JointPositions(std::initializer_list<double> joint_positions)
            {
                if (joint_positions.size() != 7)
                    throw std::invalid_argument("Invalid size of joint positions");
                std::copy(joint_positions.begin(), joint_positions.end(), q.begin());
            }
            std::array<double, 7> q{};
    };

    class JointVelocities : public Finishable {
        public:
            JointVelocities(const std::array<double, 7>& joint_velocities) noexcept : dq(joint_velocities) {}
            JointVelocities(std::initializer_list<double> joint_velocities)
            {
                if (joint_velocities.size() != 7)
                    throw std::invalid_argument("Invalid size of joint velocities");
                std::copy(joint_velocities.begin(), joint_velocities.end(), dq.begin());
            }
            std::array<double, 7> dq{};
    };

    inline Torques MotionFinished(Torques command) noexcept  // NOLINT(readability-identifier-naming)
    {
        command.motion_finished = true;
        return command;
    }

    inline JointPositions MotionFinished(JointPositions command) noexcept  // NOLINT(readability-identifier-naming)
    {
        command.motion_finished = true;
        return command;
    }

    inline JointVelocities MotionFinished(JointVelocities command) noexcept  // NOLINT(readability-identifier-naming)
    {
        command.motion_finished = true;
        return command;
    }

}  // namespace franka
