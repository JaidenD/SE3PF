// Mock of <franka/model.h> (libfranka 0.9.2): the same functions and overloads. Kinematics of the flange are
// se3pf::PandaArm's (so the driver's kinematics check passes); the dynamics are the toy arm's (mock_arm.h), with the
// end-effector load passed in (the RobotState overloads take the state's m_total, F_x_Ctotal and I_total, as the real
// ones do). Only Frame::kFlange is supported. See mock/franka/README.md.
#pragma once
#include <array>

#include "franka/exception.h"
#include "franka/mock_arm.h"
#include "franka/robot_state.h"

namespace franka {

    enum class Frame { kJoint1, kJoint2, kJoint3, kJoint4, kJoint5, kJoint6, kJoint7, kFlange, kEndEffector, kStiffness };

    class Model {
        public:
            // Mock only (the real one is built by Robot::loadModel from the robot's network connection)
            Model(const mock::ToyArm& arm, double poseError) : arm_(arm), poseError_(poseError) {}
            Model(Model&& model) noexcept = default;
            Model& operator=(Model&& model) noexcept = default;
            ~Model() noexcept = default;
            Model(const Model&) = delete;
            Model& operator=(const Model&) = delete;

            std::array<double, 16> pose(Frame frame, const franka::RobotState& robot_state) const
            {
                return pose(frame, robot_state.q, robot_state.F_T_EE, robot_state.EE_T_K);
            }
            std::array<double, 16> pose(Frame frame, const std::array<double, 7>& q,
                                        const std::array<double, 16>& /*F_T_EE*/,
                                        const std::array<double, 16>& /*EE_T_K*/) const
            {
                requireFlange(frame);
                se3pf::Pose T = se3pf::PandaArm::output(mock::fromArray(q), nullptr);
                T(0, 3) += poseError_;
                std::array<double, 16> a{};
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c)
                        a[static_cast<std::size_t>(c * 4 + r)] = T(r, c);  // column-major
                return a;
            }

            std::array<double, 42> bodyJacobian(Frame frame, const franka::RobotState& robot_state) const
            {
                return bodyJacobian(frame, robot_state.q, robot_state.F_T_EE, robot_state.EE_T_K);
            }
            std::array<double, 42> bodyJacobian(Frame frame, const std::array<double, 7>& q,
                                                const std::array<double, 16>& /*F_T_EE*/,
                                                const std::array<double, 16>& /*EE_T_K*/) const
            {
                requireFlange(frame);
                se3pf::Mat67 J = se3pf::Mat67::Zero();
                se3pf::PandaArm::output(mock::fromArray(q), &J);
                return franka6x7(J, true);
            }

            std::array<double, 42> zeroJacobian(Frame frame, const franka::RobotState& robot_state) const
            {
                return zeroJacobian(frame, robot_state.q, robot_state.F_T_EE, robot_state.EE_T_K);
            }
            std::array<double, 42> zeroJacobian(Frame frame, const std::array<double, 7>& q,
                                                const std::array<double, 16>& /*F_T_EE*/,
                                                const std::array<double, 16>& /*EE_T_K*/) const
            {
                requireFlange(frame);
                se3pf::Mat67 J = se3pf::Mat67::Zero();
                const se3pf::Pose T = se3pf::PandaArm::output(mock::fromArray(q), &J);
                // In the base frame: rotate both blocks of the body Jacobian by R
                se3pf::Mat67 J0 = se3pf::Mat67::Zero();
                J0.topRows<3>() = T.block<3, 3>(0, 0) * J.topRows<3>();
                J0.bottomRows<3>() = T.block<3, 3>(0, 0) * J.bottomRows<3>();
                return franka6x7(J0, true);
            }

            std::array<double, 49> mass(const franka::RobotState& robot_state) const noexcept
            {
                return mass(robot_state.q, robot_state.I_total, robot_state.m_total, robot_state.F_x_Ctotal);
            }
            std::array<double, 49> mass(const std::array<double, 7>& q, const std::array<double, 9>& I_total,
                                        double m_total, const std::array<double, 3>& F_x_Ctotal) const noexcept
            {
                const mock::Mat7 M = arm_.mass(mock::fromArray(q), mock::bodyFromArrays(m_total, F_x_Ctotal, I_total));
                std::array<double, 49> a{};
                for (int r = 0; r < 7; ++r)
                    for (int c = 0; c < 7; ++c)
                        a[static_cast<std::size_t>(c * 7 + r)] = M(r, c);  // column-major
                return a;
            }

            std::array<double, 7> coriolis(const franka::RobotState& robot_state) const noexcept
            {
                return coriolis(robot_state.q, robot_state.dq, robot_state.I_total, robot_state.m_total,
                                robot_state.F_x_Ctotal);
            }
            std::array<double, 7> coriolis(const std::array<double, 7>& q, const std::array<double, 7>& dq,
                                           const std::array<double, 9>& I_total, double m_total,
                                           const std::array<double, 3>& F_x_Ctotal) const noexcept
            {
                return mock::toArray(arm_.coriolis(mock::fromArray(q), mock::fromArray(dq),
                                                   mock::bodyFromArrays(m_total, F_x_Ctotal, I_total)));
            }

            std::array<double, 7> gravity(const std::array<double, 7>& q, double m_total,
                                          const std::array<double, 3>& F_x_Ctotal,
                                          const std::array<double, 3>& gravity_earth = {{0., 0., -9.81}}) const noexcept
            {
                const mock::RigidBody load = mock::bodyFromArrays(m_total, F_x_Ctotal, std::array<double, 9>{});
                return mock::toArray(arm_.gravity(mock::fromArray(q), load,
                                                  Eigen::Vector3d(gravity_earth[0], gravity_earth[1], gravity_earth[2])));
            }
            std::array<double, 7> gravity(const franka::RobotState& robot_state,
                                          const std::array<double, 3>& gravity_earth) const noexcept
            {
                return gravity(robot_state.q, robot_state.m_total, robot_state.F_x_Ctotal, gravity_earth);
            }
            std::array<double, 7> gravity(const franka::RobotState& robot_state) const noexcept
            {
                return gravity(robot_state.q, robot_state.m_total, robot_state.F_x_Ctotal);
            }

        private:
            static void requireFlange(Frame frame)
            {
                if (frame != Frame::kFlange)
                    throw ModelException("mock franka::Model: only Frame::kFlange is supported");
            }
            // libfranka's 6 x 7 Jacobians: column-major, rows (v, omega); ours are rows (omega, v)
            static std::array<double, 42> franka6x7(const se3pf::Mat67& J, bool swapRows)
            {
                std::array<double, 42> a{};
                for (int c = 0; c < 7; ++c)
                    for (int r = 0; r < 6; ++r) {
                        const int src = swapRows ? (r < 3 ? r + 3 : r - 3) : r;
                        a[static_cast<std::size_t>(c * 6 + r)] = J(src, c);
                    }
                return a;
            }

            mock::ToyArm arm_;
            double poseError_ = 0.0;
    };

}  // namespace franka
