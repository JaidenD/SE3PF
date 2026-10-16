// Mock of <franka/robot_state.h> (libfranka 0.9.2): every field of the real RobotState except the two Errors
// (current_errors, last_motion_errors), which the driver does not use. See mock/franka/README.md.
#pragma once
#include <array>
#include <ostream>

#include "franka/duration.h"

namespace franka {

    enum class RobotMode { kOther, kIdle, kMove, kGuiding, kReflex, kUserStopped, kAutomaticErrorRecovery };

    struct RobotState {
        std::array<double, 16> O_T_EE{};
        std::array<double, 16> O_T_EE_d{};
        std::array<double, 16> F_T_EE{};
        std::array<double, 16> F_T_NE{};
        std::array<double, 16> NE_T_EE{};
        std::array<double, 16> EE_T_K{};
        double m_ee{};
        std::array<double, 9> I_ee{};
        std::array<double, 3> F_x_Cee{};
        double m_load{};
        std::array<double, 9> I_load{};
        std::array<double, 3> F_x_Cload{};
        double m_total{};
        std::array<double, 9> I_total{};
        std::array<double, 3> F_x_Ctotal{};
        std::array<double, 2> elbow{};
        std::array<double, 2> elbow_d{};
        std::array<double, 2> elbow_c{};
        std::array<double, 2> delbow_c{};
        std::array<double, 2> ddelbow_c{};
        std::array<double, 7> tau_J{};
        std::array<double, 7> tau_J_d{};
        std::array<double, 7> dtau_J{};
        std::array<double, 7> q{};
        std::array<double, 7> q_d{};
        std::array<double, 7> dq{};
        std::array<double, 7> dq_d{};
        std::array<double, 7> ddq_d{};
        std::array<double, 7> joint_contact{};
        std::array<double, 6> cartesian_contact{};
        std::array<double, 7> joint_collision{};
        std::array<double, 6> cartesian_collision{};
        std::array<double, 7> tau_ext_hat_filtered{};
        std::array<double, 6> O_F_ext_hat_K{};
        std::array<double, 6> K_F_ext_hat_K{};
        std::array<double, 6> O_dP_EE_d{};
        std::array<double, 3> O_ddP_O{};
        std::array<double, 16> O_T_EE_c{};
        std::array<double, 6> O_dP_EE_c{};
        std::array<double, 6> O_ddP_EE_c{};
        std::array<double, 7> theta{};
        std::array<double, 7> dtheta{};
        double control_command_success_rate{};
        RobotMode robot_mode = RobotMode::kUserStopped;
        Duration time{};
    };

}  // namespace franka
