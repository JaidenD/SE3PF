// Mock of <franka/log.h> (libfranka 0.9.2): the records a ControlException carries. See mock/franka/README.md.
#pragma once
#include <sstream>
#include <string>
#include <vector>

#include "franka/control_types.h"
#include "franka/robot_state.h"

namespace franka {

    struct RobotCommand {
        JointPositions joint_positions{0, 0, 0, 0, 0, 0, 0};
        JointVelocities joint_velocities{0, 0, 0, 0, 0, 0, 0};
        Torques torques{0, 0, 0, 0, 0, 0, 0};
    };

    struct Record {
        RobotState state;
        RobotCommand command;
    };

    // A CSV of the records (the real one has many more columns): time, q, dq, tau_J, commanded torque
    inline std::string logToCSV(const std::vector<Record>& log)
    {
        std::ostringstream os;
        os << "time";
        for (const char* name : {"q", "dq", "tau_J", "tau_J_d_cmd"})
            for (int i = 0; i < 7; ++i)
                os << ',' << name << i;
        os << '\n';
        for (const Record& r : log) {
            os << r.state.time.toMSec();
            for (double x : r.state.q)
                os << ',' << x;
            for (double x : r.state.dq)
                os << ',' << x;
            for (double x : r.state.tau_J)
                os << ',' << x;
            for (double x : r.command.torques.tau_J)
                os << ',' << x;
            os << '\n';
        }
        return os.str();
    }

}  // namespace franka
