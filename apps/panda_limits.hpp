// The Panda's joint limits as a CSV file, for the plot scripts (tools/plot_panda.py, tools/plot_mujoco.py), so that
// they need no copy of the limits. Written by panda_sim and panda_mujoco. Header-only, for the apps.
#pragma once
#include "se3pf/panda.hpp"

#include <filesystem>
#include <fstream>

namespace se3pf {

    // One row per joint, numbered from 1: position, velocity, acceleration and torque limits of PandaArm
    // (src/panda.cpp). Columns: joint, q_min, q_max, qd_max, qdd_max, tau_max. False if the file cannot be written.
    inline bool writePandaLimits(const std::filesystem::path& file)
    {
        std::ofstream f(file);
        if (!f)
            return false;
        const Vec7 qMin = PandaArm::qMin(), qMax = PandaArm::qMax();
        const Vec7 dqMax = PandaArm::dqMax(), ddqMax = PandaArm::ddqMax(), tauMax = PandaArm::tauMax();
        f << "joint,q_min,q_max,qd_max,qdd_max,tau_max\n";
        f.precision(10);  // as the logs; the limits have at most 5 significant digits, so they are read back exactly
        for (int i = 0; i < 7; ++i)
            f << i + 1 << ',' << qMin[i] << ',' << qMax[i] << ',' << dqMax[i] << ',' << ddqMax[i] << ',' << tauMax[i]
              << '\n';
        return static_cast<bool>(f);
    }

}  // namespace se3pf
