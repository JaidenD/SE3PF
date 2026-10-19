#pragma once
// The CSV log of a run on the robot (apps/robot/robot_driver.cpp), and of its predictions in the same format: the
// preflight (robot_driver writes <log>.preflight.csv next to every run's log) and the MuJoCo closed loop
// (apps/experiments_mujoco.cpp). tools/plot_robot_log.py reads them by column name. No libfranka header.
//
// Columns, one row per control cycle (joints numbered from 1 as on the robot):
//   t, phase (0 following, 1 safety stop: damped hold), eta, eta_dot, ref (nu_d or eta_d), e_tan, xi1..xi5, xi_norm,
//   beta, beta_dot, sigma_min, ramp, newton, q1..q7, dq1..dq7, tau1..tau7 (commanded torque, without G),
//   tauG1..tauG7 (tau + G: the torque the joints deliver), qdd1..qdd7 (qdd*), tauF1..tauF7 (the friction
//   feedforward tau_f, part of tau; CTRL-1), tauJ1..tauJ7 (driver only: the measured link-side torque tau_J of the
//   robot state, gravity included, to compare with tauG), int_eta, int_xi, int_beta (the integral terms in mu, |.|
//   of the xi one; CTRL-1), compute_us, period_ms, success_rate (the last three are the driver's: time spent in the
//   callback, libfranka's period, and the robot's control_command_success_rate; 0 in predictions, as tauJ).
#include "se3pf/experiments.hpp"

#include <cstddef>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

namespace se3pf::robot {

    inline void writeLogHeader(std::ostream& os)
    {
        os << "t,phase,eta,eta_dot,ref,e_tan,xi1,xi2,xi3,xi4,xi5,xi_norm,beta,beta_dot,sigma_min,ramp,newton";
        for (const char* name : {"q", "dq", "tau", "tauG", "qdd", "tauF", "tauJ"})
            for (int i = 1; i <= 7; ++i)
                os << ',' << name << i;
        os << ",int_eta,int_xi,int_beta,compute_us,period_ms,success_rate\n";
    }

    inline void writeLogRow(std::ostream& os, const ExperimentLogRow& r)
    {
        os << r.t << ',' << r.phase << ',' << r.eta << ',' << r.etaDot << ',' << r.ref << ',' << r.eTan;
        for (int i = 0; i < 5; ++i)
            os << ',' << r.xi[i];
        os << ',' << r.xiNorm << ',' << r.beta << ',' << r.betaDot << ',' << r.sigmaMin << ',' << r.ramp << ','
           << r.newton;
        for (const Vec7* v : {&r.q, &r.dq, &r.tau, &r.tauG, &r.qdd, &r.tauF, &r.tauJ})
            for (int i = 0; i < 7; ++i)
                os << ',' << (*v)[i];
        os << ',' << r.iEta << ',' << r.iXi << ',' << r.iBeta << ',' << r.computeUs << ',' << r.periodMs << ','
           << r.successRate << '\n';
    }

    // Writes rows [0, count) to path; false if the file cannot be written.
    inline bool writeLogCsv(const std::string& path, const std::vector<ExperimentLogRow>& rows, std::size_t count)
    {
        std::ofstream f(path);
        if (!f)
            return false;
        f.precision(10);
        writeLogHeader(f);
        for (std::size_t k = 0; k < count && k < rows.size(); ++k)
            writeLogRow(f, rows[k]);
        return static_cast<bool>(f);
    }

}  // namespace se3pf::robot
