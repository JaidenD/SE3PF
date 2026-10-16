#include "se3pf/preflight.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>

namespace se3pf {

    std::string PreflightReport::summary() const
    {
        std::string s;
        char head[300];
        if (ok)
            std::snprintf(head, sizeof(head), "OK (%.2f s simulated%s)", simulated,
                          torqueChecked ? "" : ", kinematic only: torques NOT checked");
        else
            std::snprintf(head, sizeof(head), "FAILED: %s at t = %.3f s (%.2f s simulated%s)", reason.c_str(), tFail,
                          simulated, torqueChecked ? "" : ", kinematic only");
        s = head;
        s += "\n    ";
        s += peaks.summary(limits);
        return s;
    }

    PreflightReport preflight(const ExperimentSetup& setup, const DynamicsModel* dyn, const SafetyMonitor& monitor,
                              std::vector<ExperimentLogRow>* log, const Vec7* start)
    {
        PreflightReport rep;
        rep.torqueChecked = dyn != nullptr;
        rep.limits = monitor.limits();
        rep.peaks.torques = dyn != nullptr;

        ExperimentLoop loop(setup, monitor);
        const double dt = RunConstants::dt;
        // The timeout ends every run: at most (T + timeoutExtra) / dt + 2 cycles
        const long maxCycles = static_cast<long>(std::ceil(setup.timeout() / dt)) + 2;
        if (log) {
            log->clear();
            log->reserve(static_cast<std::size_t>(maxCycles));
        }

        Vec7 q = start ? *start : setup.start();
        Vec7 dq = Vec7::Zero();
        ModelTerms terms;
        ExperimentLoop::Cycle cycle;
        bool finished = false;
        for (long k = 0; k < maxCycles; ++k) {
            const double t = k * dt;
            if (dyn) {
                terms.M = dyn->mass(q);
                terms.c = dyn->coriolis(q, dq);
                terms.G = dyn->gravity(q);
            }
            loop.follow(t, q, dq, dyn ? &terms : nullptr, cycle);
            rep.simulated = t;
            if (cycle.violation == Violation::Projection || cycle.violation == Violation::NonFinite) {
                rep.violation = cycle.violation;
                rep.tFail = t;
                break;
            }
            const ExperimentLogRow row = makeLogRow(t, 0, cycle.status, q, dq, cycle.qdd, cycle.tau, cycle.tauG);
            rep.peaks.add(row);
            if (log)
                log->push_back(row);
            if (k == 0) {
                rep.startXi = cycle.status.xi.norm();
                rep.startBeta = cycle.status.beta;
            }
            if (cycle.violation != Violation::None) {
                rep.violation = cycle.violation;
                rep.tFail = t;
                break;
            }
            if (cycle.finished) {
                finished = true;
                break;
            }
            // eq:fb with an exact model gives qdd = qdd*, held over the period: exact for a constant acceleration
            q += dt * dq + (0.5 * dt * dt) * cycle.qdd;
            dq += dt * cycle.qdd;
        }

        rep.converged = rep.peaks.rows > 0 && std::abs(rep.peaks.finalXi) <= PreflightReport::kConvergence
                        && std::abs(rep.peaks.finalBeta) <= PreflightReport::kConvergence;
        if (rep.violation != Violation::None) {
            rep.reason = toString(rep.violation);
        } else if (!finished) {
            // Cannot happen (the timeout ends the loop first); kept so that a change of the loop fails closed.
            rep.violation = Violation::Timeout;
            rep.reason = toString(Violation::Timeout);
            rep.tFail = rep.simulated;
        } else if (!rep.converged) {
            rep.reason = "the outputs did not converge (final |xi| or |beta| above 1e-3)";
            rep.tFail = rep.simulated;
        }
        rep.ok = rep.reason.empty();
        return rep;
    }

}  // namespace se3pf
