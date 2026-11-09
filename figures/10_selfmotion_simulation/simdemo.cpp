// The self-motion demos in simulation: exact model (qdd = qdd*), no friction, no integral action
#include "se3pf/experiments.hpp"
#include "se3pf/preflight.hpp"
#include "robot_log.hpp"
#include <cstdio>
using namespace se3pf;
int main(int argc, char** argv)
{
    const std::string dir = argv[1];
    PandaExamples ex;
    for (const char* name : {"orbit-selfmotion-hold", "orbit-selfmotion-moving"}) {
        ExperimentSpec e = *findExperiment(name);
        e.gains.kiEta = e.gains.Ki = e.gains.KiBeta = 0.0;  // the paper's control law: PD only
        e.rampTime = 0.0;  // no fade-in of mu (a hardware measure): the designed response from t = 0
        const ExperimentSetup s(ex(e), e);
        std::vector<ExperimentLogRow> log;
        const PreflightReport r = preflight(s, nullptr, s.monitor(), &log);
        std::printf("%s: %s\n    %s\n", name, r.ok ? "OK" : r.reason.c_str(), r.summary().c_str());
        const std::string path = dir + "/" + std::string(name) + "_sim.csv";
        robot::writeLogCsv(path, log, log.size());
        robot::writeLogCsv(path + ".preflight.csv", log, log.size());  // the lift for the plot
    }
}
