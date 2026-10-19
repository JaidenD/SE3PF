#include "mujoco_experiments.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace se3pf {

    // ---------------------------------------------------------------------------------------------------------
    // MujocoDynamics

    MujocoDynamics::MujocoDynamics(const std::string& xmlPath) : arm_(std::make_unique<MujocoPanda>(xmlPath)) {}

    Mat<7, 7> MujocoDynamics::mass(const Vec7& q) const
    {
        arm_->setState(q, Vec7::Zero());
        Mat<7, 7> M = Mat<7, 7>::Zero();
        Vec7 bias = Vec7::Zero(), passive = Vec7::Zero();
        arm_->dynamics(M, bias, passive);
        return M;
    }

    Vec7 MujocoDynamics::coriolis(const Vec7& q, const Vec7& dq) const
    {
        arm_->setState(q, dq);
        Mat<7, 7> M = Mat<7, 7>::Zero();
        Vec7 bias = Vec7::Zero(), passive = Vec7::Zero();
        arm_->dynamics(M, bias, passive);
        return bias - arm_->gravity(q);  // qfrc_bias = C(q, qd) qd + G(q)
    }

    Vec7 MujocoDynamics::gravity(const Vec7& q) const { return arm_->gravity(q); }

    // ---------------------------------------------------------------------------------------------------------
    // The closed loop

    namespace {
        // libfranka's treatment of a torque command (ControlLoop::spinControl, 0.9.2): a first-order low-pass filter
        // at 100 Hz against the last command, then the rate limiter at 1000 - 1e-3 N m/s, both per 1 ms.
        Vec7 libfrankaCommand(const Vec7& tau, const Vec7& last)
        {
            const double dt = 1e-3, cutoff = 100.0, maxRate = 1000.0 - 1e-3;
            const double gain = dt / (dt + 1.0 / (2.0 * EIGEN_PI * cutoff));
            Vec7 out = Vec7::Zero();
            for (int i = 0; i < 7; ++i) {
                const double filtered = gain * tau[i] + (1.0 - gain) * last[i];
                out[i] = last[i] + std::clamp(filtered - last[i], -maxRate * dt, maxRate * dt);
            }
            return out;
        }
    }  // namespace

    MujocoExperimentRun simulateExperimentMujoco(const ExperimentSetup& setup, const SafetyMonitor& monitor,
                                                 const MujocoExperimentOptions& options)
    {
        if (options.plant.empty() || options.substeps < 1)
            throw std::invalid_argument("simulateExperimentMujoco: need a plant model and substeps >= 1");

        const double dt = RunConstants::dt;
        MujocoPanda::Options plantOptions;
        plantOptions.timestep = dt / options.substeps;
        plantOptions.integrator = MujocoPanda::Integrator::RK4;
        plantOptions.damping = options.damping;
        MujocoPanda plant(options.plant, plantOptions);
        std::unique_ptr<MujocoPanda> separate;
        if (!options.model.empty()) {
            MujocoPanda::Options modelOptions = plantOptions;
            modelOptions.damping = false;
            separate = std::make_unique<MujocoPanda>(options.model, modelOptions);
        }
        MujocoPanda& model = separate ? *separate : plant;

        MujocoExperimentRun run;
        run.peaks.torques = true;
        const long maxPeriods =
            static_cast<long>(std::ceil((setup.timeout() + RunConstants::stopTime + 1.0) / dt));
        run.log.reserve(static_cast<std::size_t>(maxPeriods));

        plant.setState(setup.start(), Vec7::Zero());
        ExperimentLoop loop(setup, monitor);
        ExperimentLoop::Cycle cycle;
        ModelTerms terms;
        Mat<7, 7> Mp = Mat<7, 7>::Zero();
        Vec7 biasP = Vec7::Zero(), passiveP = Vec7::Zero(), bias = Vec7::Zero(), passive = Vec7::Zero();
        Vec7 sent = Vec7::Zero();  // the command after libfranka's filter (RobotState::tau_J_d)
        bool stopping = false;
        SafetyStop stop;
        double cruiseSum = 0.0;
        long cruiseCount = 0;
        const SpeedProfile& speed = setup.spec().speed;
        const bool velocity = setup.spec().mode == TangentialReference::Mode::Velocity;

        for (long k = 0; k < maxPeriods; ++k) {
            const double t = k * dt;
            const Vec7 q = plant.q(), dq = plant.qd();

            // The robot's model at the measured state: M and C qd of the model, G of the plant (configured load)
            plant.dynamics(Mp, biasP, passiveP);
            run.maxContacts = std::max(run.maxContacts, plant.contacts());
            if (separate) {
                model.setState(q, dq);
                model.dynamics(terms.M, bias, passive);
                const Vec7 Gm = model.gravity(q);
                terms.c = bias - Gm;
                terms.G = options.plantGravity ? plant.gravity(q) : Gm;
            } else {
                terms.M = Mp;
                terms.G = plant.gravity(q);
                terms.c = biasP - terms.G;
            }

            Vec7 tau = Vec7::Zero(), qdd = Vec7::Zero(), tauF = Vec7::Zero();
            ExperimentLoop::Status st;
            bool done = false;
            if (!stopping) {
                loop.follow(t, q, dq, &terms, cycle);
                st = cycle.status;
                qdd = cycle.qdd;
                if (cycle.violation != Violation::None) {
                    stopping = true;
                    run.violation = cycle.violation;
                    run.tStop = t;
                    stop.begin(t, q);
                } else {
                    tau = cycle.tau;
                    tauF = cycle.tauF;
                    done = cycle.finished;
                    run.finished = done;
                }
            }
            if (stopping) {
                // The driver's safety stop: the damped hold at q_stop
                tau = stop.step(t, q, dq, done);
                if (t > run.tStop)
                    loop.observe(q, dq, st);
                run.maxStopDrift = std::max(run.maxStopDrift, (q - stop.holdConfiguration()).cwiseAbs().maxCoeff());
                run.stopSettled = stop.settled();
            }
            sent = options.filter ? libfrankaCommand(tau, sent) : tau;
            const Vec7 delivered = sent + terms.G;  // the robot adds G
            // Logged as the driver logs it: the command before libfranka's filter, and that command plus G
            ExperimentLogRow row = makeLogRow(t, stopping ? 1 : 0, st, q, dq, qdd, tau, Vec7(tau + terms.G), tauF);
            row.periodMs = 1.0;
            run.peaks.add(row);
            run.log.push_back(row);
            if (!stopping && t >= std::max(0.0, setup.duration() - 5.0)) {
                run.steadyXi = std::max(run.steadyXi, row.xiNorm);
                run.steadyBeta = std::max(run.steadyBeta, std::abs(row.beta));
            }
            if (!stopping && velocity && t >= speed.rampUp && t < speed.rampUp + speed.cruise) {
                cruiseSum += st.etaDot;
                ++cruiseCount;
            }
            if (!plant.withinTorqueLimits(delivered))
                ++run.clampedPeriods;
            run.tEnd = t;
            run.endSpeed = dq.norm();
            if (done)
                break;

            Vec7 push = Vec7::Zero();
            if (options.pushJoint >= 0 && t >= options.pushStart && t < options.pushStart + options.pushDuration)
                push[options.pushJoint] = options.pushTorque;
            plant.setAppliedTorque(push);
            plant.holdStep(delivered, options.substeps);
            if (plant.warnings() > 0)
                break;
        }
        run.cruiseEtaDot = cruiseCount > 0 ? cruiseSum / cruiseCount : 0.0;
        run.warnings = plant.warnings();
        return run;
    }

    MujocoStopRun simulateStopMujoco(const std::string& plantFile, const std::string& gravityModel, const Vec7& q0,
                                     const Vec7& qd0)
    {
        const double dt = RunConstants::dt;
        MujocoPanda::Options options;
        options.timestep = dt;
        options.integrator = MujocoPanda::Integrator::RK4;
        options.damping = true;
        MujocoPanda plant(plantFile, options);
        std::unique_ptr<MujocoPanda> separate;
        if (!gravityModel.empty()) {
            MujocoPanda::Options modelOptions = options;
            modelOptions.damping = false;
            separate = std::make_unique<MujocoPanda>(gravityModel, modelOptions);
        }
        MujocoPanda& gModel = separate ? *separate : plant;

        plant.setState(q0, qd0);
        SafetyStop stop;
        stop.begin(0.0, q0);
        MujocoStopRun run;
        run.minJointMargin = std::numeric_limits<double>::infinity();
        Vec7 sent = Vec7::Zero();
        const long maxPeriods = static_cast<long>(std::ceil(RunConstants::stopTime / dt)) + 2;
        for (long k = 0; k < maxPeriods; ++k) {
            const double t = k * dt;
            const Vec7 q = plant.q(), dq = plant.qd();
            bool done = false;
            const Vec7 tau = stop.step(t, q, dq, done);
            run.maxTorque = std::max(run.maxTorque, tau.cwiseAbs().maxCoeff());
            run.drift = run.drift.cwiseMax((q - q0).cwiseAbs());
            run.minJointMargin = std::min(run.minJointMargin, PandaArm::jointMargin(q));
            run.duration = t;
            run.endSpeed = dq.norm();
            run.settled = stop.settled();
            if (done)
                break;
            sent = libfrankaCommand(tau, sent);
            plant.holdStep(Vec7(sent + gModel.gravity(q)), 1);
            if (plant.warnings() > 0)
                break;
        }
        run.warnings = plant.warnings();
        return run;
    }

}  // namespace se3pf
