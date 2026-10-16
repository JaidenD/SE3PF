#include "se3pf/rigid_body_example.hpp"
#include "se3pf/se3.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace se3pf {

    // Members in declaration order: the path needs the metric, the tube needs both.
    RigidBodyExample::RigidBodyExample()
        : metric_(Metric::product(1.0, 1.0)), path_(1.0, metric_), tube_(path_, metric_), body_(RigidBody::example())
    {
    }

    Vec5 RigidBodyExample::xi0()
    {
        Vec5 xi;
        xi << 0.4, -0.3, 0.2, -0.1, 0.3;
        return xi;
    }

    Pose RigidBodyExample::startPose() const
    {
        return tube_.Psi(kS0, xi0());
    }

    Gains RigidBodyExample::gains()
    {
        Gains g = Gains::fromLambda(kLambda);
        g.kt = 5.0;
        g.kpEta = 9.0;
        g.kdEta = 6.0;
        return g;
    }

    TangentialReference RigidBodyExample::reference(Reference which) const
    {
        const double L = path_.length();
        if (which == Reference::Position) {
            // eta_d(t) = L/4 + (L/5) (1 - cos(w t)), w = pi / 10: from L/4 to L/4 + 2L/5 and back in 20 s
            return TangentialReference::fromPosition([L](double t, double& etaD, double& detaD, double& ddetaD) {
                const double w = EIGEN_PI / 10.0;
                etaD = L / 4.0 + (L / 5.0) * (1.0 - std::cos(w * t));
                detaD = (L / 5.0) * w * std::sin(w * t);
                ddetaD = (L / 5.0) * w * w * std::cos(w * t);
            });
        }
        // nu_d(t) = 0.5 + 0.2 sin(0.5 t), independent of eta
        return TangentialReference::fromVelocity([](double t, double, double& nu, double& dnuDt, double& dnuDeta) {
            nu = 0.5 + 0.2 * std::sin(0.5 * t);
            dnuDt = 0.1 * std::cos(0.5 * t);
            dnuDeta = 0.0;
        });
    }

    double RigidBodyExample::designedXiNorm(double t)
    {
        return xi0().norm() * (1.0 + kLambda * t) * std::exp(-kLambda * t);
    }

    RigidBodyRun simulate(const RigidBodyExample& example, RigidBodyExample::Reference which, double T, double dt,
                          int substeps)
    {
        if (substeps < 1)
            throw std::invalid_argument("simulate: need substeps >= 1");
        // Checks T and dt (finite, T >= 0, dt > 0, T / dt <= kMaxSteps) before anything is built or allocated.
        const long steps = stepCount(T, dt, RigidBodyExample::kMaxSteps);  // the instants k dt <= T

        RigidBodyRun run;
        run.reference = which;
        run.dt = dt;
        run.substeps = substeps;

        Controller<RigidBody> controller(example.tube(), RigidBodyExample::gains(), example.reference(which));
        const RigidBody& body = example.body();
        const Mat6 M = body.mass();

        // Start at rest at g0 = Psi(0.1, xi0), whose closest point is gamma(0.1).
        Pose g = example.startPose();
        Vec6 zeta = Vec6::Zero();
        controller.reset(RigidBodyExample::kS0);

        const Pose gamma0 = example.path().pose(RigidBodyExample::kS0);
        const Matrix3d Rerr = gamma0.block<3,3>(0,0).transpose() * g.block<3,3>(0,0);
        run.startAttitudeDeg = logSO3(Rerr).norm() * 180.0 / EIGEN_PI;
        run.startPosition = (g.block<3,1>(0,3) - gamma0.block<3,1>(0,3)).norm();

        run.log.reserve(static_cast<std::size_t>(steps) + 1);
        const double h = dt / substeps;
        run.ok = true;
        for (long k = 0; k <= steps; ++k) {
            const double t = k * dt;
            Vec6 qdd = Vec6::Zero();
            Controller<RigidBody>::Status st;
            if (!controller.step(g, zeta, t, dt, qdd, st)) {
                run.ok = false;
                run.failureTime = t;
                break;
            }

            // eq:fb with G = 0, held over [t, t + dt) (zero-order hold)
            const Vec6 u = M * qdd + body.coriolis(zeta);

            RigidBodyLogRow row;
            row.t = t;
            row.eta = st.eta;
            row.etaDot = st.etaDot;
            row.ref = st.ref;
            row.eTan = st.eTan;
            row.xi = st.xi;
            row.xiNorm = st.xi.norm();
            row.xiDesign = RigidBodyExample::designedXiNorm(t);
            row.u = u;
            run.log.push_back(row);

            const double deviation = std::abs(row.xiNorm - row.xiDesign);
            if (deviation > run.maxDeviation) {
                run.maxDeviation = deviation;
                run.maxDeviationTime = t;
            }
            const double uMax = u.cwiseAbs().maxCoeff();
            if (uMax > run.maxInput) {
                run.maxInput = uMax;
                run.maxInputTime = t;
            }
            run.maxIterations = std::max(run.maxIterations, st.iterations);

            for (int i = 0; i < substeps; ++i)
                body.integrate(g, zeta, u, h);
        }

        if (!run.log.empty()) {
            run.finalXiNorm = run.log.back().xiNorm;
            run.finalTanError = run.log.back().eTan;
        }
        return run;
    }

}  // namespace se3pf
