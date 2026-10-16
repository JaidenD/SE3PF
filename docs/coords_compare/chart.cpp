// Riemannian tubular chart (the paper: Psi = gamma Exp_I(N xi)) versus the screw chart (gamma exp(N xi), the matrix
// exponential): how far apart are their tangential and transverse outputs, and does |xi| still equal the distance?
#include "se3pf/panda_example.hpp"
#include "se3pf/rigid_body_example.hpp"
#include "se3pf/se3.hpp"

#include <cmath>
#include <cstdio>
#include <random>

using namespace se3pf;

namespace {
    // Matrix logarithm of SE(3) (screw coordinates): omega = log R, v = J_l(omega)^-1 p, with J_l(w) = J_r(-w).
    Vec6 screwLog(const Pose& g)
    {
        const Vector3d w = logSO3(Matrix3d(g.block<3, 3>(0, 0)));
        const Vector3d v = JrInv(Vector3d(-w)) * Vector3d(g.block<3, 1>(0, 3));
        Vec6 z;
        z << w, v;
        return z;
    }

    double wrapSym(double x, double L)
    {
        double y = std::fmod(x, L);
        if (y > 0.5 * L) y -= L;
        if (y <= -0.5 * L) y += L;
        return y;
    }

    // Projection condition of the screw chart: g = gamma(s) exp(N(s) xi)  <=>  <log(gamma(s)^-1 g), zeta(s)>_W = 0.
    double phiScrew(const Tube& tube, const Pose& g, double s)
    {
        const Vec6 e = screwLog(Pose(inverse(tube.path().pose(s)) * g));
        return tube.metric().inner(e, tube.path().velocity(s));
    }

    struct Cmp {
        double dEta, xiScrew, dist, relDist, tcpMix, tcp;
    };

    // Both charts' outputs for pose g; guess near the Riemannian foot point.
    Cmp compare(const Tube& tube, const Pose& g, double guess)
    {
        double eta;
        Vec5 xi;
        if (!tube.outputs(g, guess, eta, xi)) {
            std::printf("Riemannian projection failed\n");
            return {NAN, NAN, NAN, NAN, NAN, NAN};
        }
        // Newton with a central-difference derivative on the screw chart's projection condition, from eta
        double s = eta;
        for (int k = 0; k < 50; ++k) {
            const double h = 1e-7, f = phiScrew(tube, g, s);
            const double df = (phiScrew(tube, g, s + h) - phiScrew(tube, g, s - h)) / (2 * h);
            const double ds = f / df;
            s -= ds;
            if (std::abs(ds) < 1e-15) break;
        }
        Mat65 N;
        tube.frame(s, N);
        const Pose gErrScrew = inverse(tube.path().pose(s)) * g;
        const Vec6 eScrew = screwLog(gErrScrew);
        const Vec5 xiScrew = N.transpose() * tube.metric().matrix() * eScrew;

        // the TCP displacement from the screw foot point (product-metric translation) vs the screw chart's v
        const Vector3d tcp = gErrScrew.block<3, 1>(0, 3);
        const double dist = xi.norm();
        return {wrapSym(s - eta, tube.length()), xiScrew.norm(), dist, (xiScrew.norm() - dist) / dist,
                (eScrew.tail<3>() - tcp).norm(), tcp.norm()};
    }

    void sweep(const char* name, const Tube& tube, double radius, double s0, const Vec5& xi0)
    {
        const Cmp c = compare(tube, tube.Psi(s0, xi0), s0);
        std::printf("%s, start pose Psi(%.1f, xi0), |xi0| = %.4f:\n", name, s0, xi0.norm());
        std::printf("  foot point eta_screw - eta = %+.3e   |xi_screw| = %.6f vs distance %.6f (%+.2f %%)\n", c.dEta,
                    c.xiScrew, c.dist, 100 * c.relDist);
        std::printf("  screw chart's 'translation' differs from the TCP displacement by %.3e (|TCP disp.| %.3e, %.2f %%)\n",
                    c.tcpMix, c.tcp, 100 * c.tcpMix / c.tcp);

        std::mt19937 rng(7);
        std::normal_distribution<double> nd;
        std::uniform_real_distribution<double> ud(0.0, tube.length());
        for (double frac : {0.1, 0.25, 0.5, 0.75}) {
            double maxEta = 0, maxRel = 0;
            for (int i = 0; i < 400; ++i) {
                Vec5 d;
                for (int j = 0; j < 5; ++j) d[j] = nd(rng);
                d.normalize();
                const double s = ud(rng);
                const Cmp ci = compare(tube, tube.Psi(s, frac * radius * d), s);
                maxEta = std::max(maxEta, std::abs(ci.dEta));
                maxRel = std::max(maxRel, std::abs(ci.relDist));
            }
            std::printf("  400 random poses at |xi| = %.2f x tube radius (%.4f): max |eta_screw - eta| = %.2e (%.3f %% of L), "
                        "max ||xi_screw| - dist| / dist = %.2f %%\n",
                        frac, frac * radius, maxEta, 100 * maxEta / tube.length(), 100 * maxRel);
        }
    }
}

int main()
{
    const PandaExample panda;            // l = 0.1 m, W = diag(0.01 I, I), tube radius 0.125
    const RigidBodyExample rigid;        // W = I, tube radius 1.93
    sweep("Panda (l = 0.1 m)", panda.tube(), 0.125, PandaExample::s0, PandaExample::xi0());
    sweep("Rigid body (W = I)", rigid.tube(), 1.93, RigidBodyExample::kS0, RigidBodyExample::xi0());
}
