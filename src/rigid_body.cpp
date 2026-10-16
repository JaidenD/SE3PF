#include "se3pf/rigid_body.hpp"
#include "se3pf/se3.hpp"

#include <Eigen/Cholesky>
#include <cmath>
#include <stdexcept>

namespace se3pf {

    namespace {
        // Time derivatives of the integrator's state (phi, p, zeta) within one step (see RigidBody::integrate).
        struct Rates {
            Vector3d phi;
            Vector3d p;
            Vec6 zeta;
        };
    }  // namespace

    // J_ is the symmetric part of J. A J computed elsewhere (Q D Q^T, or read from a file) can be asymmetric in its
    // last digits, so J is accepted when symmetric to a relative 1e-12, and from then on only J_ is used: mass(),
    // coriolis() and acceleration() (through the Cholesky factor of J_) then all describe the same matrix.
    RigidBody::RigidBody(const Matrix3d& J, double m)
        : J_(0.5 * (J + J.transpose())), Jinv_(Matrix3d::Identity()), m_(m)
    {
        if (!J.allFinite() || (J - J.transpose()).cwiseAbs().maxCoeff() > 1e-12 * J.cwiseAbs().maxCoeff())
            throw std::invalid_argument("RigidBody: the inertia J must be finite and symmetric");

        // A symmetric matrix has a Cholesky factorization exactly when it is positive definite.
        const Eigen::LLT<Matrix3d> llt(J_);
        if (llt.info() != Eigen::Success)
            throw std::invalid_argument("RigidBody: the inertia J must be positive definite");

        if (!(m > 0.0) || !std::isfinite(m))  // !(m > 0) also rejects NaN
            throw std::invalid_argument("RigidBody: the mass m must be positive and finite");

        Jinv_ = llt.solve(Matrix3d::Identity());
    }

    RigidBody RigidBody::example()
    {
        return RigidBody(Matrix3d(Vector3d(1.0, 1.2, 0.8).asDiagonal()), 2.0);
    }

    Pose RigidBody::output(const Pose& g, Mat6* Jb)
    {
        if (Jb != nullptr)
            *Jb = Mat6::Identity();  // Pi = id: the body velocity of the output is zeta itself
        return g;
    }

    // g Exp_I(h zeta) with the product-type Exp_I of Example 1: (R exp(h omega^), p + h R v). A curve through g
    // with velocity zeta, which is all the controller's central difference of A needs. The closed form does not
    // depend on w_R and w_p, so no Metric is needed.
    Pose RigidBody::retract(const Pose& g, const Vec6& zeta, double h)
    {
        const Matrix3d R = g.block<3,3>(0,0);
        const Vector3d p = g.block<3,1>(0,3);
        const Vector3d w = zeta.head<3>();
        const Vector3d v = zeta.tail<3>();
        return makePose(R * expSO3(Vector3d(h * w)), p + h * R * v);
    }

    Mat6 RigidBody::mass() const
    {
        Mat6 M = Mat6::Zero();
        M.block<3,3>(0,0) = J_;
        M.block<3,3>(3,3) = m_ * Matrix3d::Identity();
        return M;
    }

    // C(g, zeta) zeta = -(J omega x omega, m v x omega): moving Euler's equations (see the header) to the left
    // side of eq:EL flips their sign. Both parts are orthogonal to zeta, so they do no work.
    Vec6 RigidBody::coriolis(const Vec6& zeta) const
    {
        const Vector3d w = zeta.head<3>();
        const Vector3d v = zeta.tail<3>();
        Vec6 c;
        c << -(J_ * w).cross(w), -m_ * v.cross(w);
        return c;
    }

    // eq:EL with G = 0 solved for zeta'. M is block diagonal, so M^-1 = diag(J^-1, I / m).
    Vec6 RigidBody::acceleration(const Vec6& zeta, const Vec6& u) const
    {
        const Vec6 rhs = u - coriolis(zeta);
        Vec6 a;
        a << Jinv_ * rhs.head<3>(), rhs.tail<3>() / m_;
        return a;
    }

    // RK4 for zeta, and for the pose in exponential coordinates around the step's starting attitude Rk:
    //     R = Rk exp(phi^),    phi' = Jr(phi)^-1 omega,    p' = R v,    phi(0) = 0.
    // RK4 updates phi and p, which live in vector spaces; the final R = Rk exp(phi^) is a rotation by
    // construction, so R never leaves SO(3), however many steps are taken.
    void RigidBody::integrate(Pose& g, Vec6& zeta, const Vec6& u, double h) const
    {
        const Matrix3d Rk = g.block<3,3>(0,0);
        const Vector3d pk = g.block<3,1>(0,3);

        // Right-hand side at a stage. p does not appear on it, so it is not an argument.
        auto rates = [&](const Vector3d& phi, const Vec6& z) -> Rates {
            const Vector3d w = z.head<3>();
            const Vector3d v = z.tail<3>();
            Rates r;
            r.phi = JrInv(phi) * w;
            r.p = Rk * expSO3(phi) * v;
            r.zeta = acceleration(z, u);
            return r;
        };

        const Rates k1 = rates(Vector3d::Zero(), zeta);
        const Rates k2 = rates(0.5 * h * k1.phi, zeta + 0.5 * h * k1.zeta);
        const Rates k3 = rates(0.5 * h * k2.phi, zeta + 0.5 * h * k2.zeta);
        const Rates k4 = rates(h * k3.phi, zeta + h * k3.zeta);

        const double c = h / 6.0;
        const Vector3d phi = c * (k1.phi + 2.0 * k2.phi + 2.0 * k3.phi + k4.phi);
        const Vector3d p = pk + c * (k1.p + 2.0 * k2.p + 2.0 * k3.p + k4.p);
        zeta += c * (k1.zeta + 2.0 * k2.zeta + 2.0 * k3.zeta + k4.zeta);
        g = makePose(Rk * expSO3(phi), p);
    }

}  // namespace se3pf
