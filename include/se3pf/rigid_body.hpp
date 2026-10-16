#pragma once
#include "se3pf/types.hpp"

namespace se3pf {

    // The free rigid body (Section III, Free rigid bodies). Configuration g = (R, p) in SE(3), generalized velocity
    // the body twist zeta = (omega, v) at the centre of mass, input the body wrench u = (tau, f). Euler's equations
    // in body coordinates,
    //     J omega' = J omega x omega + tau,        m v' = m v x omega + f,
    // are eq:EL with M = diag(J, m I), C(g, zeta) zeta = -(J omega x omega, m v x omega) and G = 0.
    //
    // As a System for the controller: the output is the pose itself (Pi = id), so the body Jacobian is J_b = I.
    // Everything except the constructor is fixed-size and allocation-free, and does not throw.
    class RigidBody {
        public:
            static constexpr int kDof = 6;
            using Config = Pose;

            // J: inertia about the centre of mass in body coordinates, m: mass. Throws std::invalid_argument
            // unless J is finite, symmetric positive definite and m is positive and finite. A J that is symmetric
            // only up to rounding (relative 1e-12) is accepted, and its symmetric part (J + J^T) / 2 is used.
            RigidBody(const Matrix3d& J, double m);
            static RigidBody example();  // the paper's example: J = diag(1.0, 1.2, 0.8) kg m^2, m = 2 kg

            // System interface (see the design doc, The system interface)
            static Pose output(const Pose& g, Mat6* Jb);                     // Pi(g) = g; *Jb = I if Jb != nullptr
            static Pose retract(const Pose& g, const Vec6& zeta, double h);  // g Exp_I(h zeta), product metric

            // Dynamics, eq:EL (used by the simulator, not by the controller)
            Mat6 mass() const;                                      // M = diag(J, m I)
            Vec6 coriolis(const Vec6& zeta) const;                  // C(g, zeta) zeta
            Vec6 acceleration(const Vec6& zeta, const Vec6& u) const;  // zeta' = M^-1 (u - C zeta), since G = 0

            // One RK4 step of length h with u held constant; updates g and zeta in place. R stays in SO(3).
            void integrate(Pose& g, Vec6& zeta, const Vec6& u, double h) const;

        private:
            Matrix3d J_;     // (J + J^T) / 2: exactly symmetric
            Matrix3d Jinv_;  // J_^-1, computed once in the constructor
            double m_;
    };

}  // namespace se3pf
