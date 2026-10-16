#pragma once
#include "se3pf/types.hpp"

namespace se3pf {
    // SO(3)
    Matrix3d hat(const Vector3d& w);
    Vector3d vee(const Matrix3d& S);    // S skew-symmetric

    Matrix3d expSO3(const Vector3d& w);
    Vector3d logSO3(const Matrix3d& R);

    Matrix3d Jr(const Vector3d& w);     // right Jacobian: d/dt exp(w) = exp(w) (Jr(w) wdot)^
    Matrix3d JrInv(const Vector3d& w);

    // SE(3), twists ordered (omega, v)
    Matrix4d hat(const Vec6& z);        // [hat(omega) v; 0 0], an element of se(3)
    Vec6 vee(const Matrix4d& X);

    Pose makePose(const Matrix3d& R, const Vector3d& p);
    inline Matrix3d rotation(const Pose& g) { return g.block<3, 3>(0, 0); }     // R of g = (R, p)
    inline Vector3d translation(const Pose& g) { return g.block<3, 1>(0, 3); }  // p of g = (R, p)
    Pose inverse(const Pose& g);        // (R^T, -R^T p), in closed form
    Mat6 Ad(const Pose& g);             // Ad_g z = (g hat(z) g^-1)^vee
    Mat6 ad(const Vec6& z);             // ad_a b = (hat(a) hat(b) - hat(b) hat(a))^vee
}
