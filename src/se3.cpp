#include "se3pf/se3.hpp"
#include <cmath>

using namespace Eigen;

namespace se3pf {

    Matrix3d hat(const Vector3d &w)
    {
        Matrix3d e1 = Matrix3d::Zero(), e2 = Matrix3d::Zero(), e3 = Matrix3d::Zero();

        e1(1,2) = -1.0;
        e1(2,1) = 1.0;

        e2(2,0) = -1.0;
        e2(0,2) = 1.0;

        e3(0,1) = -1.0;
        e3(1,0) = 1.0;

        return w[0] * e1 + w[1] * e2 + w[2] * e3;
    }

    Vector3d vee(const Matrix3d &S)
    {
        Vector3d v = Vector3d::Zero();
        v(0) = S(2,1);
        v(1) = S(0,2);
        v(2) = S(1,0);
        return v;
    }

    Matrix3d expSO3(const Vector3d &w)
    {
        const double theta = w.norm();
        const Matrix3d K = hat(w);

        // for small angles: sin(theta)/theta ~ 1, (1-cos(theta))/theta^2 ~ 1/2
        if (theta < 1e-8)
            return Matrix3d::Identity() + K + 0.5 * K * K;

        // Rodrigues' formula
        return Matrix3d::Identity() + (std::sin(theta) / theta) * K
             + ((1.0 - std::cos(theta)) / (theta * theta)) * K * K;
    }

    Vector3d logSO3(const Matrix3d &R)
    {
        const Vector3d s = 0.5 * vee(Matrix3d(R - R.transpose()));  // = sin(theta) * axis
        const double c = 0.5 * (R.trace() - 1.0);          // = cos(theta)
        const double theta = std::atan2(s.norm(), c);         // angle in [0, pi], accurate at every angle

        if (theta < 1e-6)  // small angle: theta / sin(theta) = 1 + theta^2/6 + ...
            return (1.0 + theta * theta / 6.0) * s;

        if (theta > EIGEN_PI - 1e-4) {  // near pi, sin(theta) ~ 0, so s no longer carries the axis
            // symmetric part: (R + R^T)/2 = cos(theta) I + (1 - cos(theta)) n n^T
            const Matrix3d nnT = (0.5 * (R + R.transpose()) - c * Matrix3d::Identity()) / (1.0 - c);
            Index j;
            nnT.diagonal().maxCoeff(&j);           // largest |n_j| gives the most accurate column
            Vector3d n = nnT.col(j).normalized();
            if (n.dot(s) < 0.0)                    // sign of the axis from s (at exactly pi, either sign is valid)
                n = -n;
            return theta * n;
        }

        return (theta / s.norm()) * s;  // s.norm() = sin(theta)
    }

    Matrix3d Jr(const Vector3d &w)
    {
        const double theta = w.norm();
        const Matrix3d K = hat(w);

        // small angles: (1 - cos(theta))/theta^2 ~ 1/2, (theta - sin(theta))/theta^3 ~ 1/6
        if (theta < 1e-6)
            return Matrix3d::Identity() - 0.5 * K + (1.0 / 6.0) * K * K;

        return Matrix3d::Identity()
             - ((1.0 - std::cos(theta)) / (theta * theta)) * K
             + ((theta - std::sin(theta)) / (theta * theta * theta)) * K * K;
    }

    Matrix3d JrInv(const Vector3d &w)
    {
        const double theta = w.norm();
        const Matrix3d K = hat(w);

        // small angles: 1/theta^2 - (1 + cos(theta)) / (2 theta sin(theta)) ~ 1/12
        if (theta < 1e-6)
            return Matrix3d::Identity() + 0.5 * K + (1.0 / 12.0) * K * K;

        return Matrix3d::Identity()
             + 0.5 * K
             + (1.0 / (theta * theta) - (1.0 + std::cos(theta)) / (2.0 * theta * std::sin(theta))) * K * K;
    }

    Matrix4d hat(const Vec6 &z)
    {
        Matrix4d X = Matrix4d::Zero();
        X.block<3,3>(0,0) = hat(Vector3d(z.head<3>()));
        X.block<3,1>(0,3) = z.tail<3>();
        return X;
    }

    Vec6 vee(const Matrix4d &g)
    {
        Vec6 z = Vec6::Zero();
        z.head<3>() = vee(Matrix3d(g.block<3,3>(0,0)));
        z.tail<3>() = g.block<3,1>(0,3);
        return z;
    }

    Pose makePose(const Matrix3d &R, const Vector3d &p)
    {
        Pose g = Pose::Identity();

        g.block<3,3>(0,0) = R;
        g.block<3,1>(0,3) = p;

        return g;
    }

    Pose inverse(const Pose &g)
    {
        const Matrix3d Rt = g.block<3,3>(0,0).transpose();
        const Vector3d p = g.block<3,1>(0,3);
        return makePose(Rt, -Rt * p);
    }

    Mat6 ad(const Vec6 &z)
    {
        const Vector3d w = z.head<3>();
        const Vector3d v = z.tail<3>();

        Mat6 A = Mat6::Zero();
        A.block<3,3>(0,0) = hat(w);
        A.block<3,3>(3,3) = hat(w);
        A.block<3,3>(3,0) = hat(v);

        return A;
    }

    Mat6 Ad(const Pose &g)
    {
        const Matrix3d R = g.block<3,3>(0,0);
        const Vector3d p = g.block<3,1>(0,3);

        Mat6 A = Mat6::Zero();
        A.block<3,3>(0,0) = R;
        A.block<3,3>(3,3) = R;
        A.block<3, 3>(3, 0) = hat(p) * R;

        return A;
    }

}  // namespace se3pf
