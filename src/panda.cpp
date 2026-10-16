#include "se3pf/panda.hpp"
#include "se3pf/se3.hpp"

#include <cmath>

using namespace Eigen;

namespace se3pf {
    namespace {
        // Modified DH parameters of joint i: a_{i-1}, d_i and alpha_{i-1}.
        // Every alpha is 0, -pi/2 or pi/2, so we store its cosine and sine, which are exactly 0 or +-1: no trig
        // on constants in the 1 kHz loop, and the frames at q = 0 come out exact.
        struct DHRow {
            double a, d, cosAlpha, sinAlpha;
        };

        constexpr DHRow kDH[7] = {
            //  a_{i-1}   d_i    cos   sin         alpha_{i-1}
            {0.0,     0.333, 1.0,  0.0},   // joint 1:  0
            {0.0,     0.0,   0.0, -1.0},   // joint 2: -pi/2
            {0.0,     0.316, 0.0,  1.0},   // joint 3:  pi/2
            {0.0825,  0.0,   0.0,  1.0},   // joint 4:  pi/2
            {-0.0825, 0.384, 0.0, -1.0},   // joint 5: -pi/2
            {0.0,     0.0,   0.0,  1.0},   // joint 6:  pi/2
            {0.088,   0.0,   0.0,  1.0},   // joint 7:  pi/2
        };
    }  // namespace

    // Forward kinematics T = T_1(q_1) ... T_7(q_7) Trans_z(0.107), with the modified DH transform
    //   T_i(th) = [ cos th             -sin th              0            a_{i-1}
    //               sin th cos alpha    cos th cos alpha   -sin alpha   -d_i sin alpha
    //               sin th sin alpha    cos th sin alpha    cos alpha    d_i cos alpha ],  alpha = alpha_{i-1}.
    // The frames are built once, as a running product (R, p), and give the Jacobian on the way.
    Pose PandaArm::output(const Vec7 &q, Mat67 *Jb)
    {
        Matrix3d R = Matrix3d::Identity();  // attitude of joint frame i in the base frame
        Vector3d p = Vector3d::Zero();      // its origin
        Mat67 Js = Mat67::Zero();           // spatial Jacobian, rows (omega, v)

        for (int i = 0; i < kDof; ++i) {
            const DHRow &dh = kDH[i];
            const double c = std::cos(q[i]), s = std::sin(q[i]);

            Matrix3d Ri;  // rotation of T_i(q_i), every entry set below
            Ri << c,                -s,                0.0,
                  s * dh.cosAlpha,   c * dh.cosAlpha, -dh.sinAlpha,
                  s * dh.sinAlpha,   c * dh.sinAlpha,  dh.cosAlpha;
            const Vector3d ti(dh.a, -dh.d * dh.sinAlpha, dh.d * dh.cosAlpha);  // translation of T_i

            // (R, p) <- (R, p) T_i = (R Ri, p + R ti): the translation needs the old R, so update p first.
            // R = R * Ri is safe: Eigen evaluates a matrix product into a temporary before assigning.
            p += R * ti;
            R = R * Ri;

            // Joint i turns about the z-axis z_i of its frame, through its origin p_i. A rotation with unit rate
            // about that axis has spatial velocity (omega, v) = (z_i, p_i x z_i).
            const Vector3d z = R.col(2);
            Js.col(i) << z, p.cross(z);
        }

        p += kFlangeOffset * R.col(2);  // Trans_z(0.107): the flange has joint 7's attitude
        const Pose g = makePose(R, p);

        // Body Jacobian from the spatial one: zeta_b = Ad_{g^-1} zeta_s.
        if (Jb)
            *Jb = Ad(inverse(g)) * Js;
        return g;
    }

    Vec7 PandaArm::retract(const Vec7 &q, const Vec7 &dq, double h)
    {
        return q + h * dq;
    }

    Vec7 PandaArm::restPosture()
    {
        const double pi = EIGEN_PI;
        return (Vec7() << 0.0, -pi / 4, 0.0, -3 * pi / 4, 0.0, pi / 2, pi / 4).finished();
    }

    double PandaArm::jointMargin(const Vec7 &q, int *joint)
    {
        // Distance of each joint to its nearer limit; positive inside the range, negative outside.
        const Vec7 margins = (q - qMin()).cwiseMin(qMax() - q);
        // Fail closed: the default minCoeff (PropagateFast) skips a NaN unless it is the first entry, which would
        // let a corrupted q pass a margin check. PropagateNaN returns NaN with the index of the NaN joint.
        int i = 0;
        const double margin = margins.minCoeff<PropagateNaN>(&i);
        if (joint)
            *joint = i;
        return margin;
    }

    Vec7 PandaArm::qMin()
    {
        return (Vec7() << -2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973).finished();
    }

    Vec7 PandaArm::qMax()
    {
        return (Vec7() << 2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973).finished();
    }

    Vec7 PandaArm::dqMax()
    {
        return (Vec7() << 2.1750, 2.1750, 2.1750, 2.1750, 2.6100, 2.6100, 2.6100).finished();
    }

    Vec7 PandaArm::ddqMax()
    {
        return (Vec7() << 15.0, 7.5, 10.0, 12.5, 15.0, 20.0, 20.0).finished();
    }

    Vec7 PandaArm::jerkMax()
    {
        return (Vec7() << 7500.0, 3750.0, 5000.0, 6250.0, 7500.0, 10000.0, 10000.0).finished();
    }

    Vec7 PandaArm::tauMax()
    {
        return (Vec7() << 87.0, 87.0, 87.0, 87.0, 12.0, 12.0, 12.0).finished();
    }

}  // namespace se3pf
