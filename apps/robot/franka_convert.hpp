#pragma once
// Conversions between libfranka's arrays and the se3pf types, for the robot driver (apps/robot/robot_driver.cpp).
// No libfranka header is needed: libfranka passes std::array, and these functions take and return those arrays, so
// they are unit tested without libfranka (tests/test_franka_convert.cpp).
//
// libfranka's conventions (include/franka/model.h, robot_state.h of 0.9.2): every matrix is a column-major array;
// a 6 x 7 Jacobian has the rows (v, omega), translational first, while se3pf orders twists (omega, v) (design doc,
// Conventions), so the two 3-row blocks are swapped; I_total is a column-major 3 x 3 matrix.
#include "se3pf/types.hpp"

#include <array>
#include <cstddef>

namespace se3pf::robot {

    inline Vec7 vec7(const std::array<double, 7>& a)
    {
        Vec7 v = Vec7::Zero();
        for (int i = 0; i < 7; ++i)
            v[i] = a[static_cast<std::size_t>(i)];
        return v;
    }

    inline std::array<double, 7> array7(const Vec7& v)
    {
        std::array<double, 7> a{};
        for (int i = 0; i < 7; ++i)
            a[static_cast<std::size_t>(i)] = v[i];
        return a;
    }

    inline Vector3d vec3(const std::array<double, 3>& a) { return Vector3d(a[0], a[1], a[2]); }

    // A homogeneous transform from libfranka's column-major array (model.pose, RobotState::O_T_EE, F_T_EE, ...)
    inline Pose poseFromColumnMajor(const std::array<double, 16>& a)
    {
        Pose g = Pose::Zero();
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                g(r, c) = a[static_cast<std::size_t>(c * 4 + r)];
        return g;
    }

    // M(q) from model.mass: 7 x 7, column-major (symmetric, but the layout is converted all the same)
    inline Mat<7, 7> mat7FromColumnMajor(const std::array<double, 49>& a)
    {
        Mat<7, 7> M = Mat<7, 7>::Zero();
        for (int r = 0; r < 7; ++r)
            for (int c = 0; c < 7; ++c)
                M(r, c) = a[static_cast<std::size_t>(c * 7 + r)];
        return M;
    }

    // 3 x 3 from a column-major array (I_total, I_ee, I_load)
    inline Matrix3d mat3FromColumnMajor(const std::array<double, 9>& a)
    {
        Matrix3d M = Matrix3d::Zero();
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                M(r, c) = a[static_cast<std::size_t>(c * 3 + r)];
        return M;
    }

    // The body Jacobian of model.bodyJacobian(Frame::kFlange, ...) in se3pf's convention: column-major 6 x 7 with
    // rows (v, omega) becomes rows (omega, v), so that zeta = (omega, v) = J_b qd as PandaArm::output gives it.
    inline Mat67 bodyJacobianFromFranka(const std::array<double, 42>& a)
    {
        Mat67 J = Mat67::Zero();
        for (int c = 0; c < 7; ++c)
            for (int r = 0; r < 6; ++r) {
                const int row = r < 3 ? r + 3 : r - 3;  // franka row r (v_x, v_y, v_z, w_x, w_y, w_z)
                J(row, c) = a[static_cast<std::size_t>(c * 6 + r)];
            }
        return J;
    }

}  // namespace se3pf::robot
