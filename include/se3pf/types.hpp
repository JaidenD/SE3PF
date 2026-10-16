#pragma once
#include <Eigen/Dense>

namespace se3pf {
    using Vec5 = Eigen::Matrix<double, 5, 1>;
    using Vec6 = Eigen::Matrix<double, 6, 1>;
    using Vec7 = Eigen::Matrix<double, 7, 1>;
    using Mat5 = Eigen::Matrix<double, 5, 5>;
    using Mat6 = Eigen::Matrix<double, 6, 6>;
    using Mat65 = Eigen::Matrix<double, 6, 5>;  // normal frame N
    using Mat67 = Eigen::Matrix<double, 6, 7>;  // Panda body Jacobian
    using Pose = Eigen::Matrix4d;

    // Sizes that depend on the system (n = 6 or 7)
    template <int R, int C> using Mat = Eigen::Matrix<double, R, C>;
    template <int N> using Vec = Eigen::Matrix<double, N, 1>;

    using Eigen::Matrix3d;
    using Eigen::Matrix4d;
    using Eigen::Vector3d;
}
