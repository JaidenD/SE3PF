// The simulated arm behind the mock franka::Robot and franka::Model (not a libfranka header). See
// mock/franka/README.md.
//
// Kinematics are the Panda's (se3pf::PandaArm). The dynamics are a toy arm, not the Panda's: constant diagonal joint
// inertias (rotors and links lumped together) plus rigid bodies fixed to the flange (a wrist mass, and the end-effector
// load that Desk configures). For rigid bodies at the flange the equations of motion are exact and consistent:
//     M(q) = diag(rotor) + J_b^T Lambda J_b,
//     C(q, qd) qd = J_b^T (Lambda Jdot_b qd - ad_V^T Lambda V),  V = J_b qd  (Newton-Euler in the flange frame),
//     G(q) = sum_k m_k (R (J_v - hat(c_k) J_w))^T (-g),
// with Lambda the bodies' spatial inertia at the flange (twists (omega, v)), so the torque loop, the computed torque,
// model errors and disturbances can be exercised end to end. The torque predictions that matter come from
// franka::Model on the robot and from MuJoCo.
#pragma once
#include "se3pf/panda.hpp"
#include "se3pf/se3.hpp"
#include "se3pf/types.hpp"

#include <Eigen/Dense>

#include <array>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>

namespace franka::mock {

    using Mat7 = se3pf::Mat<7, 7>;
    using se3pf::Mat6;
    using se3pf::Mat67;
    using se3pf::Vec6;
    using se3pf::Vec7;

    inline Vec7 fromArray(const std::array<double, 7>& a) { return Vec7(a.data()); }

    inline std::array<double, 7> toArray(const Vec7& v)
    {
        std::array<double, 7> a{};
        for (int i = 0; i < 7; ++i)
            a[static_cast<std::size_t>(i)] = v[i];
        return a;
    }

    // A rigid body fixed to the flange: mass, centre of mass and inertia about the centre of mass, in the flange
    // frame (Desk's end-effector parameters).
    struct RigidBody {
        double mass = 0.0;
        Eigen::Vector3d com = Eigen::Vector3d::Zero();
        Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();

        // Spatial inertia at the flange origin, twists (omega, v): X^T diag(I_c, m I) X with X = [I 0; -hat(c) I]
        Mat6 spatialInertia() const
        {
            const Eigen::Matrix3d C = se3pf::hat(com);
            Mat6 L = Mat6::Zero();
            L.topLeftCorner<3, 3>() = inertia - mass * C * C;
            L.topRightCorner<3, 3>() = mass * C;
            L.bottomLeftCorner<3, 3>() = -mass * C;
            L.bottomRightCorner<3, 3>() = mass * Eigen::Matrix3d::Identity();
            return L;
        }
    };

    // libfranka's load arrays: I_total is a column-major 3 x 3 matrix, F_x_Ctotal the centre of mass in the flange
    // frame.
    inline RigidBody bodyFromArrays(double m, const std::array<double, 3>& c, const std::array<double, 9>& I)
    {
        RigidBody b;
        b.mass = m;
        b.com = Eigen::Vector3d(c[0], c[1], c[2]);
        for (int r = 0; r < 3; ++r)
            for (int k = 0; k < 3; ++k)
                b.inertia(r, k) = I[static_cast<std::size_t>(k * 3 + r)];
        return b;
    }

    // The Franka Hand as Desk configures it by default: 0.73 kg, centre of mass (-0.01, 0, 0.03) m
    inline RigidBody frankaHand()
    {
        RigidBody b;
        b.mass = 0.73;
        b.com = Eigen::Vector3d(-0.01, 0.0, 0.03);
        b.inertia = Eigen::Vector3d(0.001, 0.0025, 0.0017).asDiagonal();
        return b;
    }

    class ToyArm {
        public:
            // Lumped joint inertias, kg m^2, chosen so that diag M near the start configurations is of the order of the
            // MuJoCo Menagerie Panda (0.91, 1.60, 1.04, 0.91, 0.13, 0.13, 0.11 with its 0.1 kg m^2 armature), and a
            // wrist mass at the flange
            static Vec7 rotor() { return (Vec7() << 0.6, 1.2, 0.6, 0.6, 0.12, 0.12, 0.1).finished(); }
            static RigidBody wrist()
            {
                RigidBody b;
                b.mass = 1.0;
                b.com = Eigen::Vector3d(0.0, 0.0, -0.05);
                b.inertia = 0.004 * Eigen::Matrix3d::Identity();
                return b;
            }

            // inertiaScale multiplies M (and with it C qd): a model with scale != 1 is a wrong model of the plant.
            explicit ToyArm(double inertiaScale = 1.0) : scale_(inertiaScale) {}

            // The equations of motion with the end-effector load `load` at the flange
            Mat7 mass(const Vec7& q, const RigidBody& load) const
            {
                Mat67 J = Mat67::Zero();
                se3pf::PandaArm::output(q, &J);
                const Mat6 L = lambda(load);
                return scale_ * (Mat7(rotor().asDiagonal()) + J.transpose() * L * J);
            }

            Vec7 coriolis(const Vec7& q, const Vec7& dq, const RigidBody& load) const
            {
                if (dq.cwiseAbs().maxCoeff() == 0.0)
                    return Vec7::Zero();
                Mat67 J = Mat67::Zero(), Jp = Mat67::Zero(), Jm = Mat67::Zero();
                se3pf::PandaArm::output(q, &J);
                const double h = 1e-6;
                se3pf::PandaArm::output(q + h * dq, &Jp);
                se3pf::PandaArm::output(q - h * dq, &Jm);
                const Vec6 Jdotdq = ((Jp - Jm) * dq) / (2.0 * h);
                const Vec6 V = J * dq;
                const Mat6 L = lambda(load);
                return scale_ * (J.transpose() * (L * Jdotdq - se3pf::ad(V).transpose() * (L * V)));
            }

            // G(q) for gravity g (libfranka's gravity_earth, default (0, 0, -9.81))
            Vec7 gravity(const Vec7& q, const RigidBody& load, const Eigen::Vector3d& g) const
            {
                Mat67 J = Mat67::Zero();
                const se3pf::Pose T = se3pf::PandaArm::output(q, &J);
                const Eigen::Matrix3d R = T.block<3, 3>(0, 0);
                Vec7 G = Vec7::Zero();
                for (const RigidBody& b : {wrist(), load}) {
                    // d p_com / dq = R (J_v - hat(c) J_w): the velocity of the centre of mass is v + w x c
                    const Eigen::Matrix<double, 3, 7> dp = R * (J.bottomRows<3>() - se3pf::hat(b.com) * J.topRows<3>());
                    G += b.mass * (dp.transpose() * (-g));
                }
                return G;  // weights are not scaled: inertiaScale models an inertia error only
            }

        private:
            Mat6 lambda(const RigidBody& load) const { return wrist().spatialInertia() + load.spatialInertia(); }
            double scale_;
    };

    // The settings of the simulated robot, from environment variables (mock/franka/README.md).
    struct Config {
        int pushJoint = -1;         // zero-based; -1: no push
        double pushTorque = 0.0, pushStart = 0.0, pushDuration = 0.0;
        double friction = 0.0;      // viscous, N m s/rad, every joint, plant only
        Vec7 coulomb = Vec7::Zero();  // Coulomb friction, N m, per joint, plant only (FRANKA_MOCK_COULOMB)
        double plantLoad = 0.0;     // extra mass at the flange, plant only, kg
        bool configuredLoad = false;  // the end effector configured in Desk is not the default Franka Hand ...
        RigidBody load;               // ... but this one (and the plant carries it too)
        Vec7 q0 = se3pf::PandaArm::restPosture();
        std::string mode = "idle";  // idle, reflex or user-stopped
        double poseError = 0.0;     // offset of the model's flange position, m
        double holdError = 1e-5;    // measured q - q_d while the robot holds the arm (Idle), every joint, rad
        double approachError = 0.0; // a joint position motion settles this far short of its goal, every joint, rad
        int lateEvery = 0;          // every lateEvery-th torque period lasts lateMs ms (0: never)
        int lateMs = 1;
        bool noRealtime = false;    // the kernel has no realtime support: control() with kEnforce throws
        int userStopAfterReads = 0; // the user stop is pressed right after the n-th readOnce (0: never)
        int nanJoint = -1;          // zero-based joint whose measured qd is NaN in one torque period ...
        double nanTime = 0.0;       // ... at this time after the torque loop starts, s
        double maxSeconds = 600.0;

        static Config fromEnvironment()
        {
            Config c;
            auto numbers = [](const char* name, int count, int alternative = -1) {
                std::array<double, 7> v{};
                const char* s = std::getenv(name);
                if (s == nullptr)
                    return std::pair<bool, std::array<double, 7>>(false, v);
                std::stringstream ss(s);
                std::string item;
                int n = 0;
                while (std::getline(ss, item, ',') && n < 7)
                    v[static_cast<std::size_t>(n++)] = std::stod(item);
                if (n != count && n != alternative)
                    throw std::invalid_argument(std::string("mock: ") + name + " needs " + std::to_string(count)
                                                + " comma-separated numbers");
                return std::pair<bool, std::array<double, 7>>(true, v);
            };
            if (auto p = numbers("FRANKA_MOCK_PUSH", 4); p.first) {
                c.pushJoint = static_cast<int>(p.second[0]) - 1;
                c.pushTorque = p.second[1];
                c.pushStart = p.second[2];
                c.pushDuration = p.second[3];
                if (c.pushJoint < 0 || c.pushJoint > 6)
                    throw std::invalid_argument("mock: FRANKA_MOCK_PUSH joint must be 1 to 7");
            }
            if (auto p = numbers("FRANKA_MOCK_FRICTION", 1); p.first)
                c.friction = p.second[0];
            if (auto p = numbers("FRANKA_MOCK_COULOMB", 7, 1); p.first) {
                // one value for every joint, or seven
                const bool perJoint = std::string(std::getenv("FRANKA_MOCK_COULOMB")).find(',') != std::string::npos;
                c.coulomb = perJoint ? fromArray(p.second) : Vec7::Constant(p.second[0]);
                if ((c.coulomb.array() < 0.0).any())
                    throw std::invalid_argument("mock: FRANKA_MOCK_COULOMB needs values >= 0");
            }
            if (auto p = numbers("FRANKA_MOCK_PLANT_LOAD", 1); p.first)
                c.plantLoad = p.second[0];
            if (auto p = numbers("FRANKA_MOCK_Q0", 7); p.first)
                c.q0 = fromArray(p.second);
            if (auto p = numbers("FRANKA_MOCK_POSE_ERROR", 1); p.first)
                c.poseError = p.second[0];
            if (auto p = numbers("FRANKA_MOCK_MAX_SECONDS", 1); p.first)
                c.maxSeconds = p.second[0];
            if (const char* m = std::getenv("FRANKA_MOCK_MODE"))
                c.mode = m;
            if (auto p = numbers("FRANKA_MOCK_CONFIGURED_LOAD", 4, 1); p.first) {
                // m, or m,cx,cy,cz: the hand's inertia, and its centre of mass unless given
                c.configuredLoad = true;
                c.load = frankaHand();
                c.load.mass = p.second[0];
                if (std::string(std::getenv("FRANKA_MOCK_CONFIGURED_LOAD")).find(',') != std::string::npos)
                    c.load.com = Eigen::Vector3d(p.second[1], p.second[2], p.second[3]);
            }
            if (auto p = numbers("FRANKA_MOCK_HOLD_ERROR", 1); p.first)
                c.holdError = p.second[0];
            if (auto p = numbers("FRANKA_MOCK_APPROACH_ERROR", 1); p.first)
                c.approachError = p.second[0];
            if (auto p = numbers("FRANKA_MOCK_LATE", 2); p.first) {
                c.lateEvery = static_cast<int>(p.second[0]);
                c.lateMs = static_cast<int>(p.second[1]);
                if (c.lateEvery < 2 || c.lateMs < 1 || c.lateMs > 10)
                    throw std::invalid_argument("mock: FRANKA_MOCK_LATE=n,ms needs n >= 2 and 1 <= ms <= 10");
            }
            if (const char* r = std::getenv("FRANKA_MOCK_NO_RT"))
                c.noRealtime = std::string(r) == "1";
            if (auto p = numbers("FRANKA_MOCK_USER_STOP_AFTER_READS", 1); p.first)
                c.userStopAfterReads = static_cast<int>(p.second[0]);
            if (auto p = numbers("FRANKA_MOCK_NAN", 2); p.first) {
                c.nanJoint = static_cast<int>(p.second[0]) - 1;
                c.nanTime = p.second[1];
                if (c.nanJoint < 0 || c.nanJoint > 6)
                    throw std::invalid_argument("mock: FRANKA_MOCK_NAN joint must be 1 to 7");
            }
            return c;
        }
    };

}  // namespace franka::mock
