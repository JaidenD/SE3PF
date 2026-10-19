#pragma once
#include "se3pf/types.hpp"

namespace se3pf {

    // The Franka Emika Panda (FER) as a system of the form eq:EL (Problem formulation, Systems: robot arms). The
    // configuration q is the vector of the seven joint angles, the output Pi(q) is the pose of the flange (Franka
    // frame F) in the base frame, and J_b is defined there by (J_b(q) qdot)^ = Pi(q)^-1 T_q Pi(qdot). n = 7, so the
    // arm has one redundant joint (r = 1). Kinematics and limits are Franka's control parameters for the Panda; see
    // the design doc, appendix "Panda kinematics: modified DH (Craig)" and "Panda limits".
    //
    // The controller is a template over the system, so everything it calls is static. output and retract run in
    // every 1 kHz step: fixed-size types only, no allocation, no exceptions.
    class PandaArm {
        public:
            static constexpr int kDof = 7;
            using Config = Vec7;

            // The flange (Franka frame F) is joint 7's frame moved by Trans_z(kFlangeOffset), with its attitude.
            static constexpr double kFlangeOffset = 0.107;  // m

            // Flange pose g = Pi(q). If Jb is not null, also the body Jacobian, rows (omega, v) in the flange frame,
            // so that zeta = (g^-1 gdot)^vee = J_b(q) qdot.
            static Pose output(const Vec7& q, Mat67* Jb);

            // A curve through q with velocity dq, for the controller's central difference of A: q + h dq.
            static Vec7 retract(const Vec7& q, const Vec7& dq, double h);

            // Rest posture q_r = (0, -pi/4, 0, -3pi/4, 0, pi/2, pi/4): approach axis pointing down.
            static Vec7 restPosture();

            // Distance from q to the nearest joint position limit, negative if a joint is outside its range.
            // If joint is not null, it receives the (zero-based) index of that joint. NaN if any joint is NaN (and
            // joint is then a NaN joint), so a safety check such as margin >= 0.10 fails on a corrupted q.
            static double jointMargin(const Vec7& q, int* joint = nullptr);

            // Limits per joint. Functions rather than static Vec7 members: a static Eigen object is built at
            // program start, in an unspecified order across files, so another file's static could read it unset.
            static Vec7 qMin();     // position, rad
            static Vec7 qMax();     // position, rad
            static Vec7 dqMax();    // velocity, rad/s
            static Vec7 ddqMax();   // acceleration, rad/s^2
            static Vec7 jerkMax();  // jerk, rad/s^3
            static Vec7 tauMax();   // torque, N m

            static constexpr double tauRateMax = 1000.0;     // torque rate, N m/s, every joint
            static constexpr double flangeSpeedMax = 1.7;    // translational speed of the flange, m/s
    };

    // The margins every motion on the robot must keep (design doc, Hardware; the SafetyMonitor margins of the earlier
    // implementation). Fractions are of each joint's own limit in PandaArm. The MuJoCo verification checks them
    // (apps/panda_mujoco.cpp, tests_mujoco); the robot driver's preflight and safety monitor will use them too.
    struct HardwareBounds {
        static constexpr double jointMargin = 0.10;  // rad from every position limit, at least
        static constexpr double dqFrac = 0.5;        // |qd_i| <= 50 % of dqMax_i
        static constexpr double ddqFrac = 0.3;       // commanded |qdd*_i| <= 30 % of ddqMax_i
        static constexpr double tauFrac = 0.6;       // delivered |tau_i| (G included) <= 60 % of tauMax_i
        static constexpr double tauRateFrac = 0.6;   // commanded |d(tau_i - G_i)/dt| <= 60 % of tauRateMax
        static constexpr double flangeSpeed = 0.5;   // m/s, translational
        static constexpr double xi = 0.04;           // |xi|, in the units of the metric (Panda: m, 1 rad = 0.1 m)
        static constexpr double beta = 0.1;         // |beta|, rad
        static constexpr double sigmaMin = 1e-3;     // sigma_min(A) at least
    };

}  // namespace se3pf
