// The tubes of the paper's two examples, built once for the whole test program and shared by test_tube and
// test_controller. Built here from their parts (not taken from PandaExample or RigidBodyExample), so that the tube
// tests do not depend on the lift or the rigid body. Function-local statics, so each path outlives its tube.
#pragma once
#include "se3pf/metric.hpp"
#include "se3pf/path.hpp"
#include "se3pf/tube.hpp"
#include "test_helpers.hpp"

namespace testing {

    // Panda path (coning circle, p0 = (0.40, 0, 0.45), r = 0.12, alpha = 0.25), l = 0.1 m: W = diag(0.01 I, I)
    inline const se3pf::Tube& pandaTube()
    {
        static const se3pf::Metric W = se3pf::Metric::product(0.01, 1.0);
        static const se3pf::ConingCircle path(restAttitude(), Vector3d(0.40, 0.0, 0.45), 0.12, 0.25, W);
        static const se3pf::Tube tube(path, W);
        return tube;
    }

    // The same path with l = 1 m: W = I
    inline const se3pf::Tube& pandaTube1()
    {
        static const se3pf::Metric W = se3pf::Metric::product(1.0, 1.0);
        static const se3pf::ConingCircle path(restAttitude(), Vector3d(0.40, 0.0, 0.45), 0.12, 0.25, W);
        static const se3pf::Tube tube(path, W);
        return tube;
    }

    // Rigid-body path (tumbling circle, r = 1 m), W = I
    inline const se3pf::Tube& rigidTube()
    {
        static const se3pf::Metric W = se3pf::Metric::product(1.0, 1.0);
        static const se3pf::TumblingCircle path(1.0, W);
        static const se3pf::Tube tube(path, W);
        return tube;
    }

}  // namespace testing
