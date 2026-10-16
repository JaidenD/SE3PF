// The paper's Panda example, built once for the whole test program (about 4 s in Debug, 0.03 s in Release) and
// shared by every test file that needs it: test_lift, test_controller, test_panda_example and tests_mujoco.
#pragma once
#include "se3pf/panda_example.hpp"

namespace testing {

    // An inline function, so its static is one object across all translation units. Built on first use, so a test
    // program that never asks for it does not build it.
    inline const se3pf::PandaExample& pandaExample()
    {
        static const se3pf::PandaExample ex;
        return ex;
    }

}  // namespace testing
