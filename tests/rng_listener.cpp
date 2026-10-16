// Reseeds the tests' shared random generator (testing::rng(), test_helpers.hpp) before every run of a test case,
// that is before each leaf section, so the random inputs of a test do not depend on which tests ran before it. A
// failure then reproduces when the test runs alone, with a name or section filter, or in random order
// (--order rand). Compiled into both test programs (tests and tests_mujoco); a listener must be registered once.
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

#include "test_helpers.hpp"

#include <cstdint>

namespace {

    class ReseedTestRng : public Catch::EventListenerBase {
        public:
            using Catch::EventListenerBase::EventListenerBase;

            void testCasePartialStarting(Catch::TestCaseInfo const&, std::uint64_t) override
            {
                testing::rng().seed(testing::kRngSeed);
            }
    };

}  // namespace

CATCH_REGISTER_LISTENER(ReseedTestRng)
