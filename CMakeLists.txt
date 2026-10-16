cmake_minimum_required(VERSION 3.18)
project(se3pf LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release unless asked otherwise: a build without -O runs the tests about 100 times slower (unoptimized Eigen), and
# the closed-loop tests tagged SLOW_IN_DEBUG are hidden in Debug. Use -DCMAKE_BUILD_TYPE=Debug for Eigen's checks.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type: Debug or Release" FORCE)
endif()

# Eigen 3.4: minCoeff<PropagateNaN> (src/panda.cpp) and the 5- and 6-argument vector constructors are 3.4 only.
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# Warnings for every target of ours: -Wall -Wextra, and with -DSE3PF_WERROR=ON (tools/ci.sh) they are errors, so a
# warning fails the build however often it is repeated (a grep of the build log misses warnings of files that an
# incremental build did not recompile).
option(SE3PF_WERROR "Treat compiler warnings as errors" OFF)
add_library(se3pf_warnings INTERFACE)
target_compile_options(se3pf_warnings INTERFACE -Wall -Wextra $<$<BOOL:${SE3PF_WERROR}>:-Werror>)

# The library: one .cpp per non-template header (types.hpp, spline.hpp, redundancy.hpp and controller.hpp are
# header-only). safety, experiments and preflight are the hardware layer's library part (docs/HARDWARE.md); the
# library still includes no libfranka header: only apps/robot/robot_driver.cpp does.
add_library(se3pf
    src/se3.cpp
    src/metric.cpp
    src/path.cpp
    src/tube.cpp
    src/rigid_body.cpp
    src/panda.cpp
    src/rigid_body_example.cpp
    src/lift.cpp
    src/panda_example.cpp
    src/safety.cpp
    src/experiments.cpp
    src/preflight.cpp
)
target_include_directories(se3pf PUBLIC include)
target_link_libraries(se3pf PUBLIC Eigen3::Eigen PRIVATE se3pf_warnings)
# Lets a test forbid Eigen's heap allocations (Eigen::internal::set_is_malloc_allowed(false)) around the 1 kHz
# control step. Checked by an assert, so only in Debug builds; no effect unless a test switches it on.
target_compile_definitions(se3pf PUBLIC EIGEN_RUNTIME_NO_MALLOC)

# Apps
add_executable(rigid_body_sim apps/rigid_body_sim.cpp)
target_link_libraries(rigid_body_sim PRIVATE se3pf se3pf_warnings)
add_executable(panda_sim apps/panda_sim.cpp)
target_link_libraries(panda_sim PRIVATE se3pf se3pf_warnings)

# Tests (Catch2 v3, downloaded the first time you configure; -DFETCHCONTENT_SOURCE_DIR_CATCH2=<dir> uses a local copy)
option(SE3PF_BUILD_TESTS "Build the tests" ON)
if(SE3PF_BUILD_TESTS)
    enable_testing()
    include(FetchContent)
    FetchContent_Declare(Catch2 GIT_REPOSITORY https://github.com/catchorg/Catch2.git GIT_TAG v3.5.2 GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(Catch2)
    add_executable(tests
        tests/rng_listener.cpp
        tests/test_se3.cpp
        tests/test_metric.cpp
        tests/test_path.cpp
        tests/test_spline.cpp
        tests/test_tube.cpp
        tests/test_rigid_body.cpp
        tests/test_panda.cpp
        tests/test_controller.cpp
        tests/test_rigid_body_example.cpp
        tests/test_lift.cpp
        tests/test_panda_example.cpp
        tests/test_safety.cpp
        tests/test_preflight.cpp
        tests/test_franka_convert.cpp
    )
    target_include_directories(tests PRIVATE apps/robot)  # franka_convert.hpp: no libfranka needed
    target_link_libraries(tests PRIVATE se3pf Catch2::Catch2WithMain se3pf_warnings)
    add_test(NAME tests COMMAND tests)
endif()

# The robot driver (docs/HARDWARE.md): apps/robot/robot_driver.cpp is the only code that includes libfranka.
# robot_driver_mock is the same source against mock/franka (a simulated Panda behind the libfranka 0.9.2 API), built
# with the tests and run end to end by tests_sil. robot_driver links the real libfranka 0.9.2: the Panda (FER) needs
# 0.9.x, and libfranka 0.10 and later support only the FR3.
if(SE3PF_BUILD_TESTS)
    add_executable(robot_driver_mock apps/robot/robot_driver.cpp)
    target_include_directories(robot_driver_mock PRIVATE mock apps/robot)
    target_compile_definitions(robot_driver_mock PRIVATE SE3PF_MOCK_FRANKA)
    target_link_libraries(robot_driver_mock PRIVATE se3pf se3pf_warnings)

    # Software in the loop: runs robot_driver_mock (a separate process) for every experiment and with faults
    add_executable(tests_sil tests/rng_listener.cpp tests/sil/test_robot_driver.cpp)
    target_include_directories(tests_sil PRIVATE tests apps/robot)
    target_link_libraries(tests_sil PRIVATE se3pf Catch2::Catch2WithMain se3pf_warnings)
    target_compile_definitions(tests_sil PRIVATE SE3PF_ROBOT_DRIVER_MOCK="$<TARGET_FILE:robot_driver_mock>")
    add_dependencies(tests_sil robot_driver_mock)
    add_test(NAME tests_sil COMMAND tests_sil)
endif()

option(SE3PF_WITH_LIBFRANKA "Build robot_driver against libfranka 0.9.2 (-DFranka_DIR=<libfranka build dir>)" OFF)
option(SE3PF_ALLOW_LIBFRANKA_08 "Accept libfranka 0.8.x (a Panda on system 4.0 to 4.2.0; never linked or tested)" OFF)
if(SE3PF_WITH_LIBFRANKA)
    find_package(Franka REQUIRED)
    if(NOT DEFINED Franka_VERSION OR Franka_VERSION VERSION_LESS 0.8.0 OR NOT Franka_VERSION VERSION_LESS 0.10.0)
        message(FATAL_ERROR "robot_driver needs libfranka 0.9.x for the Panda (FER); found '${Franka_VERSION}' "
                            "(libfranka 0.10 and later support only the FR3). Build libfranka 0.9.2: docs/HARDWARE.md.")
    endif()
    if(Franka_VERSION VERSION_LESS 0.9.0)
        # Only for a Panda whose system version is older than libfranka 0.9 accepts (robot server version 4,
        # docs/HARDWARE.md). The driver compiles against the 0.8.0 headers (the API it uses is the same), but it was
        # never linked or run with 0.8: it must be asked for explicitly.
        if(NOT SE3PF_ALLOW_LIBFRANKA_08)
            message(FATAL_ERROR "libfranka ${Franka_VERSION}: the driver was linked and tested with libfranka 0.9.2 "
                                "only. For a Panda on system 4.0 to 4.2.0 (robot server version 4) configure with "
                                "-DSE3PF_ALLOW_LIBFRANKA_08=ON (only compiled against the 0.8.0 headers, never run).")
        endif()
        message(WARNING "libfranka ${Franka_VERSION} accepted by SE3PF_ALLOW_LIBFRANKA_08: never linked or tested")
    endif()
    message(STATUS "libfranka ${Franka_VERSION}: ${Franka_DIR}")
    add_executable(robot_driver apps/robot/robot_driver.cpp)
    target_include_directories(robot_driver PRIVATE apps/robot)
    target_link_libraries(robot_driver PRIVATE se3pf Franka::Franka se3pf_warnings)
endif()

# MuJoCo simulation (optional; fetch MuJoCo and the Menagerie Panda first with  bash tools/fetch_mujoco.sh).
# See docs/MUJOCO.md. The se3pf library never includes MuJoCo: everything that does is below.
option(SE3PF_WITH_MUJOCO "Build the MuJoCo simulation" OFF)
if(SE3PF_WITH_MUJOCO)
    include(cmake/mujoco.cmake)   # mujoco::mujoco and MENAGERIE_PANDA_DIR

    # Raw C API checks, without se3pf
    add_executable(mujoco_smoke apps/mujoco_smoke.cpp)
    target_link_libraries(mujoco_smoke PRIVATE mujoco::mujoco se3pf_warnings)

    # The Panda in MuJoCo (torque motors) and the closed loop of the Panda example on it
    add_library(se3pf_mujoco sim/mujoco_panda.cpp sim/mujoco_experiments.cpp)
    target_include_directories(se3pf_mujoco PUBLIC sim)
    target_link_libraries(se3pf_mujoco PUBLIC se3pf mujoco::mujoco PRIVATE se3pf_warnings)
    target_compile_definitions(se3pf_mujoco PRIVATE SE3PF_MENAGERIE_PANDA_DIR="${MENAGERIE_PANDA_DIR}")

    add_executable(panda_mujoco apps/panda_mujoco.cpp)
    target_link_libraries(panda_mujoco PRIVATE se3pf_mujoco se3pf_warnings)

    # The hardware experiments on the MuJoCo Panda (docs/HARDWARE.md): preflights with MuJoCo's model and the closed
    # loops at the torque level, the best prediction of the robot
    add_executable(experiments_mujoco apps/experiments_mujoco.cpp)
    target_include_directories(experiments_mujoco PRIVATE apps/robot)
    target_link_libraries(experiments_mujoco PRIVATE se3pf_mujoco se3pf_warnings)

    # Separate test program: the default tests never need MuJoCo
    if(SE3PF_BUILD_TESTS)
        add_executable(tests_mujoco tests/rng_listener.cpp tests/mujoco/test_mujoco_panda.cpp
                                    tests/mujoco/test_mujoco_experiments.cpp)
        target_include_directories(tests_mujoco PRIVATE tests)
        target_link_libraries(tests_mujoco PRIVATE se3pf_mujoco Catch2::Catch2WithMain se3pf_warnings)
        add_test(NAME tests_mujoco COMMAND tests_mujoco)
    endif()
endif()
