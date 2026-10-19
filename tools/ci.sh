#!/usr/bin/env bash
# Local CI: the checks a pull request must pass before review.
#   1. Debug and Release builds with -DSE3PF_WERROR=ON: any compiler warning in our code is an error (and the build
#      log is also grepped for warnings, ignoring third-party code in _deps).
#   2. The whole test suite passes in both builds (in Release including the tests tagged [slow], which Debug hides).
#   3. Release: the apps rigid_body_sim and panda_sim run and exit 0 (their exit code checks the runs are valid).
#   4. Release: software in the loop (docs/HARDWARE.md): tests_sil runs robot_driver_mock, the robot driver built
#      against mock/franka, end to end for every experiment and with injected faults. (Debug builds both, but does
#      not run them: the simulated runs take minutes there.)
#   5. Only if MuJoCo was fetched (third_party/mujoco-*, see tools/fetch_mujoco.sh): a Release build with
#      -DSE3PF_WITH_MUJOCO=ON, mujoco_smoke passes, tests_mujoco passes, every check of panda_mujoco's summary
#      passes (docs/MUJOCO.md), and experiments_mujoco runs every hardware experiment at the torque level (all must
#      end at rest within every margin). Skipped otherwise.
#   6. The real libfranka 0.9.2, if given: LIBFRANKA_INCLUDE=<libfranka source>/include compiles the driver against
#      its headers (-fsyntax-only); FRANKA_DIR=<libfranka build dir> builds robot_driver with
#      -DSE3PF_WITH_LIBFRANKA=ON, links it, runs 'robot_driver list', and checks that the test-only options (--yes,
#      --inject, --compensation-scale) do not exist in it and that it refuses --no-rt for a run (no robot needed).
#      Skipped otherwise.
#   7. Prints the hash of the validated source tree (tools/tree_hash.sh): run the robot session from a tree with the
#      same hash (docs/HARDWARE.md, section 1).
# Run from anywhere:  bash tools/ci.sh
# Without network access, point CMake at a local Catch2 v3.5.2 source tree:  CATCH2_SRC=/path/to/Catch2 bash tools/ci.sh
# With libfranka:  FRANKA_DIR=~/libfranka/build LIBFRANKA_INCLUDE=~/libfranka/include bash tools/ci.sh
set -euo pipefail
cd "$(dirname "$0")/.."

extra=()
if [[ -n "${CATCH2_SRC:-}" ]]; then
    extra+=("-DFETCHCONTENT_SOURCE_DIR_CATCH2=${CATCH2_SRC}")
fi

# configure, build and check for warnings: stage <build dir> <build type> [extra cmake arguments...]
stage() {
    local dir="$1" cfg="$2"
    shift 2
    echo "== ${dir}: configure (${cfg})"
    cmake -S . -B "$dir" -DCMAKE_BUILD_TYPE="$cfg" -DSE3PF_WERROR=ON ${extra[@]+"${extra[@]}"} "$@" \
        > "$dir.configure.log" 2>&1 \
        || { cat "$dir.configure.log"; echo "CI FAILED: configure (${dir})"; exit 1; }

    echo "== ${dir}: build"
    if ! cmake --build "$dir" -j > "$dir.build.log" 2>&1; then
        grep -E "error|warning" "$dir.build.log" | head -40 || tail -40 "$dir.build.log"
        echo "CI FAILED: build (${dir})"
        exit 1
    fi
    if grep -E "warning:" "$dir.build.log" | grep -v "/_deps/"; then
        echo "CI FAILED: compiler warnings (${dir}); the definition of done allows none"
        exit 1
    fi
}

# run an app and check its exit code: app <build dir> <name> [arguments...]
app() {
    local dir="$1" name="$2"
    shift 2
    echo "== ${dir}: ${name} (output in ${dir}.${name}.log)"
    "./$dir/$name" "$@" > "$dir.$name.log" 2>&1 \
        || { grep -E "FAIL|error" "$dir.$name.log" || tail -20 "$dir.$name.log"; echo "CI FAILED: ${name}"; exit 1; }
}

for cfg in Debug Release; do
    dir="build-ci-${cfg,,}"
    stage "$dir" "$cfg"
    echo "== ${dir}: tests"
    "./$dir/tests" || { echo "CI FAILED: tests (${cfg})"; exit 1; }
done

# The apps write their logs into the build directory, so out/ is not touched.
app build-ci-release rigid_body_sim build-ci-release/out
app build-ci-release panda_sim build-ci-release/out

# Software in the loop: the robot driver against the mock libfranka, every experiment end to end and with faults
echo "== build-ci-release: tests_sil (robot_driver_mock)"
./build-ci-release/tests_sil || { echo "CI FAILED: tests_sil"; exit 1; }
app build-ci-release robot_driver_mock list

if compgen -G "third_party/mujoco-*" > /dev/null; then
    dir="build-ci-mujoco"
    stage "$dir" Release -DSE3PF_WITH_MUJOCO=ON
    app "$dir" mujoco_smoke third_party/mujoco_menagerie/franka_emika_panda/panda_nohand.xml
    echo "== ${dir}: tests_mujoco"
    "./$dir/tests_mujoco" || { echo "CI FAILED: tests_mujoco"; exit 1; }
    # The full runs (20 s each, about 15 s in all) and every check of their summary
    app "$dir" panda_mujoco "$dir/out"
    # The hardware experiments at the torque level (about 25 s): every one must end at rest within every margin
    app "$dir" experiments_mujoco "$dir/out"
else
    echo "== MuJoCo: skipped (no third_party/mujoco-*; bash tools/fetch_mujoco.sh adds it)"
fi

# The real libfranka 0.9.2 (the driver is the only code that includes it)
if [[ -n "${LIBFRANKA_INCLUDE:-}" ]]; then
    echo "== robot_driver.cpp against the libfranka headers in ${LIBFRANKA_INCLUDE}"
    eigen_inc="$(pkg-config --cflags-only-I eigen3 2>/dev/null || echo -I/usr/include/eigen3)"
    # shellcheck disable=SC2086
    g++ -std=c++17 -fsyntax-only -Wall -Wextra -Werror -I "${LIBFRANKA_INCLUDE}" -I include -I apps/robot ${eigen_inc} \
        apps/robot/robot_driver.cpp || { echo "CI FAILED: robot_driver.cpp against the real libfranka headers"; exit 1; }
else
    echo "== libfranka headers: skipped (set LIBFRANKA_INCLUDE=<libfranka 0.9.2>/include)"
fi
if [[ -n "${FRANKA_DIR:-}" ]]; then
    dir="build-ci-franka"
    stage "$dir" Release -DSE3PF_WITH_LIBFRANKA=ON -DFranka_DIR="${FRANKA_DIR}"
    app "$dir" robot_driver list
    # The real driver: no test-only options, and no motion without realtime scheduling (exit 1 before connecting)
    for args in "0.0.0.0 hold --yes" "0.0.0.0 hold --inject=tube@1" "0.0.0.0 hold --compensation-scale=0" \
                "0.0.0.0 hold --no-rt"; do
        code=0
        # shellcheck disable=SC2086
        "./$dir/robot_driver" $args > /dev/null 2>&1 || code=$?
        if [[ "$code" != 1 ]]; then
            echo "CI FAILED: 'robot_driver ${args}' exited ${code}, not 1 (usage)"
            exit 1
        fi
    done
    if strings "./$dir/robot_driver" | grep -q -e "--inject" -e "(--yes)" -e "--compensation-scale"; then
        echo "CI FAILED: robot_driver contains the mock-only options"
        exit 1
    fi
else
    echo "== libfranka: skipped (set FRANKA_DIR=<libfranka 0.9.2 build directory> to build and link robot_driver)"
fi

echo "CI PASSED for source tree $(bash tools/tree_hash.sh)"
