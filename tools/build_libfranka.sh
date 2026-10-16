#!/usr/bin/env bash
# Builds libfranka 0.9.2 (the last libfranka for the Panda / FER; robot system >= 4.2.1) without sudo, as it was
# built for the checks of the robot driver (docs/HARDWARE.md): Poco 1.11.8 (Foundation and Net only) from source,
# then libfranka 0.9.2 against it, both Release, with the system compiler. Nothing is installed system-wide.
#   Usage:  bash tools/build_libfranka.sh [work directory, default ~/opt/franka]
#   Result: <work dir>/libfranka/build, to pass as  -DFranka_DIR=<work dir>/libfranka/build
#           (and the libfranka examples, e.g. <work dir>/libfranka/build/examples/communication_test)
# Needs git, cmake >= 3.14, g++ and Eigen 3 (libeigen3-dev). With sudo, the simpler route is
#   sudo apt install build-essential cmake git libpoco-dev libeigen3-dev
#   git clone --recursive --branch 0.9.2 https://github.com/frankaemika/libfranka && cd libfranka
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=OFF && cmake --build build -j
set -euo pipefail
W="${1:-$HOME/opt/franka}"
J="${JOBS:-$(nproc)}"
mkdir -p "$W"
cd "$W"

if [ ! -d poco-src ]; then
    git clone --depth 1 --branch poco-1.11.8-release https://github.com/pocoproject/poco poco-src
fi
cmake -S poco-src -B poco-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$W/poco-install" \
    -DBUILD_SHARED_LIBS=ON -DENABLE_FOUNDATION=ON -DENABLE_NET=ON \
    -DENABLE_ENCODINGS=OFF -DENABLE_XML=OFF -DENABLE_JSON=OFF -DENABLE_MONGODB=OFF -DENABLE_DATA_SQLITE=OFF \
    -DENABLE_REDIS=OFF -DENABLE_PDF=OFF -DENABLE_UTIL=OFF -DENABLE_SEVENZIP=OFF -DENABLE_ZIP=OFF \
    -DENABLE_CPPPARSER=OFF -DENABLE_POCODOC=OFF -DENABLE_PAGECOMPILER=OFF -DENABLE_PAGECOMPILER_FILE2PAGE=OFF \
    -DENABLE_ACTIVERECORD=OFF -DENABLE_ACTIVERECORD_COMPILER=OFF -DENABLE_TESTS=OFF -DENABLE_NETSSL=OFF \
    -DENABLE_CRYPTO=OFF -DENABLE_DATA=OFF -DENABLE_DATA_MYSQL=OFF -DENABLE_DATA_POSTGRESQL=OFF \
    -DENABLE_DATA_ODBC=OFF -DENABLE_JWT=OFF -DENABLE_APACHECONNECTOR=OFF
cmake --build poco-build -j"$J"
cmake --install poco-build

if [ ! -d libfranka ]; then
    git clone --recursive --branch 0.9.2 https://github.com/frankaemika/libfranka libfranka
fi
# The 0.9.2 tag, exactly
test "$(git -C libfranka rev-parse HEAD)" = "f3b8d775a9c847cab32684c8a316f67867761674"
cmake -S libfranka -B libfranka/build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=ON \
    -DCMAKE_PREFIX_PATH="$W/poco-install"
cmake --build libfranka/build -j"$J"
echo "Franka_DIR=$W/libfranka/build"
