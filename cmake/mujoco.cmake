# Finds the MuJoCo C library (a release tarball unpacked by tools/fetch_mujoco.sh) and defines the imported
# target mujoco::mujoco. Use:  include(cmake/mujoco.cmake)  then  target_link_libraries(app PRIVATE mujoco::mujoco).
# Also sets MENAGERIE_PANDA_DIR, the directory of the Menagerie Panda models (see the end of this file).
#
# Where it looks: the cache variable MUJOCO_DIR if set (cmake -DMUJOCO_DIR=/path/to/mujoco-3.15.0), otherwise
# third_party/mujoco-* in the project root. A release directory holds include/mujoco/mujoco.h and lib/libmujoco.so.

set(MUJOCO_DIR "" CACHE PATH "MuJoCo release directory (contains include/mujoco/mujoco.h and lib/libmujoco.so)")

if(NOT MUJOCO_DIR)
    file(GLOB _mujoco_candidates LIST_DIRECTORIES true "${PROJECT_SOURCE_DIR}/third_party/mujoco-*")
    if(NOT _mujoco_candidates)
        message(FATAL_ERROR
            "MuJoCo not found in ${PROJECT_SOURCE_DIR}/third_party. Run  bash tools/fetch_mujoco.sh  "
            "or pass -DMUJOCO_DIR=/path/to/mujoco-<version>.")
    endif()
    # Several versions side by side: take the newest (natural order, so 3.15.0 comes after 3.9.0) and say so (set
    # MUJOCO_DIR to choose).
    list(SORT _mujoco_candidates COMPARE NATURAL)
    list(GET _mujoco_candidates -1 _mujoco_dir)
    list(LENGTH _mujoco_candidates _mujoco_count)
    if(_mujoco_count GREATER 1)
        message(STATUS "Several MuJoCo versions in third_party/; using ${_mujoco_dir} (set MUJOCO_DIR to choose)")
    endif()
else()
    set(_mujoco_dir "${MUJOCO_DIR}")
endif()

# Check the files directly instead of find_library: a cached search result would survive a change of MUJOCO_DIR.
set(MUJOCO_INCLUDE_DIR "${_mujoco_dir}/include")
set(MUJOCO_LIBRARY "${_mujoco_dir}/lib/libmujoco.so")
if(NOT EXISTS "${MUJOCO_INCLUDE_DIR}/mujoco/mujoco.h" OR NOT EXISTS "${MUJOCO_LIBRARY}")
    message(FATAL_ERROR "${_mujoco_dir} is not a MuJoCo release directory (no include/mujoco/mujoco.h or lib/libmujoco.so)")
endif()
get_filename_component(MUJOCO_LIBRARY_DIR "${MUJOCO_LIBRARY}" DIRECTORY)

if(NOT TARGET mujoco::mujoco)
    add_library(mujoco::mujoco SHARED IMPORTED)
    # SYSTEM include directory: -Wall -Wextra in our targets should not report warnings from MuJoCo's headers.
    set_target_properties(mujoco::mujoco PROPERTIES
        IMPORTED_LOCATION "${MUJOCO_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MUJOCO_INCLUDE_DIR}"
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${MUJOCO_INCLUDE_DIR}")
    # Run without LD_LIBRARY_PATH: every executable that links MuJoCo records its lib/ directory as a runtime
    # search path. (CMake adds it to the build-tree RPATH by itself; the explicit option also covers installed
    # copies and builds with CMAKE_SKIP_BUILD_RPATH.)
    set_property(TARGET mujoco::mujoco PROPERTY INTERFACE_LINK_OPTIONS "LINKER:-rpath,${MUJOCO_LIBRARY_DIR}")
endif()

message(STATUS "MuJoCo: ${_mujoco_dir}")

# The MuJoCo Menagerie Panda (panda_nohand.xml, panda.xml, scene.xml), also fetched by tools/fetch_mujoco.sh. The
# apps and tests_mujoco use it as their default model directory.
set(MENAGERIE_PANDA_DIR "${PROJECT_SOURCE_DIR}/third_party/mujoco_menagerie/franka_emika_panda" CACHE PATH
    "MuJoCo Menagerie franka_emika_panda directory (contains panda_nohand.xml)")
if(NOT EXISTS "${MENAGERIE_PANDA_DIR}/panda_nohand.xml")
    message(FATAL_ERROR
        "No panda_nohand.xml in ${MENAGERIE_PANDA_DIR}. Run  bash tools/fetch_mujoco.sh  "
        "or pass -DMENAGERIE_PANDA_DIR=/path/to/mujoco_menagerie/franka_emika_panda.")
endif()
message(STATUS "Menagerie Panda: ${MENAGERIE_PANDA_DIR}")
