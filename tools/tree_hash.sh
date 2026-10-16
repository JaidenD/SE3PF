#!/usr/bin/env bash
# Prints a short hash of the project's sources (code, build files, tests and tools; not the docs, build directories,
# downloads or logs), so that the tree run on the lab PC can be checked against the tree that tools/ci.sh validated:
# ci.sh prints it as its last line ("CI PASSED for source tree <hash>"); on the lab PC, `bash tools/tree_hash.sh` must
# print the same hash (docs/HARDWARE.md, section 1). Independent of where the tree is.
#   Usage:  bash tools/tree_hash.sh
set -euo pipefail
cd "$(dirname "$0")/.."
find CMakeLists.txt cmake include src apps sim mock tests tools -type f \
    \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' -o -name 'CMakeLists.txt' -o -name '*.cmake' -o -name '*.sh' \
       -o -name '*.py' \) -not -path '*/__pycache__/*' -print0 \
    | LC_ALL=C sort -z | xargs -0 sha256sum | sha256sum | cut -c1-16
