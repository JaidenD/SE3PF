#!/usr/bin/env bash
# Fetch the third-party pieces of the MuJoCo simulation (nothing here is needed for the library or the unit tests):
#   1. the MuJoCo C library (official release tarball, linux x86_64) -> <dest>/mujoco-<version>
#   2. the Franka Emika Panda model of MuJoCo Menagerie, pinned to one commit
#                                                       -> <dest>/mujoco_menagerie/franka_emika_panda
# Usage:  bash tools/fetch_mujoco.sh [dest]        (default dest: third_party/ in the project root)
# Safe to re-run: whatever is already in place and complete is kept.
# Override the pins with MUJOCO_VERSION=... or MENAGERIE_COMMIT=... in the environment.
set -euo pipefail

# --- Pins (update together with the Python package in .venv: pip install mujoco==<version>) ---
MUJOCO_VERSION="${MUJOCO_VERSION:-3.15.0}"
# sha256 of mujoco-3.15.0-linux-x86_64.tar.gz, as published with the release; checked in addition to the
# release's own .sha256 file. Only used for the default version.
MUJOCO_SHA256_PINNED="319943b332fd84c968f655f695d8f4b03cc6a68f658bacd0591f0a552744e497"
# google-deepmind/mujoco_menagerie main on 2026-10-07; franka_emika_panda last changed in 71f066a (2026-07-04).
MENAGERIE_COMMIT="${MENAGERIE_COMMIT:-0059d4335f8156206f63a35662313385f7ad6d74}"

MUJOCO_URL="https://github.com/google-deepmind/mujoco/releases/download/${MUJOCO_VERSION}"
MUJOCO_TARBALL="mujoco-${MUJOCO_VERSION}-linux-x86_64.tar.gz"
MENAGERIE_REPO="https://github.com/google-deepmind/mujoco_menagerie.git"

dest="${1:-"$(dirname "$0")/../third_party"}"
mkdir -p "$dest"
dest="$(cd "$dest" && pwd)"

# --- 1. MuJoCo C library ---
mujoco_dir="$dest/mujoco-${MUJOCO_VERSION}"
if [[ -f "$mujoco_dir/lib/libmujoco.so" && -f "$mujoco_dir/include/mujoco/mujoco.h" ]]; then
    echo "MuJoCo ${MUJOCO_VERSION}: already in $mujoco_dir"
else
    tmp="$(mktemp -d)"
    trap 'rm -rf "$tmp"' EXIT
    echo "MuJoCo ${MUJOCO_VERSION}: downloading $MUJOCO_URL/$MUJOCO_TARBALL"
    curl -fsSL --retry 3 -o "$tmp/$MUJOCO_TARBALL" "$MUJOCO_URL/$MUJOCO_TARBALL"
    actual="$(sha256sum "$tmp/$MUJOCO_TARBALL" | cut -d' ' -f1)"

    # The release publishes <tarball>.sha256 ("<hash>  <name>"); verify against it when it exists.
    if curl -fsSL --retry 3 -o "$tmp/$MUJOCO_TARBALL.sha256" "$MUJOCO_URL/$MUJOCO_TARBALL.sha256"; then
        published="$(cut -d' ' -f1 < "$tmp/$MUJOCO_TARBALL.sha256")"
        if [[ "$actual" != "$published" ]]; then
            echo "error: sha256 mismatch: got $actual, the release publishes $published" >&2
            exit 1
        fi
        echo "MuJoCo ${MUJOCO_VERSION}: sha256 matches the published $MUJOCO_TARBALL.sha256"
    else
        echo "warning: the release has no $MUJOCO_TARBALL.sha256; relying on the pinned hash only" >&2
    fi
    if [[ "$MUJOCO_VERSION" == "3.15.0" && "$actual" != "$MUJOCO_SHA256_PINNED" ]]; then
        echo "error: sha256 $actual differs from the pinned $MUJOCO_SHA256_PINNED" >&2
        exit 1
    fi

    # Extract into a staging directory first, so an interrupted run never leaves a half-extracted library.
    mkdir -p "$tmp/extract"
    tar -xzf "$tmp/$MUJOCO_TARBALL" -C "$tmp/extract"
    if [[ ! -d "$tmp/extract/mujoco-${MUJOCO_VERSION}" ]]; then
        echo "error: the tarball does not contain mujoco-${MUJOCO_VERSION}/" >&2
        exit 1
    fi
    rm -rf "$mujoco_dir"
    mv "$tmp/extract/mujoco-${MUJOCO_VERSION}" "$mujoco_dir"
    echo "MuJoCo ${MUJOCO_VERSION}: extracted to $mujoco_dir"
fi

# --- 2. Menagerie: only franka_emika_panda, shallow, at the pinned commit ---
menagerie_dir="$dest/mujoco_menagerie"
current=""
if [[ -d "$menagerie_dir/.git" ]]; then
    current="$(git -C "$menagerie_dir" rev-parse HEAD 2>/dev/null || true)"
fi
if [[ "$current" == "$MENAGERIE_COMMIT" && -f "$menagerie_dir/franka_emika_panda/panda.xml" ]]; then
    echo "Menagerie: franka_emika_panda already at ${MENAGERIE_COMMIT:0:12} in $menagerie_dir"
else
    if [[ ! -d "$menagerie_dir/.git" ]]; then
        rm -rf "$menagerie_dir"
        git init -q "$menagerie_dir"
        git -C "$menagerie_dir" remote add origin "$MENAGERIE_REPO"
    fi
    # Sparse checkout: the repository holds dozens of robots; we need one directory (cone mode also checks out
    # the few files at the repository root, such as README.md and LICENSE). init + set works with git >= 2.27.
    git -C "$menagerie_dir" sparse-checkout init --cone
    git -C "$menagerie_dir" sparse-checkout set franka_emika_panda
    echo "Menagerie: fetching franka_emika_panda at ${MENAGERIE_COMMIT:0:12}"
    git -C "$menagerie_dir" fetch -q --depth 1 --filter=blob:none origin "$MENAGERIE_COMMIT"
    git -C "$menagerie_dir" -c advice.detachedHead=false checkout -q FETCH_HEAD
    echo "Menagerie: checked out $(git -C "$menagerie_dir" rev-parse HEAD)"
fi

echo "Done. MuJoCo: $mujoco_dir"
echo "      Panda:  $menagerie_dir/franka_emika_panda"
