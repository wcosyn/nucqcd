#!/usr/bin/env bash
# Installs the standalone Minuit2 (the ROOT sources of Minuit2 with their own CMake project,
# https://github.com/GooFit/Minuit2, tag v6-40-02 = ROOT 6.40.02) into a prefix (default ~/.local),
# where nucqcd's CMake finds it. Not needed if you let CMake download and build Minuit2 in the build
# tree (the default when no installed Minuit2 is found, see cmake/Minuit2.cmake).
# Usage: scripts/install_minuit2.sh [prefix] [path/to/Minuit2-v6-40-02.tar.gz]
set -euo pipefail
PREFIX="${1:-$HOME/.local}"
TARBALL="${2:-}"
SHA256=72049214d0c6cdc5389ab7d8acf474638b68eaa1e66dd01389b5e7d70657fef8
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK:?}"' EXIT
if [ -z "$TARBALL" ]; then
  TARBALL="$WORK/Minuit2-v6-40-02.tar.gz"
  curl -fsSL -o "$TARBALL" https://github.com/GooFit/Minuit2/archive/refs/tags/v6-40-02.tar.gz
fi
echo "$SHA256  $TARBALL" | sha256sum -c -
tar xzf "$TARBALL" -C "$WORK"
cmake -S "$WORK/Minuit2-6-40-02" -B "$WORK/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_CXX_STANDARD=17 -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cmake --build "$WORK/build" --target Minuit2 Minuit2Math --parallel
cmake --install "$WORK/build"
echo "Minuit2 installed in $PREFIX (lib/libMinuit2.a, lib/libMinuit2Math.a, include/Minuit2)"
