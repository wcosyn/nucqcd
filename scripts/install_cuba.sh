#!/usr/bin/env bash
# Installs Cuba 4.2.1 with the nucqcd patches (cmake/patches/cuba-4.2.1-nucqcd.patch) into a prefix
# (default ~/.local), where nucqcd's CMake finds it. Not needed if you let CMake download and build
# Cuba in the build tree (the default when no installed Cuba is found).
# Usage: scripts/install_cuba.sh [prefix] [path/to/Cuba-4.2.1.tar.gz]
set -euo pipefail
PREFIX="${1:-$HOME/.local}"
TARBALL="${2:-}"
SHA256=3789c458193440e73a603f6649acbe8f1903d61203905e1b63c00ef0db8580a8
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATCH="$HERE/../cmake/patches/cuba-4.2.1-nucqcd.patch"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
if [ -z "$TARBALL" ]; then
  TARBALL="$WORK/Cuba-4.2.1.tar.gz"
  curl -fsSL -o "$TARBALL" https://feynarts.de/cuba/Cuba-4.2.1.tar.gz
fi
echo "$SHA256  $TARBALL" | sha256sum -c -
tar xzf "$TARBALL" -C "$WORK"
cd "$WORK/Cuba-4.2.1"
patch -p1 -i "$PATCH"
# -fPIC: libcuba.a is linked into shared libraries; gnu17: the code does not compile as C23 (gcc >= 15)
CFLAGS="-O3 -fPIC -std=gnu17" ./configure --prefix="$PREFIX"
make -j1 lib
make -j1 install
echo "Cuba 4.2.1 installed in $PREFIX (lib/libcuba.a, include/cuba.h)"
