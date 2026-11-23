#!/usr/bin/env bash
# Install PARTONS v5 and its dependencies so nucqcd can be built with -DPARTONS_ENABLED=yes.
#
#   Usage: scripts/install_partons.sh
#
# What it does (everything except the apt step installs without sudo into $PREFIX):
#   1. apt-installs the system dependencies that are missing (asks for sudo only if needed)
#   2. builds SFML 2.6 (system module only) -- Ubuntu >= 26.04 ships SFML 3, which removed
#      sf::Thread / sf::Mutex / sf::Lock that PARTONS uses
#   3. builds LHAPDF with its Python bindings (module `lhapdf`, installed in the user site-packages of
#      python3), patching lhapdf-config so it accepts --incdir, which PARTONS' FindLHAPDF calls, and
#      installs the PDF sets CT18LO and NNPDFpol11_100 (used by nucqcd's LHAPDF tests); nucqcd itself
#      uses LHAPDF when it is found, also without PARTONS ("LHAPDF:<set>" structure functions)
#   4. builds APFEL++
#   5. clones/updates and builds elementary-utils, numa and partons (in that order)
#   6. writes partons.properties / logger.properties into GPD_Partons/bin (where the GPD programs
#      are installed; PARTONS needs them next to the executable) unless they already exist
#   7. test-configures nucqcd with PARTONS enabled and builds libGPD
#
# Everything is linked with an RPATH to $PREFIX/lib, so no LD_LIBRARY_PATH is needed.
# Re-running is safe: existing clones are updated and rebuilt.
#
# Override defaults with environment variables, e.g.  PREFIX=/opt/partons JOBS=4 scripts/install_partons.sh
set -euo pipefail

PREFIX=${PREFIX:-$HOME/.local}          # install prefix (nucqcd's CMake searches ~/.local by default)
SRC=${SRC:-$HOME/Git}                   # where the apfelxx and PARTONS repos are cloned
JOBS=${JOBS:-$(nproc)}
SFML_VERSION=${SFML_VERSION:-2.6.2}
LHAPDF_VERSION=${LHAPDF_VERSION:-6.5.6}

REPO_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
WORK=$(mktemp -d)                       # temporary sources for SFML/LHAPDF, removed on exit
trap 'rm -rf "$WORK"' EXIT

CMAKE_COMMON=(-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX"
              -DCMAKE_INSTALL_RPATH="$PREFIX/lib" -DCMAKE_PREFIX_PATH="$PREFIX")

step() { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }

# build_cmake <source dir> [extra cmake args...]: configure, build and install a CMake project
build_cmake() {
  local src=$1; shift
  cmake -S "$src" -B "$src/build" "${CMAKE_COMMON[@]}" "$@" > "$src/build.log" 2>&1 \
    || { echo "cmake configure failed, see $src/build.log"; exit 1; }
  cmake --build "$src/build" -j "$JOBS" >> "$src/build.log" 2>&1 \
    || { echo "build failed, see $src/build.log"; exit 1; }
  cmake --install "$src/build" >> "$src/build.log" 2>&1
}

# clone_or_update <url> <dir>
clone_or_update() {
  if [[ -d $2/.git ]]; then git -C "$2" pull --ff-only -q; else git clone -q "$1" "$2"; fi
}

mkdir -p "$PREFIX" "$SRC"

# ---- 1. system packages --------------------------------------------------------------------
step "System packages"
PKGS=(build-essential cmake git curl libgsl-dev libxml2-dev libeigen3-dev libcln-dev python3-dev cython3)
MISSING=()
for p in "${PKGS[@]}"; do dpkg -s "$p" > /dev/null 2>&1 || MISSING+=("$p"); done
if (( ${#MISSING[@]} )); then
  echo "Installing: ${MISSING[*]}"
  sudo apt-get install -y "${MISSING[@]}"
else
  echo "all present"
fi

# ---- 2. SFML 2.6 (system module only) -------------------------------------------------------
step "SFML $SFML_VERSION"
git clone -q --depth 1 -b "$SFML_VERSION" -c advice.detachedHead=false https://github.com/SFML/SFML.git "$WORK/sfml"
build_cmake "$WORK/sfml" -DSFML_BUILD_WINDOW=OFF -DSFML_BUILD_GRAPHICS=OFF \
                         -DSFML_BUILD_AUDIO=OFF -DSFML_BUILD_NETWORK=OFF

# ---- 3. LHAPDF (with Python bindings) ---------------------------------------------------------
# - cython3 is needed: the lhapdf.cpp shipped in the tarball is too old for recent Pythons, and
#   configure only regenerates it when it finds Cython.
# - LHAPDF's default module path (<prefix>/lib/python3.X/dist-packages) is not on sys.path on
#   Debian/Ubuntu, so install into the user site-packages of python3 instead (not seen by venvs
#   created without --system-site-packages).
# - the RPATH lets the module find libLHAPDF in $PREFIX/lib without LD_LIBRARY_PATH.
step "LHAPDF $LHAPDF_VERSION"
curl -sSfL "https://lhapdf.hepforge.org/downloads/LHAPDF-$LHAPDF_VERSION.tar.gz" | tar xz -C "$WORK"
PYSITE=$(python3 -c 'import site; print(site.getusersitepackages())')
# (steps chained with &&: set -e does not apply inside a subshell followed by ||;
#  the log is printed on failure because $WORK is removed on exit)
(
  cd "$WORK/LHAPDF-$LHAPDF_VERSION" &&
  ./configure --prefix="$PREFIX" PYTHON=python3 LHAPDF_PYTHONPATH="$PYSITE" \
              LDFLAGS="-Wl,-rpath,$PREFIX/lib" > build.log 2>&1 &&
  make -j "$JOBS" >> build.log 2>&1 &&
  make install >> build.log 2>&1
) || { tail -30 "$WORK/LHAPDF-$LHAPDF_VERSION/build.log"; echo "LHAPDF build failed"; exit 1; }
python3 -c 'import lhapdf; print("lhapdf python module", lhapdf.__version__, "in", lhapdf.__file__)' \
  || { echo "LHAPDF python module not importable"; exit 1; }
# lhapdf-config advertises --incdir but only accepts --includedir; PARTONS calls --incdir
sed -i 's/^\(\s*\)--includedir)\(\s*\)show_incdir=1 ;;/\1--includedir|--incdir)\2show_incdir=1 ;;/' "$PREFIX/bin/lhapdf-config"
"$PREFIX/bin/lhapdf-config" --incdir > /dev/null || { echo "lhapdf-config --incdir patch failed"; exit 1; }
# PDF sets for nucqcd's tests (tests/unit/unit_lhapdf.cpp), into $PREFIX/share/LHAPDF
for set in CT18LO NNPDFpol11_100; do
  [ -d "$PREFIX/share/LHAPDF/$set" ] || "$PREFIX/bin/lhapdf" install "$set" \
    || { echo "lhapdf install $set failed"; exit 1; }
done

# ---- 4. APFEL++ -----------------------------------------------------------------------------
step "APFEL++"
clone_or_update https://github.com/vbertone/apfelxx.git "$SRC/apfelxx"
build_cmake "$SRC/apfelxx"

# ---- 5. PARTONS: elementary-utils -> numa -> partons ----------------------------------------
# lhapdf-config / apfelxx-config are found on PATH by the PARTONS find modules
export PATH="$PREFIX/bin:$PATH"
for r in elementary-utils numa partons; do
  step "PARTONS: $r"
  clone_or_update "https://github.com/3d-partons/$r.git" "$SRC/$r"
  build_cmake "$SRC/$r" -DSFML_ROOT="$PREFIX"
done

# ---- 6. PARTONS runtime configuration for nucqcd -----------------------------------------
# PARTONS refuses to start without partons.properties in the directory of the executable.
# The GPD programs are installed into GPD_Partons/bin (git-ignored), so put the config there.
# Existing files are left alone, so local edits survive re-running this script.
step "PARTONS configuration in GPD_Partons/bin"
BIN="$REPO_DIR/GPD_Partons/bin"
mkdir -p "$BIN/tmp"
if [[ ! -f $BIN/xmlSchema.xsd ]]; then
  # XML schema (only used by PARTONS' XML automation), taken from the partons-example templates
  clone_or_update https://github.com/3d-partons/partons-example.git "$SRC/partons-example"
  cp "$SRC/partons-example/data/xmlSchema.xsd" "$BIN/"
fi
if [[ ! -f $BIN/partons.properties ]]; then
  cat > "$BIN/partons.properties" <<EOF
# PARTONS v5 configuration (must sit next to the executable)
log.file.path = $BIN/logger.properties
xml.schema.file.path = $BIN/xmlSchema.xsd
computation.nb.processor = 1
gpd.service.batch.size = 1000
collinear_distribution.service.batch.size = 1000
ccf.service.batch.size = 1000
observable.service.batch.size = 1000
EOF
  echo "wrote $BIN/partons.properties"
else
  echo "$BIN/partons.properties exists, left unchanged"
fi
if [[ ! -f $BIN/logger.properties ]]; then
  cat > "$BIN/logger.properties" <<EOF
enable = true
default.level = INFO
print.mode = COUT
log.folder.path = $BIN/tmp
EOF
  echo "wrote $BIN/logger.properties"
else
  echo "$BIN/logger.properties exists, left unchanged"
fi

# ---- 7. test nucqcd ---------------------------------------------------------------------
step "Test: configure nucqcd with PARTONS and build libGPD"
if command -v root-config > /dev/null || [[ -x /opt/root/bin/root-config ]]; then
  cmake -S "$REPO_DIR" -B "$WORK/nucqcd" -DPARTONS_ENABLED=yes > "$WORK/test.log" 2>&1 \
    && cmake --build "$WORK/nucqcd" -j "$JOBS" --target GPD >> "$WORK/test.log" 2>&1 \
    && echo "libGPD builds fine" \
    || { tail -30 "$WORK/test.log"; echo "nucqcd test build FAILED"; exit 1; }
else
  echo "ROOT not found (no root-config on PATH or in /opt/root), skipping the nucqcd test"
fi

step "Done. PARTONS and its dependencies are installed in $PREFIX"
