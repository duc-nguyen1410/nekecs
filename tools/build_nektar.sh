#!/usr/bin/env bash
# Configure, build and install Nektar++ for nekecs (run inside WSL Ubuntu).
# Usage: bash tools/build_nektar.sh [src_dir] [prefix]
set -euo pipefail

SRC=${1:-$HOME/src/nektar}
PREFIX=${2:-$HOME/opt/nektar++}
BUILD=$SRC/build
JOBS=${JOBS:-$(nproc)}

# Only build the solvers nekecs needs. Besides saving compile time, this avoids
# a v5.10.0 bug: with NEKTAR_BUILD_SOLVER_LIBS=ON, PulseWaveSolver's CMakeLists
# lists a header (Utilities/PulseWaveSystemOutput.h) that does not exist.
# Patch solvers/CMakeLists.txt (idempotently) to honour NEKTAR_SOLVERS_ONLY.
KEEP_SOLVERS="IncNavierStokesSolver;ADRSolver"
SOLVERS_CML=$SRC/solvers/CMakeLists.txt
if ! grep -q NEKTAR_SOLVERS_ONLY "$SOLVERS_CML"; then
  sed -i 's|^    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)$|    GET_FILENAME_COMPONENT(_solver ${dir} NAME)\n    IF(NEKTAR_SOLVERS_ONLY AND NOT _solver IN_LIST NEKTAR_SOLVERS_ONLY)\n        CONTINUE()\n    ENDIF()\n&|' "$SOLVERS_CML"
  grep -q NEKTAR_SOLVERS_ONLY "$SOLVERS_CML" || { echo "patching $SOLVERS_CML failed"; exit 1; }
fi

mkdir -p "$BUILD"
cd "$BUILD"

# NEKTAR_BUILD_SOLVER_LIBS: builds the IncNavierStokes solver classes as a
# library, so nekecs can link them out-of-tree (see DESIGN.md §4).
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DNEKTAR_BUILD_SOLVERS=ON \
  -DNEKTAR_BUILD_SOLVER_LIBS=ON \
  "-DNEKTAR_SOLVERS_ONLY=$KEEP_SOLVERS" \
  -DNEKTAR_BUILD_UTILITIES=ON \
  -DNEKTAR_BUILD_DEMOS=OFF \
  -DNEKTAR_BUILD_TESTS=ON \
  -DNEKTAR_BUILD_UNIT_TESTS=OFF \
  -DNEKTAR_USE_MPI=ON \
  -DNEKTAR_USE_FFTW=ON \
  -DNEKTAR_USE_ARPACK=ON \
  -DNEKTAR_USE_SCOTCH=OFF \
  -DNEKTAR_USE_METIS=ON \
  -DTHIRDPARTY_BUILD_METIS=ON \
  -DNEKTAR_USE_HDF5=ON \
  -DNEKTAR_USE_SYSTEM_BLAS_LAPACK=ON \
  -DNEKTAR_USE_MESHGEN=OFF \
  > cmake.log 2>&1 || { grep -A6 "CMake Error" cmake.log; echo NEKTAR_CONFIGURE_FAILED; exit 1; }

if ! make -j"$JOBS" install > make.log 2>&1; then
  grep -nE 'error|Error' make.log | head -40
  echo "NEKTAR_BUILD_FAILED (see $BUILD/make.log)"
  exit 1
fi
echo "NEKTAR_BUILD_DONE prefix=$PREFIX"
