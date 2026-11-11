#!/usr/bin/env bash
# Own-geometry example, start to finish (see docs/tutorial.md):
#   gmsh geometry -> Nektar++ mesh -> DNS guess -> nekecs equilibrium
#   [-> continuation in Re] -> ParaView files
#
# Usage: bash cases/cyl_gmsh/run.sh [run_dir] [--continuation]
# Environment: NEKECS_BIN (default ~/build/nekecs), NEKTAR_BIN
# (default ~/opt/nektar++/bin).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
RUN=$HOME/runs/cyl_gmsh
CONT=0
for a in "$@"; do
  case $a in
    --continuation) CONT=1 ;;
    *) RUN=$a ;;
  esac
done
BIN=${NEKECS_BIN:-$HOME/build/nekecs}
NEK=${NEKTAR_BIN:-$HOME/opt/nektar++/bin}
export OPENBLAS_NUM_THREADS=1

mkdir -p "$RUN" && cp "$HERE"/cylinder.geo "$HERE"/*.xml "$RUN"/ && cd "$RUN"

echo "[1/5] mesh: gmsh, then NekMesh"
gmsh -2 cylinder.geo -o cylinder.msh > gmsh.log
"$NEK/NekMesh" cylinder.msh mesh.xml > nekmesh.log

echo "[2/5] starting guess: DNS from a uniform stream to t = 40 (~3 min)"
"$NEK/IncNavierStokesSolver" mesh.xml dns.xml > dns.log

echo "[3/5] nekecs: equilibrium and its stability (~15 min)"
"$BIN/nekecs" mesh.xml session.xml > nekecs.log
grep -E "^\[nekecs\] (it|status|\|F\||stability)|^\[nekecs\]   [0-9]" nekecs.log

if [ $CONT = 1 ]; then
  echo "[4/5] continuation in Re, 60 -> 40 (~1-2 h)"
  "$BIN/nekecs" mesh.xml cont.xml > cont.log
  grep -E "^\[nekecs\] (point|\*\*\*)" cont.log
else
  echo "[4/5] continuation skipped (add --continuation)"
fi

echo "[5/5] ParaView files: mesh.vtu, session_eq.vtu (with vorticity)"
"$NEK/FieldConvert" -f mesh.xml mesh.vtu > /dev/null
"$NEK/FieldConvert" -f -m vorticity mesh.xml session.xml session_eq.fld session_eq.vtu > /dev/null
echo "done: results in $RUN"
