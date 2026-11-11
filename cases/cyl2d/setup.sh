#!/usr/bin/env bash
# Prepare a run directory for the 2D cylinder case: the mesh and initial
# fields come from the Nektar++ regression test Cyl_AdaptiveSFD (746 quads),
# so this needs the Nektar++ source tree.
# Usage: bash cases/cyl2d/setup.sh [run_dir] [nektar_src]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
RUN=${1:-$HOME/runs/cyl2d}
NEK=${2:-$HOME/src/nektar}
T=$NEK/solvers/IncNavierStokesSolver/Tests

[ -d "$T" ] || { echo "Nektar++ source not found: $NEK (pass it as 2nd argument)"; exit 1; }
mkdir -p "$RUN"
cp "$T/Cyl_AdaptiveSFD.xml.gz" "$RUN/mesh.xml.gz"
cp "$T/Cyl_AdaptiveSFD.rst"    "$RUN/init.rst"   # near-steady SFD field
cp "$T/Cyl_AdaptiveSFD_Eig.rst" "$RUN/stab_init.rst" # Arnoldi start vector
cp "$HERE"/*.xml "$RUN"/            # session files
echo "run directory: $RUN"
echo "next: bash $HERE/run_all.sh $RUN            # EQ, stability, DNS, PO (~20 min)"
echo "  or: bash $HERE/run_all.sh $RUN --continuation  # plus the Re branches (hours)"
