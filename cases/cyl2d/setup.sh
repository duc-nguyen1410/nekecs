#!/usr/bin/env bash
# Prepare a run directory for the 2D cylinder case.
# Usage: bash cases/cyl2d/setup.sh [run_dir] [nektar_src]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
RUN=${1:-$HOME/runs/cyl2d}
NEK=${2:-$HOME/src/nektar}
T=$NEK/solvers/IncNavierStokesSolver/Tests

mkdir -p "$RUN"
cp "$T/Cyl_AdaptiveSFD.xml.gz" "$RUN/mesh.xml.gz"
cp "$T/Cyl_AdaptiveSFD.rst"    "$RUN/init.rst"   # near-steady SFD field
cp "$T/Cyl_AdaptiveSFD_Eig.rst" "$RUN/stab_init.rst" # Arnoldi start vector
cp "$HERE/session.xml"         "$RUN/session.xml"
cp "$HERE/stability.xml"       "$RUN/stability.xml"
cp "$HERE/shed.xml"            "$RUN/shed.xml"
cp "$HERE/po.xml"              "$RUN/po.xml"
cp "$HERE/cont_eq.xml"         "$RUN/cont_eq.xml"
cp "$HERE/cont_po.xml"         "$RUN/cont_po.xml"
echo "run directory: $RUN"
echo "  OPENBLAS_NUM_THREADS=1 nekecs mesh.xml.gz session.xml"
echo "  OPENBLAS_NUM_THREADS=1 IncNavierStokesSolver mesh.xml.gz stability.xml"
echo "  OPENBLAS_NUM_THREADS=1 shed mesh.xml.gz shed.xml     # DNS -> shed_shed.fld"
echo "  OPENBLAS_NUM_THREADS=1 nekecs mesh.xml.gz po.xml     # periodic orbit"
echo "  OPENBLAS_NUM_THREADS=1 nekecs mesh.xml.gz cont_eq.xml # EQ branch in Re (needs session_eq.fld)"
echo "  OPENBLAS_NUM_THREADS=1 nekecs mesh.xml.gz cont_po.xml # PO branch in Re (needs po_po.fld)"
