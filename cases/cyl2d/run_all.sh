#!/usr/bin/env bash
# Run the whole cyl2d pipeline in order, skipping steps whose output exists.
#
#   1. equilibrium + its stability   session.xml   -> session_eq.fld
#   2. Nektar++ Arnoldi cross-check  stability.xml -> mesh_eig_0.fld (figures)
#   3. DNS to the shedding cycle     shed.xml      -> shed_shed.fld
#   4. shedding periodic orbit       po.xml        -> po_po.fld
#   5. (--continuation) EQ and PO branches in Re, run in parallel
#                                    cont_eq.xml, cont_po.xml -> *_branch.dat
#
# Steps 1-4 take ~20 min on one core; step 5 several hours.
#
# Usage: bash cases/cyl2d/run_all.sh [run_dir] [--continuation]
# Environment: NEKECS_BIN (default ~/build/nekecs), NEKTAR_BIN
# (default ~/opt/nektar++/bin).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
RUN=$HOME/runs/cyl2d
CONT=0
for a in "$@"; do
  case $a in
    --continuation) CONT=1 ;;
    *) RUN=$a ;;
  esac
done
BIN=${NEKECS_BIN:-$HOME/build/nekecs}
NEK=${NEKTAR_BIN:-$HOME/opt/nektar++/bin}
export OPENBLAS_NUM_THREADS=1   # threaded BLAS is ~8x slower here

[ -f "$RUN/mesh.xml.gz" ] || bash "$HERE/setup.sh" "$RUN"
cd "$RUN"

step() { # step <output> <log> <command...>
  local out=$1 log=$2; shift 2
  if [ -f "$out" ]; then echo "[run_all] $out exists, skipping"; return; fi
  echo "[run_all] $* > $log"
  "$@" > "$log" 2>&1 || { echo "[run_all] FAILED, see $RUN/$log"; tail -20 "$log"; exit 1; }
  grep -E "^\[(nekecs|shed)\] (status|T =|stability|wrote)" "$log" || true
}

step session_eq.fld session.log "$BIN/nekecs" mesh.xml.gz session.xml
step mesh_eig_0.fld stability.log "$NEK/IncNavierStokesSolver" mesh.xml.gz stability.xml
step shed_shed.fld shed.log "$BIN/shed" mesh.xml.gz shed.xml
step po_po.fld po.log "$BIN/nekecs" mesh.xml.gz po.xml

if [ $CONT = 1 ]; then
  echo "[run_all] continuation in Re (EQ and PO branches in parallel; hours)"
  "$BIN/nekecs" mesh.xml.gz cont_eq.xml > cont_eq.log 2>&1 &
  "$BIN/nekecs" mesh.xml.gz cont_po.xml > cont_po.log 2>&1 &
  wait
  grep -h "\*\*\*\|continuation " cont_eq.log cont_po.log || true
fi
echo "[run_all] done: results in $RUN"
