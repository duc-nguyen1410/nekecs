#!/usr/bin/env bash
# Make the README figures for the cyl2d case.
# Needs session_eq.fld (from nekecs) and mesh_eig_0.fld (from stability.xml)
# in RUN_DIR. Usage: bash tools/make_cyl2d_figures.sh [run_dir] [out_dir]
# Environment: NEKTAR_BIN (FieldConvert location, default ~/opt/nektar++/bin).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
RUN=${1:-$HOME/runs/cyl2d}
OUT=${2:-$HERE/../docs/figures}
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)   # absolute: we cd into RUN below
export PATH=${NEKTAR_BIN:-$HOME/opt/nektar++/bin}:$PATH OPENBLAS_NUM_THREADS=1
cd "$RUN"

# Regular grid x in [-2, 14], y in [-3.5, 3.5].
python3 - <<'PY'
import numpy as np
X, Y = np.meshgrid(np.linspace(-2, 14, 801), np.linspace(-3.5, 3.5, 351))
with open("grid.pts", "w") as f:
    f.write('<?xml version="1.0" encoding="utf-8" ?>\n<NEKTAR>\n<POINTS DIM="2" FIELDS="">\n')
    for a, b in zip(X.ravel(), Y.ravel()):
        f.write(f"{a:.6f} {b:.6f}\n")
    f.write("</POINTS>\n</NEKTAR>\n")
PY

pairs=("session_eq.fld eq_vort" "mesh_eig_0.fld eig0_vort")
[ -f po_po.fld ] && pairs+=("po_po.fld po_vort")   # periodic orbit, if computed
for pair in "${pairs[@]}"; do
  set -- $pair
  FieldConvert -f -m vorticity mesh.xml.gz session.xml "$1" "$2.fld" > /dev/null
  FieldConvert -f -m interppoints:fromxml=mesh.xml.gz,session.xml:fromfld="$2.fld":topts=grid.pts \
    "$2_grid.pts" > /dev/null
done

python3 "$HERE/plot_cyl2d.py" "$RUN" "$OUT"
