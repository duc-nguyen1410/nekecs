#!/usr/bin/env python3
"""Bifurcation diagram of the 2D cylinder wake from nekecs continuation runs.

Inputs (in RUN_DIR), written by `nekecs mesh.xml.gz cont_eq.xml / cont_po.xml`:
  cont_eq_branch.dat   steady-wake branch with leading eigenvalue per point
  cont_po_branch.dat   shedding-orbit branch (period per point)

Columns (see the header line of each file):
  point Re arclength dmu/ds newton_its |F|/|x| verify |u| u_probe v_probe T St
  [indicator mode_re mode_im]

Usage: python3 tools/plot_branch.py RUN_DIR OUT_DIR
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    rows = [l.split() for l in Path(path).read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    return np.array(rows, dtype=float) if rows else np.zeros((0, 15))


def hopf(eq):
    """Linear interpolation of the zero of Re(lambda) along the EQ branch."""
    re, sig, om = eq[:, 1], eq[:, 12], np.abs(eq[:, 14])
    for i in range(1, len(re)):
        if sig[i - 1] * sig[i] < 0:
            a = sig[i - 1] / (sig[i - 1] - sig[i])
            return re[i - 1] + a * (re[i] - re[i - 1]), om[i - 1] + a * (om[i] - om[i - 1])
    return None, None


def main(run_dir, out_dir):
    run_dir, out_dir = Path(run_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    eq = load(run_dir / "cont_eq_branch.dat")
    po_path = run_dir / "cont_po_branch.dat"
    po = load(po_path) if po_path.exists() else np.zeros((0, 12))

    Rec, omc = hopf(eq) if eq.shape[1] >= 15 else (None, None)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)

    ax = axes[0]
    ax.axhline(0, color="0.6", lw=0.8)
    ax.plot(eq[:, 1], eq[:, 12], "o-", color="C0", ms=4, label="steady wake: leading Re($\\lambda$)")
    if Rec is not None:
        ax.axvline(Rec, color="C3", ls="--", lw=1)
        ax.annotate(f"Hopf  Re$_c$ = {Rec:.2f}", (Rec, 0), xytext=(8, 12),
                    textcoords="offset points", color="C3", fontsize=9)
    ax.set_xlabel("Re")
    ax.set_ylabel("growth rate  Re($\\lambda$) D / U")
    ax.set_title("Steady wake: leading eigenvalue along the branch", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")

    ax = axes[1]
    ax.plot(eq[:, 1], np.abs(eq[:, 14]) / (2 * np.pi), "o-", color="C0", ms=4,
            label="linear mode of the steady wake, $\\omega / 2\\pi$")
    if len(po):
        ax.plot(po[:, 1], po[:, 11], "s-", color="C1", ms=4,
                label="shedding periodic orbit, 1 / T")
    r = np.linspace(50, 160, 100)
    ax.plot(r, -3.3265 / r + 0.1816 + 1.6e-4 * r, ":", color="0.4", lw=1,
            label="Williamson (1998) fit to experiments")
    if Rec is not None:
        ax.plot([Rec], [omc / (2 * np.pi)], "*", color="C3", ms=11,
                label=f"Hopf point, St$_c$ = {omc / (2 * np.pi):.4f}")
    ax.set_xlabel("Re")
    ax.set_ylabel("Strouhal number  f D / U")
    ax.set_title("Frequencies along the branches", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")

    out = out_dir / "cyl2d_bifurcation.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)
    if Rec is not None:
        print(f"Hopf: Re_c = {Rec:.4f}, omega_c = {omc:.5f}, St_c = {omc / (2 * np.pi):.5f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
