#!/usr/bin/env python3
"""Bifurcation diagram of the 2D cylinder wake from nekecs continuation runs.

Inputs: <session>_branch.dat files written by `nekecs` continuation runs:
  --eq FILES   steady-wake branch(es), with the leading eigenvalue per point
  --po FILES   shedding-orbit branch(es), with the period per point
Several files per branch (e.g. a run resumed from its last point) are
joined in the order given; repeated parameter values are dropped.

Columns (see the header line of each file):
  point Re arclength dmu/ds newton_its |F|/|x| verify |u| u_probe v_probe T St
  [indicator mode_re mode_im]

Usage: python3 tools/plot_branch.py OUT_DIR --eq cont_eq_branch.dat [...] [--po cont_po_branch.dat ...]
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(paths, ncol):
    rows, seen = [], set()
    for path in paths:
        for l in Path(path).read_text().splitlines():
            if not l.strip() or l.startswith("#"):
                continue
            r = [float(v) for v in l.split()]
            key = round(r[1], 6)
            if key not in seen:  # a resumed leg repeats its starting point
                seen.add(key)
                rows.append(r)
    return np.array(rows) if rows else np.zeros((0, ncol))


def hopf(eq):
    """Linear interpolation of the zero of Re(lambda) along the EQ branch."""
    re, sig, om = eq[:, 1], eq[:, 12], np.abs(eq[:, 14])
    for i in range(1, len(re)):
        if sig[i - 1] * sig[i] < 0:
            a = sig[i - 1] / (sig[i - 1] - sig[i])
            return re[i - 1] + a * (re[i] - re[i - 1]), om[i - 1] + a * (om[i] - om[i - 1])
    return None, None


def main(out_dir, eq_files, po_files, hopf_point=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    eq = load(eq_files, 15)
    po = load(po_files, 12)

    Rec, omc = hopf(eq) if eq.shape[1] >= 15 else (None, None)
    if hopf_point:  # refined by direct solves near the crossing
        Rec, omc = hopf_point

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
    r = np.linspace(47, 105, 100)  # fit valid for ~47 < Re < 180
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
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--eq", nargs="+", required=True)
    ap.add_argument("--po", nargs="*", default=[])
    ap.add_argument("--hopf", nargs=2, type=float, metavar=("RE_C", "OMEGA_C"),
                    help="Hopf point from a refined computation (else interpolated)")
    a = ap.parse_args()
    main(a.out_dir, a.eq, a.po, a.hopf)
