#!/usr/bin/env python3
"""Plot the cyl2d equilibrium and its leading eigenmode.

Inputs (in RUN_DIR), produced by cases/cyl2d (see README.md):
  eq_vort_grid.pts    equilibrium u, v, p, W_z on a regular grid
  eig0_vort_grid.pts  leading linear eigenmode (real part), same grid

Usage: python3 tools/plot_cyl2d.py RUN_DIR OUT_DIR
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

R = 0.5  # cylinder radius (diameter 1)


def load_grid(path):
    """Read a Nektar++ .pts file written on a tensor grid; return x, y, fields."""
    text = Path(path).read_text()
    names = re.search(r'FIELDS="([^"]*)"', text).group(1).split(",")
    rows = [l.split() for l in text.splitlines() if re.match(r"^\s*-?[0-9]", l)]
    data = np.array(rows, dtype=float)
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    # points were written with x fastest (meshgrid order)
    fields = {n: data[:, 2 + i].reshape(len(y), len(x)) for i, n in enumerate(names)}
    X, Y = np.meshgrid(x, y)
    inside = X**2 + Y**2 < R**2
    for f in fields.values():
        f[inside] = np.nan
    return x, y, fields


def add_cylinder(ax):
    ax.add_patch(Circle((0, 0), R, facecolor="0.55", edgecolor="k", lw=0.8, zorder=5))
    ax.set_aspect("equal")
    ax.set_xlabel("x / D")
    ax.set_ylabel("y / D")


def main(run_dir, out_dir):
    run_dir, out_dir = Path(run_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- equilibrium -------------------------------------------------------
    x, y, f = load_grid(run_dir / "eq_vort_grid.pts")
    u, v, w = f["u"], f["v"], f["W_z"]

    fig, axes = plt.subplots(2, 1, figsize=(10, 6.6), constrained_layout=True)

    ax = axes[0]
    cf = ax.contourf(x, y, u, levels=np.linspace(-0.25, 1.25, 31), cmap="RdBu_r", extend="both")
    ax.contour(x, y, u, levels=[0.0], colors="k", linewidths=1.2)
    ax.streamplot(x, y, np.nan_to_num(u), np.nan_to_num(v), density=(2.2, 1.0),
                  color="k", linewidth=0.4, arrowsize=0.6)
    add_cylinder(ax)
    ax.set_xlim(x[0], x[-1]); ax.set_ylim(y[0], y[-1])
    ax.set_title("Steady (unstable) equilibrium, Re = 100: streamwise velocity u, "
                 "streamlines, u = 0 contour (black)", fontsize=10)
    fig.colorbar(cf, ax=ax, label="u / U", shrink=0.9)

    ax = axes[1]
    lev = np.linspace(-3, 3, 41)
    cf = ax.contourf(x, y, np.clip(w, -3, 3), levels=lev, cmap="RdBu_r")
    add_cylinder(ax)
    ax.set_xlim(x[0], x[-1]); ax.set_ylim(y[0], y[-1])
    ax.set_title("Spanwise vorticity $\\omega_z$ (clipped at ±3)", fontsize=10)
    fig.colorbar(cf, ax=ax, label="$\\omega_z D / U$", shrink=0.9)

    fig.savefig(out_dir / "cyl2d_re100_equilibrium.png", dpi=150)
    plt.close(fig)

    # --- leading eigenmode ---------------------------------------------------
    x, y, g = load_grid(run_dir / "eig0_vort_grid.pts")
    we = g["W_z"]
    we = we / np.nanmax(np.abs(we))
    fig, ax = plt.subplots(figsize=(10, 3.4), constrained_layout=True)
    cf = ax.contourf(x, y, np.clip(we, -0.5, 0.5), levels=np.linspace(-0.5, 0.5, 41), cmap="PuOr_r")
    add_cylinder(ax)
    ax.set_xlim(x[0], x[-1]); ax.set_ylim(y[0], y[-1])
    ax.set_title("Leading eigenmode about the equilibrium, λ = 0.127 ± 0.741i "
                 "(vorticity, real part, normalised, clipped at ±0.5)", fontsize=10)
    fig.colorbar(cf, ax=ax, shrink=0.9)
    fig.savefig(out_dir / "cyl2d_re100_eigenmode.png", dpi=150)
    plt.close(fig)

    print("wrote", out_dir / "cyl2d_re100_equilibrium.png")
    print("wrote", out_dir / "cyl2d_re100_eigenmode.png")

    # --- periodic orbit (optional) -------------------------------------------
    po, probe = run_dir / "po_vort_grid.pts", run_dir / "shed_probe.dat"
    if po.exists() and probe.exists():
        x, y, h = load_grid(po)
        t, vp = np.loadtxt(probe, unpack=True)
        fig = plt.figure(figsize=(10, 6.4), constrained_layout=True)
        gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1])
        ax = fig.add_subplot(gs[0])
        cf = ax.contourf(x, y, np.clip(h["W_z"], -3, 3), levels=np.linspace(-3, 3, 41),
                         cmap="RdBu_r")
        add_cylinder(ax)
        ax.set_xlim(x[0], x[-1]); ax.set_ylim(y[0], y[-1])
        ax.set_title("Shedding periodic orbit, Re = 100: vorticity at the Newton phase "
                     "(T = 5.98417, St = 0.1671)", fontsize=10)
        fig.colorbar(cf, ax=ax, label="$\\omega_z D / U$", shrink=0.9)

        ax = fig.add_subplot(gs[1])
        ax.plot(t, vp, color="C0", lw=0.8)
        ax.set_xlabel("t U / D")
        ax.set_ylabel("v at (2, 0)")
        ax.set_xlim(t[0], t[-1])
        ax.set_title("DNS from the perturbed equilibrium (growth rate 0.127) to the "
                     "limit cycle; the final state is the Newton initial guess",
                     fontsize=10)
        fig.savefig(out_dir / "cyl2d_re100_periodic_orbit.png", dpi=150)
        plt.close(fig)
        print("wrote", out_dir / "cyl2d_re100_periodic_orbit.png")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
