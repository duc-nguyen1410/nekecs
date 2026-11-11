# nekecs

Exact coherent structures of the incompressible Navier–Stokes equations in
complex geometries. The solver runs Newton–Krylov iterations with a hookstep
trust region around the flow map of the [Nektar++](https://www.nektar.info)
spectral/hp time-stepper. It follows the same approach as ECSAct (Dedalus) and
Channelflow's `findsoln`, but works on unstructured meshes, such as flow past a
cylinder.

| Solution type | Status |
|---|---|
| Equilibria (EQ) | **working**: unstable steady wake of a 2D cylinder at Re = 100 |
| Periodic orbits (PO) | **working**: vortex-shedding cycle at Re = 100 |
| Stability / Floquet multipliers | **working**: Arnoldi on the same flow-map Jacobian |
| Traveling waves / relative POs | planned, in quasi-3D runs (Fourier in span) |

## Cylinder at Re = 100: the steady wake and the shedding cycle

### Equilibrium: the unstable steady wake

At Re = 100 the steady wake is unstable: it sheds vortices above Re ≈ 47, so
time stepping alone cannot reach it. `nekecs` finds it directly. The run starts
from a nearly steady field and converges in **4 Newton steps** (251 flow-map
calls, about 4 min on one core):

| Newton step | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| ‖Φ_T(u) − u‖ / ‖u‖ | 6.7e-5 | 2.1e-6 | 3.8e-8 | 6.1e-10 | **1.0e-11** |

The state also stays put under continued time stepping, not only over the
T = 1 used in the residual: ‖Φ_kT(u) − u‖ / ‖u‖ is 1.0e-11, 2.1e-11, 4.0e-11
and 7.3e-11 for k = 1, 2, 4, 8.

![Re = 100 cylinder equilibrium: streamwise velocity with streamlines, and vorticity](docs/figures/cyl2d_re100_equilibrium.png)

Checks on the converged field that don't use the Newton solver:

- **Symmetry:** the wake is symmetric even though nothing imposes it. The
  centreline cross-stream velocity |v| is at most 3.6e-6.
- **Recirculation bubble:** 6.19 D long, with a minimum centreline velocity of −0.203 U.
- **Linear stability:** the leading eigenvalue pair is

  | Method | Leading λ |
  |---|---|
  | nekecs, two-level state (exact time-stepper Jacobian), T = 1 | **0.127049 ± 0.740684i** |
  | nekecs, restarted integrator, extrapolated T → ∞ from T = 1, 2 | 0.127050 ± 0.740685i |
  | Nektar++ Modified Arnoldi, linearised NS, T = 1 | 0.12709 ± 0.740665i |

  The state is unstable, with St = 0.118, below the nonlinear shedding value,
  as expected for the steady base flow. Nektar++'s Arnoldi restarts its
  integrator at every interval, which biases λ by about 4e-5 at T = 1. The
  nekecs two-level state has no such restart, so it gives the extrapolated
  value directly.

![Leading eigenmode of the equilibrium](docs/figures/cyl2d_re100_eigenmode.png)

### Periodic orbit: vortex shedding

The shedding cycle starts from a short DNS (`shed`): the equilibrium plus
1e-3 noise is integrated until the instability saturates, which gives an
initial guess for the state and the period. Newton then solves for
x = (u_n, u_{n−1}, T) with 600 steps per period. It converges in **2 Newton
steps** (36 flow-map calls, about 4 min on one core):

| | |
|---|---|
| ‖F‖ / ‖x‖ | 6.9e-4 → 1.1e-6 → **5.4e-10** |
| Period | **T = 5.98417083, St = 0.167108** (DNS estimate 5.98416) |
| Closure after 2 periods in one call | 6.6e-10 |
| Floquet multipliers | **1.00000** (trivial mode, \|μ − 1\| ≈ 2e-7), 0.649, 0.579 ± 0.146i: a stable limit cycle |

![Shedding periodic orbit: vorticity snapshot and the DNS probe signal that produced the initial guess](docs/figures/cyl2d_re100_periodic_orbit.png)

### What makes the flow map exact

Nektar++'s second-order IMEX scheme is multistep. A naive flow map restarts it
at every call, and the first-order start-up step leaves an artefact of about
2e-8: Φ_2T ≠ Φ_T∘Φ_T. nekecs instead runs its own stepping loop and carries the
previous time level u_{n−1} in the state, so each call is exactly N regular
steps. Composing two half-period calls then matches one full call to 2e-14,
even for a strongly unsteady state. Equilibria may also use a constant-history
start instead, which is cheaper. Details and all measurements are in
[DESIGN.md](DESIGN.md).

## Layout

```
core/      header-only C++17: GMRES, hookstep, Newton trust region, Arnoldi (no dependencies)
nektar/    Nektar++ backend: vector space over velocity fields, exact flow map, variable dt
apps/      nekecs (EQ/PO solver + stability), shed (DNS for initial guesses), diag_startup, probe
cases/     cyl2d: equilibrium, stability, shedding DNS and periodic-orbit sessions, setup script
tools/     build_nektar.sh, make_cyl2d_figures.sh, plot_cyl2d.py
docs/      figures
DESIGN.md  formulation, algorithms, every finding and measurement so far
```

## Build and run (WSL / Linux)

1. Build Nektar++ v5.10.0 with the solver libraries enabled. It installs to `~/opt/nektar++`.

   ```bash
   bash tools/build_nektar.sh
   ```

2. Build nekecs.

   ```bash
   cmake -S . -B ~/build/nekecs -DCMAKE_BUILD_TYPE=Release -DNektar++_DIR=$HOME/opt/nektar++/lib64/nektar++/cmake
   ```

   ```bash
   make -C ~/build/nekecs -j8 && (cd ~/build/nekecs && ctest)
   ```

3. Set up the cylinder case. Run serially with single-threaded OpenBLAS:
   threaded BLAS is about 8× slower here, and MPI doesn't pay off on this small
   mesh (see DESIGN.md).

   ```bash
   bash cases/cyl2d/setup.sh ~/runs/cyl2d
   ```

   ```bash
   cd ~/runs/cyl2d && export OPENBLAS_NUM_THREADS=1
   ```

4. Equilibrium, followed by its stability (writes `session_eq.fld`):

   ```bash
   ~/build/nekecs/nekecs mesh.xml.gz session.xml
   ```

   Optional cross-check with Nektar++'s own Arnoldi on the linearised equations:

   ```bash
   ~/opt/nektar++/bin/IncNavierStokesSolver mesh.xml.gz stability.xml
   ```

5. Periodic orbit: run the DNS to the limit cycle, which prints the period
   estimate, then Newton with Floquet multipliers. `po.xml` already contains
   the period estimate as `ECS_Period`.

   ```bash
   ~/build/nekecs/shed mesh.xml.gz shed.xml
   ```

   ```bash
   ~/build/nekecs/nekecs mesh.xml.gz po.xml
   ```

6. Make the figures (FieldConvert for vorticity and grid interpolation, then matplotlib):

   ```bash
   bash tools/make_cyl2d_figures.sh ~/runs/cyl2d docs/figures
   ```

Solver settings are `ECS_*` parameters in the session file; they are listed in
[apps/nekecs.cpp](apps/nekecs.cpp).
