# nekecs — design

A toolkit for computing exact coherent structures (ECS) of the incompressible
Navier–Stokes equations in complex geometries: equilibria (EQ) and periodic
orbits (PO), their linear stability, and their continuation in parameters.
Modeled on ECSAct (Dedalus-based), with Nektar++ as the time-stepper.

Scope (decided 2026-11-10): EQ and PO only. Traveling waves and relative
periodic orbits need a continuous symmetry (a periodic direction), which
complex inflow/outflow domains do not have; discrete symmetries
(shift-reflect) are equally geometry-specific. Neither is planned.

## 1. Problem formulation

Let `Φ_T(u)` be the Nektar++ flow map: integrate the incompressible NS from
velocity field `u` for time `T`. We seek

    F(u, T) = Φ_T(u) − u = 0

| Solution | Unknowns    | Phase condition                          |
|----------|-------------|------------------------------------------|
| EQ       | u           | none (T fixed, chosen for conditioning)  |
| PO       | u, T        | temporal: ⟨f(u_ref), u − u_ref⟩ = 0      |

where `f(u)` is the NS right-hand side (≈ (Φ_dt(u) − u)/dt) at a reference
state. In practice the state also carries the previous time level u_{n−1}
(see "Two-level history state").

Inner product: `⟨a, b⟩ = ∫_Ω a·b dΩ`, evaluated with the spectral/hp quadrature
(Nektar++ `Integral`), never a plain coefficient dot product. Perturbations
vanish on Dirichlet boundaries because every state is projected with the
boundary data imposed.

## 2. Algorithm

Newton–Krylov–hookstep (Viswanath 2007):

1. GMRES on `J δx = −F` with matrix-free `J v ≈ [F(x + εv) − F(x)]/ε`,
   with `ε = √eps_mach · (1 + ‖x‖)/‖v‖`. The linearised NS option (Nektar++
   `EvolutionOperator = Direct` around a time-dependent base flow) comes later.
2. After m GMRES iterations, `J V_m ≈ V_{m+1} H_m`. Hookstep: minimise
   `‖β e₁ − H_m y‖` subject to `‖y‖ ≤ Δ` via the SVD of `H_m` (bisect on μ in
   `y = Σ_i s_i b_i /(s_i² + μ) v_i`).
3. Accept or reject by comparing actual and predicted residual reduction; grow
   or shrink `Δ` (standard trust-region rules, as in Channelflow `findsoln`).
4. Converged when `‖F‖/‖u‖ < tol` (target 1e-10 for EQ, 1e-8 for PO).

Stability: Arnoldi on the same `J v` operator (or the existing Nektar++
`DriverModifiedArnoldi` with a converged base flow) → eigenvalues (EQ) and
Floquet multipliers (PO).

## 3. Making Φ_T differentiable

- Fixed `dt` with `T = N dt`; re-set `dt = T/N` whenever T changes. No CFL adaptivity.
- Every call to Φ starts from the same state: clear the IMEX multistep history
  and restart with the same start-up scheme. Otherwise Φ depends on hidden history.
- Pressure/velocity solver tolerances (`IterativeSolverTolerance`) at least ~1e-3
  below the Newton target, or use direct solvers in 2D.
- Mesh, polynomial order, and `dt` are frozen for the whole continuation.

## 4. Code layout

```
core/      Backend-agnostic, header-only C++17 (core/include/nekecs/): vector-space
           concept, GMRES, hookstep, Newton, Arnoldi, pseudo-arclength
           continuation. Tested on Lorenz and fold problems, no Nektar++ link.
nektar/    Nektar++ adapter (nektar/include/nekecs/nektar/flowmap.hpp): field
           vector and inner product, exact flow map on the velocity-correction
           scheme with its own time loop, variable dt, parameter changes.
apps/      nekecs (EQ/PO, stability, continuation), shed (DNS for PO guesses),
           diag_startup and diag_param (flow-map diagnostics).
cases/     cyl2d: session files, setup.sh, run_all.sh.
tools/     build_nektar.sh, figure scripts (make_cyl2d_figures.sh, plot_cyl2d.py,
           plot_branch.py).
```

### core/ (implemented 2026-11-10)
Header-only C++17, `core/include/nekecs/`, no dependencies:

| Header | Contents |
|---|---|
| `space.hpp` | Vector-space concept (`zeros_like`, `copy`, `dot`, `axpy`, `scale`) and `StdVectorSpace` |
| `dense.hpp` | `Matrix`, one-sided Jacobi `svd` |
| `gmres.hpp` | Unrestarted GMRES, MGS twice; returns basis, (k+1)×k Hessenberg, LS solution |
| `hookstep.hpp` | `HookstepSolver(H, beta).solve(delta)`: SVD once per Krylov space, bisection on μ |
| `newton.hpp` | `newton_hookstep(space, F, x0, options[, jacobian])`: FD or user J·v, trust-region accept/reject/grow, `on_accept` callback (return true to re-evaluate F after moving a phase reference) |

Tests (`core/tests`, run by `ctest`):
- `test_linalg`: SVD on cond ≈ 1e10 matrix, GMRES vs dense solve (err 4e-13),
  Arnoldi relation, unconstrained hookstep = GMRES solution, constrained
  hookstep on the trust-region boundary, predicted residual exact.
- `test_lorenz`: EQ via Φ_T − I (4 iterations); POs from recurrence guesses with
  T = 1.55865 (AB, matches Viswanath 2003), 2.30591 and 3.1173, each in
  4 iterations, closing to < 1e-8 at 4× finer dt; perturbed AB guess
  (|F0| = 9.3) converges; far guess (|F0| = 31.8) decreases monotonically with
  18 rejected/retried steps (it stalls at 0.17, as a local method should).

Eigensolver (added 2026-11-10):

| Header | Contents |
|---|---|
| `eig.hpp` | Complex shifted-QR eigenvalues of a small Hessenberg matrix (Wilkinson shifts, Givens), inverse-iteration eigenvectors |
| `arnoldi.hpp` | `arnoldi(space, op, v0, {nev, max_dim, tol})`: unrestarted Arnoldi (MGS twice), leading Ritz pairs by \|θ\|, residual \|h_{k+1,k} y_k\|, Ritz vectors as re/im Space vectors |

`test_arnoldi`: dense trace/det to 1e-14 and eigvector residuals ≤ 1e-9 on 20
random Hessenberg matrices; 80×80 A = QDQᵀ: 0.9 ± 0.3i, 0.8, 0.5 recovered to
2e-15 at k = 30; Lorenz C⁺ multipliers e^{λT} vs the characteristic cubic to
2e-12; AB orbit Floquet multipliers 4.713, 1 (to 1e-13), 1.19e-10, product =
e^{−(σ+1+β)T} to 4e-6 (limited by ε‖A‖ on the smallest multiplier). The
Lorenz eigen-tests use the exact tangent-linear flow (`lorenz::flow_tangent`),
so they test Arnoldi, not finite differences.

Not yet in core: restarted Arnoldi. (Continuation: see "Parameter continuation".)

### Integration with Nektar++ (verified 2026-11-10, v5.10.0)
Out-of-tree works, no fork needed. (The first probe of this, `apps/probe/`, was
removed once the adapter replaced it; it is in the git history.)

- Build Nektar++ with `NEKTAR_BUILD_SOLVER_LIBS=ON`; it installs
  `libIncNavierStokesSolver.so`, its headers under
  `include/nektar++/IncNavierStokesSolver/`, and `Nektar++Config.cmake`
  (which puts `IncNavierStokesSolverLib` in `NEKTAR++_LIBRARIES`).
- Link with `LINKER:--no-as-needed`: the solver classes register with the
  factory from static initialisers, so `--as-needed` silently drops the library
  ("No such module" at run time).
- Factory key is `SolverType` (`VelocityCorrectionScheme`, `VCSImplicit`,
  `CoupledLinearisedNS`, ...), not `EqType`.
- **Hidden history.** `UnsteadySystem::DoSolve` re-initialises the IMEX scheme
  each call, but the VCS high-order pressure-BC history (`Extrapolate::
  m_pressureHBCs`, `m_iprodnormvel`, `m_pressureCalls`) persists. Measured on the
  Re=100 cylinder (P=6, dt=0.01, 20 steps): repeat-call difference 1.1e-9
  without reset and 8.7e-15 with `m_extrapolation->GenerateHOPBCMap(m_session)`
  before each call. The adapter therefore subclasses `VelocityCorrectionScheme`
  (protected access to `m_extrapolation`) and constructs it directly, not
  through the factory. Other schemes (VCSImplicit, mapping, moving bodies) need
  the same check before use.

### nektar/ adapter and `nekecs` app (implemented 2026-11-10)
- `nektar/include/nekecs/nektar/flowmap.hpp`: `FieldVector` (velocity phys
  arrays + scalar extras), `NektarSpace` (quadrature L2 via `ExpList::Integral`,
  MPI-reduced), `VCSFlowMap::advance(u, nsteps, out)` (project u with Dirichlet
  values imposed, zero p, reset pressure-BC history, `DoSolve`).
- `apps/nekecs.cpp`: equilibria, F(u) = Φ_T(u) − u; `ECS_*` session parameters.
- `nektar/tests/test_flowmap.cpp`: determinism 1.3e-14; FD Jacobian converges at
  first order (10.0× per decade) to 6.7e-7 at ε = 1e-7‖u‖ (default `ECS_FDEps`),
  2e-6 at 1e-8; projection idempotent 1.4e-15.
- **cyl2d Re=100 equilibrium** (numbers below are from the first run, with the
  standard start-up; see "Own time loop" for the constant-history re-run and
  the bias-corrected eigenvalue) (`cases/cyl2d`, T = 1, 60 Krylov vectors, start
  from the Nektar++ SFD restart file, which is 2.5e-4 away): ‖F‖/‖u‖
  3.6e-3/53.96 → 2.1e-6 → 3.7e-8 → 5.6e-10 → 1.0e-11 in 4 Newton steps, 251
  flow-map calls, about 4 min serial. GMRES used all 60 vectors every step
  (linear residual 2–3e-2, not 1e-3); Newton still contracts ~50× per step.
  Independent checks on `session_eq.fld`:
  - wake symmetric without imposing it: max |v| on centreline 3.6e-6,
    |u(y=1) − u(y=−1)| ≤ 1.3e-4 (unstructured mesh);
  - recirculation length 6.19 D, min centreline u = −0.203;
  - Nektar++ Modified Arnoldi (`stability.xml`, linearised NS about our field):
    leading pair λ = 0.1271 ± 0.7407i, i.e. unstable with St = 0.118, below the
    nonlinear shedding St ≈ 0.164 as expected for the base flow.
  - nekecs' own Arnoldi (`ECS_Eigs`, FD Jacobian of Φ_T, T = 1): leading
    multipliers 0.838036 ± 0.766225i, λ = log μ / T = **0.127090 ± 0.740665i**,
    Ritz residual 1e-5 after 80 flow-map calls; identical to Nektar++'s result
    in every printed digit. The next pair (|μ| = 0.9927, λ ≈ −0.0074 ±
    0.760i) only reached residual 3.5e-3 at k = 80: weakly damped, clustered
    modes need a larger Krylov space or a restarted Arnoldi.
- **Run with `OPENBLAS_NUM_THREADS=1`.** Threaded OpenBLAS on the small elemental
  matrices: 94 s for 10 flow-map calls (sys time 2× user); single-threaded:
  11.6 s. Parallelise with MPI.

### Start-up artefact (IMEX order ≥ 2)
Every flow-map call restarts the multistep integrator, whose first step uses
the order-1 scheme; the VCS steady state of an order-1 step differs from that of
an order-2 step. So Φ_T is not a semigroup: on cyl2d, ‖Φ_2T(u) − Φ_T(Φ_T(u))‖/‖u‖
= 2.1e-8 with IMEX order 2 and 1.1e-14 with order 1 (`apps/diag_startup.cpp`).
Pressure zeroing and pressure-BC history make no difference (identical numbers
with both kept). Consequences:
- A Newton fixed point of Φ_T is invariant under *further* time stepping only to
  ~3e-8 (Φ_4T check), however small ‖F‖ is.
- For scale: the order-1 and order-2 VCS steady states differ by ~4e-6, i.e. the
  scheme's own splitting error is ~100× larger than the artefact.
- A burn-in residual F = Φ_{s+T}(u) − Φ_s(u) does not fix it: ‖F‖ floors at
  1.3e-8 and GMRES cannot reduce it (the start-up leaves an unreachable
  component). Removed.
- Fixed (2026-11-10) by our own time loop, below.

### Own time loop and start modes (implemented 2026-11-10)
`FlowMapVCS::Integrate(nsteps, StartMode)` replaces `DoSolve` in the flow map:
`InitializeScheme`, then per step `v_PreIntegrate`, `TimeIntegrate`, copy-back
with `FwdTransLocalElmt`, `v_PostIntegrate` (no filters, IO or CFL control).
- `StartMode::Standard`: Nektar++'s phases (order-1 first step). Matches
  `DoSolve` to 1.3e-14.
- `StartMode::ConstantHistory`: after `InitializeScheme`, the IMEX solution
  vector `[y_n..y_{n-k+1}, dtN_n..dtN_{n-k+1}]` gets y_{n-j} = u and
  dtN_{n-j} = dtN(u); the ODE RHS is evaluated at u k−1 more times so the VCS
  pressure-BC history (`m_pressureHBCs`, `m_iprodnormvel`, `m_pressureCalls`)
  is full and constant; every step uses the final order-k phase.

Results on cyl2d (fresh run from the SFD restart file, `ECS_Start = 1`, now the
default): Newton 4 steps to 1.05e-11 as before; invariance under continued
stepping ‖Φ_kT(u*) − u*‖/‖u*‖ = 1.0e-11, 2.1e-11, 4.0e-11, 7.3e-11 for
k = 1, 2, 4, 8 (was 1e-11, 2e-8, 3e-8, 4e-8); semigroup defect at u* 5.7e-14
(was 2.1e-8). So u* is now a steady state of the regular order-2 scheme.

**Stability needs the opposite choice.** Perturbations grow in time, so a
constant history is an O(dt λ) error for them, while the standard start-up is
second-order consistent. A per-call start-up error scales like 1/T, while the
time-stepper's true eigenvalue does not depend on the Arnoldi interval T:

| Jacobian start | λ (T = 1) | λ (T = 2) | 2λ(2) − λ(1) |
|---|---|---|---|
| standard | 0.127090 + 0.740665i | 0.127070 + 0.740675i | 0.127050 + 0.740685i |
| constant history | 0.126434 + 0.736971i | 0.126742 + 0.738827i | 0.127050 + 0.740683i |

Both extrapolate to **λ = 0.12705 ± 0.74068i**. The standard start is off by
4e-5 at T = 1 (Nektar++'s Modified Arnoldi has the same start-up and the same
bias, which is why it agreed with us to 6 digits); constant history is off by
6e-4 / 4e-3. Defaults: `ECS_Start = 1` (Newton), `ECS_EigStart = 0` (Arnoldi),
`ECS_EigSteps` sets the Arnoldi interval separately (longer T, smaller bias).

**Periodic orbits:** neither mode is exact (constant history is first-order
inconsistent for a time-dependent state; the standard start-up puts one order-1
step in every period). Resolved by the two-level state below.

### Two-level history state (implemented 2026-11-10)
`StartMode::History`: the state carries all k IMEX levels, x = (u_n, u_{n−1})
for order 2 (`FieldVector::u` = [u_n components, u_{n−1} components]).
`Integrate` writes u_{n−1} and dt N(u_{n−1}) into the IMEX solution vector and
evaluates the ODE RHS at u_{n−1} then u_n, which leaves the VCS pressure-BC
history exactly as after regular stepping; on exit it returns the new
(u_{n+N}, u_{n+N−1}). A call is then exactly N regular steps:
- semigroup for an unsteady state (‖Φ_T(x) − x‖/‖x‖ = 5e-2): one call vs two
  half calls 2.0e-14; 1 + 1 steps vs 2 steps 3.3e-15; History(w, w) vs
  ConstantHistory(w) 2.0e-14 (`apps/diag_startup.cpp`).
- **Bug found on the way:** projecting u_{n−1} with `ContField::FwdTrans` into a
  zero-initialised coefficient array imposed *homogeneous* Dirichlet data
  (u = 0 at the inflow); `FwdTrans` takes the Dirichlet values from its output
  array. Seeding it with the field's coefficients fixed it (1-step error 6e-2
  → 1.6e-15).
- Stability with the two-level state (`ECS_EigStart = 2`, now the default):
  λ = 0.127049 ± 0.740684i at T = 1, equal to the T → ∞ extrapolation above
  to ~1e-6. No start-up bias.

Variable time step (`FlowMapVCS::SetDt`): sets the session `TimeStep`,
`m_timestep`, rebuilds the pressure-BC extrapolation (`SetUpExtrapolation()`,
protected), because `Extrapolate` keeps its own copy of TimeStep for the du/dt
boundary term, and clears the velocity Helmholtz matrix cache
(λ = 1/(a_ii dt ν)). Changing dt by 2% and back reproduces results to 2e-14.

### Periodic orbits (implemented 2026-11-10)
`nekecs` with `ECS_Period > 0`: x = (u_n, u_{n−1}, T), dt = T/N,
F = (Φ_T(x) − x, ⟨b, u_n − u_ref⟩), b = unit flow direction at the guess
(fixed Poincaré hyperplane). User Jacobian: FD in the state at fixed T (no dt
change, so no Helmholtz refactorisation per GMRES vector) plus the period
column ∂Φ/∂T, computed by FD in T once per Newton iterate.

cyl2d Re = 100 shedding orbit (`cases/cyl2d/shed.xml`, `po.xml`):
- initial guess: `apps/shed.cpp` DNS from the equilibrium + 1e-3 noise to
  t = 300 (two-level state, chunked calls = continuous stepping). Saturates by
  t ≈ 90; last 12 periods 5.98416 ± 1e-4 (crossing interpolation on a 0.1 grid).
- Newton, N = 600 steps/period, 60 Krylov max: ‖F‖/‖x‖ 6.9e-4 → 1.1e-6 →
  **5.4e-10** in 2 steps (GMRES k = 8, 17), 36 flow-map calls, ~4 min serial.
- **T = 5.98417083, St = 0.167108**, dt = 9.97362e-3; closure after 2 periods
  in one call 6.6e-10; DNS estimate agrees to 2e-6.
- Floquet multipliers (two-level Jacobian over one period): **1.00000**
  (|μ − 1| ≈ 2e-7, residual 2.6e-13: the time-translation mode), 0.649
  (residual 2.5e-4), 0.579 ± 0.146i (not converged at k = 60, residual 2e-2).
  Stable limit cycle, as expected for 2D shedding at Re = 100.

### Parameter continuation (implemented 2026-11-10)
`core/continuation.hpp`: pseudo-arclength continuation of G(x; μ) = 0 for
any problem the Newton solver handles (EQ: x = u; PO: x = (u_n, u_{n−1}, T)).
- z = (x, μ) in X × R with ⟨z, z'⟩ = ⟨x, x'⟩ + w μ μ'. nekecs sets
  w = ECS_ContWeight · (‖x0‖/μ0)², so a 1% change of the parameter weighs like
  a 1% change of the state; ECS_ContDs is then roughly the parameter step
  away from folds, in parameter units.
- Start: Newton at μ0, then a natural-parameter step to μ0 + Δμ0. Then secant
  predictor z1 + ds t and corrector Newton on (G(x; μ), ⟨t, z − z1⟩ − ds),
  which goes round folds. ds grows ×1.5 after ≤ 2 Newton iterations and
  shrinks ×0.5 after ≥ 5 or a failed corrector (retry from z1).
- Jacobian (J_x v_x + v_μ ∂G/∂μ, ⟨t, v⟩): ∂G/∂μ by FD once per Newton
  iterate, J_x at fixed μ (the user's operator, e.g. the PO Jacobian with its
  period column, or FD). GMRES never changes μ, so Nektar++ never refactorises
  the viscous matrices inside a Krylov solve.
- PO: the Poincaré hyperplane moves to each accepted point (on_point →
  rebase); the previous point still lies on it, so the secant stays valid.
- Core tests (`test_continuation`): saddle-node x' = μ − x² (round the fold
  at μ = 0 to the lower branch, all points equilibria to 1e-12, one fold);
  Lorenz AB orbit ρ = 28 → 32.3 (T 1.55865 → 1.41551, closure under 4× finer
  integration 6e-10, last point = direct fixed-ρ Newton to 1e-8).

Parameters in the Nektar++ adapter (`VCSFlowMap::set_parameter`):
- **Re / Kinvis** (`FlowMapVCS::SetKinvis`): ν enters λ = 1/(a_ii dt ν), the
  forcing scale 1/ν (`m_diffCoeff`) and the pressure BCs (`m_kinvis`, passed
  to the extrapolation at each step, not cached). Set both, drop the velocity
  matrix cache, update the field metadata. Refused with SVV or a
  DiffusionCoefficient function.
- **Boundary-data parameters** (inflow angle, rotation rate, blowing, ...):
  Nektar++ binds session parameters at expression *evaluation* (interpreter
  StorePrm), and `EvaluateBoundaryConditions(0.0)` re-evaluates every BC, so
  `SetParameter` + re-evaluation changes the BC in place. Limits: names are
  case-sensitive as written in the session; parameters derived from it in the
  session (P2 = 2*P1) were evaluated once at load and do not follow; forcing
  and other functions are not re-evaluated.
- **Geometric parameters** (gap ratio, ...): not yet. Plan: a mesh family
  x(ξ; μ) with fixed topology and order (vertices and curved nodes moved by a
  smooth displacement field), one solver per μ (μ and μ + h cached), state
  arrays carried over by DOF index. Alternative: Nektar++'s coordinate
  mappings (VCSMapping), which would need the same hidden-history audit as
  the VCS.
- Check (`apps/diag_param.cpp`, two-level Φ, 100 steps): in-place change vs a
  solver built from scratch with `-P` at the new value, and change-and-back
  vs the original:

  | Parameter | Effect on Φ | In place vs fresh | Back vs original |
  |---|---|---|---|
  | Re 100 → 80 | 8.0e-4 | 1.5e-14 | 1.3e-14 |
  | AoA 0 → 0.1 (far-field u = cos AoA, v = sin AoA) | 8.9e-2 | 1.5e-14 | 1.4e-14 |

  Note `-P Re=80` alone would not do: `Kinvis = 1.0/Re` is evaluated at load,
  before command-line overrides. The fresh solver gets `-P Kinvis=0.0125`.

### MPI
In parallel Nektar++ defaults to an iterative static-condensation solver: the
flow map is then not deterministic (6.6e-8 repeat-call difference), projection
is not idempotent (1.8e-8) and the FD Jacobian is useless. With
`GlobalSysSoln = XxtStaticCond` (parallel direct) everything passes again
(determinism 1.9e-14). But on this 746-element mesh XXT set-up dominates: 10
flow-map calls take 2 min (2 ranks), 34 s (4), 13 s (8) vs ~12 s serial. Run
small 2D cases serially; revisit for 3D.

## 5. Validation ladder and roadmap

1. `core/`: hookstep on Lorenz POs; GMRES vs dense solve on random matrices. **Done.**
2. 2D cylinder, Re = 100: unstable steady wake (EQ), leading eigenvalue
   unstable. **Done:** λ = 0.127049 ± 0.740684i.
3. 2D cylinder, Re = 100: shedding PO and Floquet multipliers. **Done:**
   St = 0.1671, trivial multiplier 1 to 2e-7.
4. Continuation in Re: EQ branch through the Hopf point (Re_c ≈ 46.6,
   ω ≈ 0.74 in the literature for unconfined domains), PO branch towards it.
   **In progress.**
5. Geometric-parameter continuation (gap ratio of two cylinders, ...).
6. An asymmetric 2D case (staggered cylinders, cylinder near a wall).
7. 3D: MPI with a deterministic solver, restarted Arnoldi, a GMRES
   preconditioner or linearised J·v.

## 6. Environment

Develop in WSL Ubuntu 24.04. System packages: build-essential, gfortran,
cmake, Boost, OpenBLAS/LAPACK, FFTW, OpenMPI, parallel HDF5, ARPACK (full
`apt install` line in the README).
Nektar++ v5.10.0 from source via `tools/build_nektar.sh`
(src `~/src/nektar`, install `~/opt/nektar++`). Notes:
- Scotch is off: Ubuntu's libscotch lacks PT-Scotch for MPI builds, and building
  it needs flex. METIS (bundled) does partitioning.
- Only IncNavierStokesSolver and ADRSolver are built (`NEKTAR_SOLVERS_ONLY`,
  a small idempotent patch to `solvers/CMakeLists.txt`). This also avoids a
  v5.10.0 bug where PulseWaveSolver's library target lists a missing header.
- 36 IncNS regression tests pass (cylinder SFD, quasi-3D cylinder, Modified
  Arnoldi, channel).

## 7. Open questions

- Nektar++ version to pin (latest v5.x release vs master).
- Finite-difference vs linearised Jacobian-vector products as the default.
- Geometric parameters: deformed mesh family vs Nektar++ coordinate mappings.
- Cheaper EQ correctors: at T = 1 GMRES stalls at 60 vectors (linear
  residual ~5e-2), so Newton converges only linearly. Options: larger T,
  Krylov recycling between Newton steps, a preconditioner.
