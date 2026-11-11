# nekecs tutorial: from a geometry to equilibria, orbits and branches

This tutorial assumes no experience with Nektar++. It takes one example, a
cylinder in a channel made with gmsh (`cases/cyl_gmsh/`), from the geometry
to a converged unstable steady state, its stability and its continuation in
Re, and explains every input file and every output file on the way. All
commands were run as written; the outputs shown are from those runs.

You need nekecs installed (see the [README](../README.md#installation)). Run
everything inside Linux (WSL on Windows), in a terminal.

1. [How the pieces fit together](#1-how-the-pieces-fit-together)
2. [The geometry and mesh](#2-the-geometry-and-mesh)
3. [The session file](#3-the-session-file)
4. [A starting guess](#4-a-starting-guess)
5. [Finding an equilibrium](#5-finding-an-equilibrium)
6. [Reading and viewing the results](#6-reading-and-viewing-the-results)
7. [Continuation in a parameter](#7-continuation-in-a-parameter)
8. [Periodic orbits](#8-periodic-orbits)
9. [ECS parameters: what to change when](#9-ecs-parameters-what-to-change-when)
10. [Troubleshooting](#10-troubleshooting)

The whole example runs with one script; the sections below go through it
step by step:

```bash
bash cases/cyl_gmsh/run.sh ~/runs/cyl_gmsh
```

## 1. How the pieces fit together

A Nektar++ run always has two inputs:

| File | Contains | Made with |
|---|---|---|
| **mesh file** (`mesh.xml`) | element vertices, edges, curved boundaries, and numbered groups of boundary edges ("composites") | gmsh + `NekMesh`, or any mesher NekMesh reads |
| **session file** (`session.xml`) | polynomial order, equations, time step, parameters such as Re, boundary conditions per composite, initial condition | a text editor |

and you always give both on the command line, mesh first:

```
<program> mesh.xml session.xml
```

The programs used here:

| Program | From | Does |
|---|---|---|
| `gmsh` | apt | geometry and mesh |
| `NekMesh` | Nektar++ | converts the gmsh mesh to `mesh.xml` |
| `IncNavierStokesSolver` | Nektar++ | ordinary time stepping (DNS) |
| `FieldConvert` | Nektar++ | post-processing: vorticity, `.vtu` for ParaView, interpolation |
| `nekecs` | this project | equilibria, periodic orbits, stability, continuation |

**What nekecs does.** Nektar++ can advance a velocity field u by a time T:
call that Φ_T(u). A steady state (equilibrium) is a field with
Φ_T(u) = u; a periodic orbit with period T is one with Φ_T(u) = u for that
T. nekecs solves these equations with Newton's method, so it also finds
**unstable** states, which time stepping never settles on. It needs a
starting guess near the answer, and it reports the residual
‖Φ_T(u) − u‖/‖u‖ at each Newton step: that number going to about 1e-10 means
success.

Fields (velocity and pressure) are stored in **`.fld` files**, Nektar++'s
own format. You view them by converting them to `.vtu` for ParaView
([section 6](#6-reading-and-viewing-the-results)).

## 2. The geometry and mesh

### 2.1 The gmsh geometry file

[`cases/cyl_gmsh/cylinder.geo`](../cases/cyl_gmsh/cylinder.geo) describes a
cylinder of diameter 1 at the origin in the box −10 ≤ x ≤ 25,
−10 ≤ y ≤ 10. Its four parts:

**Shapes.** A rectangle minus a disk:

```
Rectangle(1) = {xin, -ylim, 0, xout - xin, 2 * ylim};
Disk(2) = {0, 0, 0, 0.5};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
```

Any shape gmsh can build works the same way: several disks for several
cylinders, an airfoil from a spline through its points, a CAD file loaded
with `Merge "part.step";`.

**Element sizes.** Small on the body, larger downstream, large far away
(`h_cyl`, `h_wake`, `h_far` at the top of the file). Change these to refine
or coarsen the mesh.

**Element type.** Quadrilaterals (`Mesh.RecombineAll = 1`) of second order
(`Mesh.ElementOrder = 2`): the mid-edge nodes lie on the true circle, which
NekMesh turns into curved element edges. Spectral/hp methods need curved
boundaries: straight-sided elements on a circle cost accuracy at any
polynomial order. Triangles work too (drop `RecombineAll`).

**Physical groups.** These name the boundaries, and they are what the
session file refers to:

```
Physical Curve(1) = {wall[]};                     // cylinder wall
Physical Curve(2) = {inflow[], bottom[], top[]};  // inflow and sides
Physical Curve(3) = {outflow[]};                  // outflow
Physical Surface(10) = {3};                       // the fluid
```

The numbers are yours to choose. NekMesh keeps them: physical group 1
becomes composite `C[1]` in `mesh.xml`, and so on. Every boundary curve must
be in some physical group, and every boundary group gets its own boundary
condition later.

### 2.2 Making and converting the mesh

```bash
mkdir -p ~/runs/cyl_gmsh && cd ~/runs/cyl_gmsh
```

```bash
cp ~/nekecs/cases/cyl_gmsh/* .
```

```bash
gmsh -2 cylinder.geo -o cylinder.msh
```

```bash
~/opt/nektar++/bin/NekMesh cylinder.msh mesh.xml
```

`gmsh -2` makes the 2D mesh (it reports `5744 nodes 1516 elements`: the
1411 quadrilaterals plus the boundary edges). `NekMesh` writes `mesh.xml`.
To check it worked, look at the composites:

```bash
grep -o '<C ID="[0-9]*"> [A-Z]' mesh.xml
```

```
<C ID="1"> E
<C ID="2"> E
<C ID="3"> E
<C ID="10"> Q
```

Three groups of edges (E) and one of quadrilaterals (Q), as defined in the
`.geo` file. To see the mesh, convert it for ParaView:

```bash
~/opt/nektar++/bin/FieldConvert mesh.xml mesh.vtu
```

and open `mesh.vtu` in ParaView (from Windows: `\\wsl$\Ubuntu\home\<you>\runs\cyl_gmsh`),
with "Surface With Edges".

## 3. The session file

[`cases/cyl_gmsh/session.xml`](../cases/cyl_gmsh/session.xml) is commented
line by line; here is what each block means and what you would change.

### 3.1 Expansion: polynomial order

```xml
<E COMPOSITE="C[10]" NUMMODES="5" FIELDS="u,v,p" TYPE="MODIFIED"/>
```

On every element of composite 10 (the fluid), each field is a polynomial of
order `NUMMODES − 1` = 4. Raising `NUMMODES` is how you refine a spectral/hp
calculation without remeshing; check that your results do not change when
you do (a resolution study). Keep `TYPE="MODIFIED"`.

### 3.2 Solver: which equations

```xml
<I PROPERTY="SolverType"        VALUE="VelocityCorrectionScheme"/>
<I PROPERTY="EQTYPE"            VALUE="UnsteadyNavierStokes"/>
<I PROPERTY="EvolutionOperator" VALUE="Nonlinear"/>
<I PROPERTY="Projection"        VALUE="Galerkin"/>
```

plus `<METHOD> IMEX </METHOD> <ORDER> 2 </ORDER>`. **Keep these exactly:**
nekecs is built on (and checked for) this scheme only.

### 3.3 Parameters

```xml
<P> Re       = 60      </P>
<P> Kinvis   = 1.0/Re  </P>
<P> TimeStep = 0.005   </P>
```

Everything is non-dimensional: lengths in body diameters D, velocities in
the inflow speed U, time in D/U. `Kinvis` (the kinematic viscosity ν) is
what Nektar++ uses; with U = D = 1 it is 1/Re. Parameters can be
expressions and can use earlier parameters.

`TimeStep` must be small enough for stability (a CFL condition). If a run
produces `NaN` or blows up, halve it. Smaller elements or a higher
`NUMMODES` need a smaller step.

The `ECS_*` parameters are for nekecs; Nektar++ ignores them. They are
explained where they are used ([section 5](#5-finding-an-equilibrium),
[7](#7-continuation-in-a-parameter), and the table in
[section 9](#9-ecs-parameters-what-to-change-when)).

### 3.4 Boundary regions and conditions

First, **regions**: give each mesh composite a region number.

```xml
<B ID="0"> C[1] </B>  <!-- cylinder wall -->
<B ID="1"> C[2] </B>  <!-- inflow and the two sides -->
<B ID="2"> C[3] </B>  <!-- outflow -->
```

Then one **condition** per region and per variable: `D` = Dirichlet (the
value is prescribed), `N` = Neumann (the normal derivative is prescribed).

```xml
<REGION REF="0">  <!-- no slip -->
    <D VAR="u" VALUE="0"/>
    <D VAR="v" VALUE="0"/>
    <N VAR="p" USERDEFINEDTYPE="H" VALUE="0"/>
</REGION>
```

`VALUE` can be an expression in `x`, `y` and the parameters. Common cases:

| Boundary | u | v | p |
|---|---|---|---|
| no-slip wall | `D 0` | `D 0` | `N 0`, `USERDEFINEDTYPE="H"` |
| uniform inflow | `D 1` | `D 0` | `N 0`, `USERDEFINEDTYPE="H"` |
| inflow at angle `AoA` (radians, a parameter) | `D cos(AoA)` | `D sin(AoA)` | `N 0`, `USERDEFINEDTYPE="H"` |
| cylinder rotating at rate `Omega` (counter-clockwise) | `D -Omega*y` | `D Omega*x` | `N 0`, `USERDEFINEDTYPE="H"` |
| slip wall / symmetry line along x | `N 0` | `D 0` | `N 0`, `USERDEFINEDTYPE="H"` |
| outflow | `N 0` | `N 0` | `D 0` |

The `USERDEFINEDTYPE="H"` on the pressure marks the high-order pressure
condition that the velocity-correction scheme uses wherever the velocity is
prescribed; always use it there. The rotating-wall line assumes the body is
centred at the origin (wall speed Omega × radius).

### 3.5 Initial condition

```xml
<FUNCTION NAME="InitialConditions">
    <F VAR="u,v,p" FILE="mesh.fld"/>
</FUNCTION>
```

reads the fields from a file. An expression works too
(`<E VAR="u" VALUE="1"/>`, as in `dns.xml`). **For nekecs the initial
condition is the Newton starting guess.**

## 4. A starting guess

Newton needs a guess close to the answer. For a steady wake, a short time
integration from a uniform stream is enough: the flow near the body settles
quickly, while the instability (if Re is above the critical value) has not
grown yet. [`dns.xml`](../cases/cyl_gmsh/dns.xml) is the same physics with
`FinTime = 40` and a uniform initial condition:

```bash
export OPENBLAS_NUM_THREADS=1
```

```bash
~/opt/nektar++/bin/IncNavierStokesSolver mesh.xml dns.xml
```

It prints a short summary of the setup, then a progress line every
`IO_InfoSteps` steps:

```
Steps: 1000     Time: 5            CPU Time: 21.2596s
Steps: 2000     Time: 10           CPU Time: 20.7634s
...
Steps: 8000     Time: 40           CPU Time: 23.1653s
```

(`CPU Time` is per 1000 steps), then a timing table. This takes about 3 min
and writes `mesh.fld`: Nektar++ names output after the **first** file on the
command line. `session.xml` reads
`mesh.fld` as its initial condition.

`OPENBLAS_NUM_THREADS=1` matters: multithreaded OpenBLAS makes these runs
about 8 times slower. Put the `export` in your `~/.bashrc`.

## 5. Finding an equilibrium

```bash
~/build/nekecs/nekecs mesh.xml session.xml
```

The settings that matter (in `session.xml`):

```xml
<P> ECS_Steps     = 200    </P>  <!-- T = 200 x 0.005 = 1 -->
<P> ECS_Tol       = 1e-10  </P>
<P> ECS_Eigs      = 2      </P>
```

nekecs looks for u with Φ_T(u) = u, where T = `ECS_Steps` × `TimeStep`.
The answer does not depend on T (a steady state is steady for any T); T only
affects how fast Newton converges. T of order 1 to 5 D/U is a good range.

### Reading the log

The log (on screen; `run.sh` saves it in `nekecs.log`) for this run, which
took 16 min on one core:

```
[nekecs] equilibrium search: nvel = 2, steps/call = 200, dt = 5.000e-03, T = 1.000e+00, start = constant-history, |x0| = 2.671e+01, |F(x0)|/|x0| = 1.011e-02
[newton] it 1  gmres k = 37  lin.res = 9.351e-04  |dx_newton| = 5.273e-01
[nekecs] it 1  |F|/|x| = 5.254e-04  flow-map calls = 41
[newton] it 2  gmres k = 44  lin.res = 9.046e-04  |dx_newton| = 8.251e-02
[nekecs] it 2  |F|/|x| = 6.925e-06  flow-map calls = 87
[newton] it 3  gmres k = 41  lin.res = 8.498e-04  |dx_newton| = 1.226e-03
[nekecs] it 3  |F|/|x| = 6.085e-09  flow-map calls = 130
[newton] it 4  gmres k = 38  lin.res = 8.605e-04  |dx_newton| = 4.413e-07
[nekecs] it 4  |F|/|x| = 5.235e-12  flow-map calls = 170
[nekecs] status = converged, newton its = 4, rejected = 0, flow-map calls = 172
[nekecs] |F|/|x| = 5.235e-12, |Phi_4T(x) - x|/|x| = 1.408e-11, |u - u0|/|u0| = 1.891e-02
[nekecs] wrote session_eq.fld
[nekecs] stability: converged, Krylov dim = 55, flow-map calls = 55, T = 1.0e+00, levels = 2
[nekecs]   i            mu            |mu|    lambda = log(mu)/T        residual
[nekecs]   0  +7.41956e-01 -7.48377e-01i  1.05383e+00  +5.24355e-02 -7.89706e-01i  2.7e-05
[nekecs]   1  +7.41956e-01 +7.48377e-01i  1.05383e+00  +5.24355e-02 +7.89706e-01i  2.7e-05
```

| Line | Meaning |
|---|---|
| `\|F(x0)\|/\|x0\|` | how far the guess is from steady: 1e-2 here |
| `[nekecs] it n  \|F\|/\|x\|` | residual after Newton step n. It should fall fast: 1e-2, 5e-4, 7e-6, 6e-9, ... (the exponent roughly doubling is the sign of Newton converging) |
| `gmres k` | Krylov vectors used to solve the Newton step (each costs one flow-map call) |
| `lin.res` | how well that linear solve succeeded; it aims for 1e-3 (`ECS_GmresTol`) |
| `flow-map calls` | total cost so far: each call is `ECS_Steps` time steps |
| `status = converged` | reached `ECS_Tol` |
| `\|Phi_4T(x) - x\|/\|x\|` | independent check: the state is still steady after 4 more periods |
| `stability:` and the eigenvalue table | leading eigenvalues λ (`lambda = log(mu)/T`): **Re(λ) > 0 means unstable**; Im(λ) is the angular frequency of the growing mode. Here λ = 0.0524 ± 0.790i: the steady wake at Re = 60 is unstable, so a DNS would leave it and start shedding vortices. `residual` is the accuracy of each eigenvalue |

## 6. Reading and viewing the results

### 6.1 Output files

| File | What |
|---|---|
| `session_eq.fld` | the equilibrium (u, v, p) |
| `session_eig0_re.fld`, `session_eig0_im.fld`, `..._eig1_...` | real and imaginary parts of the eigenvectors (a complex pair gives eig0 and eig1 as conjugates) |
| `session_ecs_latest.fld` | the last accepted Newton iterate, updated every step: restart from it if a run stops early (point the initial condition at it) |

Names come from the session file name: `cont.xml` gives `cont_...`.

### 6.2 Viewing in ParaView

Convert to `.vtu`, adding vorticity on the way (`-m vorticity`; FieldConvert
needs the mesh and session files to know the expansion):

```bash
~/opt/nektar++/bin/FieldConvert -m vorticity mesh.xml session.xml session_eq.fld session_eq.vtu
```

Open `session_eq.vtu` in ParaView and colour by `u` or `W_z` (the
vorticity). Eigenvectors the same way:

```bash
~/opt/nektar++/bin/FieldConvert -m vorticity mesh.xml session.xml session_eig0_re.fld eig0.vtu
```

Other useful FieldConvert modules: `-m interppoints` (values on a grid or a
line, for plotting in Python; see `tools/make_cyl2d_figures.sh`), `-m
wss` (wall shear stress). `FieldConvert -l` lists them all.

### 6.3 Values in Python

A `.fld` is easiest to read in Python after interpolating onto points. For
example the centreline velocity behind the body:

```bash
python3 -c "import numpy as np; np.savetxt('line.csv', np.c_[np.linspace(0.6, 15, 200), np.zeros(200)], delimiter=',', header='x,y', comments='')"
```

```bash
~/opt/nektar++/bin/FieldConvert -m interppoints:fromxml=mesh.xml,session.xml:fromfld=session_eq.fld:topts=line.csv line_values.csv
```

`line_values.csv` then has columns x, y, u, v, p (header line `# x,y,u,v,p`).
For example, the length of the recirculation bubble (where u < 0 behind the
body):

```bash
python3 -c "import numpy as np; d = np.loadtxt('line_values.csv', delimiter=',', comments='#'); x, u = d[:, 0], d[:, 2]; print('bubble length', x[u < 0].max() - 0.5, 'D')"
```

```
bubble length 3.50101 D
```

## 7. Continuation in a parameter

Continuation follows a solution while a parameter changes: start from a
converged solution, step the parameter, re-converge, repeat. nekecs uses
pseudo-arclength continuation, so it also follows a branch round a fold
(where it turns back in the parameter).

[`cases/cyl_gmsh/cont.xml`](../cases/cyl_gmsh/cont.xml) is `session.xml`
with three changes:

1. which parameter to vary, in `SOLVERINFO`:

   ```xml
   <I PROPERTY="ECS_ContParam" VALUE="Re"/>
   ```

2. how far and how fast, in `PARAMETERS`:

   ```xml
   <P> ECS_ContEnd   = 40  </P>  <!-- stop after passing Re = 40 -->
   <P> ECS_ContDmu   = -2  </P>  <!-- first step: Re 60 -> 58 (sign = direction) -->
   <P> ECS_ContDs    = 4   </P>  <!-- then steps of about 4 in Re -->
   <P> ECS_ContDsMax = 8   </P>  <!-- never more than about 8 -->
   ```

   The start value is the parameter's value in the file (Re = 60). The step
   grows when Newton converges easily and shrinks (and retries) when it does
   not.

3. the initial condition is the converged equilibrium `session_eq.fld`.

```bash
~/build/nekecs/nekecs mesh.xml cont.xml
```

**Which parameters can be continued:** `Re`; `Kinvis`; and any parameter
that appears **directly** in a boundary-condition expression, such as
`AoA` in `VALUE="cos(AoA)"` or `Omega` for a rotating cylinder. Use the name
exactly as written in the session file. Geometric parameters (gap ratio,
body shape) are not supported yet: they change the mesh.

### Reading the branch file

`cont_branch.dat` gets one line per converged point:

| Column | Name | Meaning |
|---|---|---|
| 1 | point | 0, 1, 2, ... (matches `cont_c000.fld`, `cont_c001.fld`, ...) |
| 2 | Re | the parameter value |
| 3, 4 | arclength, dmu/ds | position along the branch; dmu/ds changes sign at a fold |
| 5 | newton_its | Newton iterations for this point |
| 6, 7 | \|F\|/\|x\|, verify | residual and the independent check |
| 8 | \|u\| | norm of the velocity |
| 9, 10 | u_probe, v_probe | velocity at the probe point (`ECS_ProbeX`, `ECS_ProbeY`, default (2, 0)) |
| 11, 12 | T, St | period and Strouhal number (periodic orbits; 0 for equilibria) |
| 13 | stab_indicator | equilibria: largest Re(λ); orbits: largest \|μ\| − 1. **> 0: unstable** |
| 14, 15 | mode_re, mode_im | that eigenvalue (equilibria: λ; orbits: multiplier μ) |

The run of `cont.xml` (about 1.5 h on one core) prints one line per point
and, when the indicator changes sign, where the stability changed:

```
[nekecs] point 0: Re = 60, |F|/|x| = 5.235e-12, verify = 1.408e-11, indicator = 5.2436e-02, flow-map calls = 60, wrote cont_c000.fld
[nekecs] point 1: Re = 58, |F|/|x| = 1.022e-12, verify = 2.595e-12, indicator = 4.6064e-02, flow-map calls = 260, wrote cont_c001.fld
[nekecs] point 2: Re = 54.00374, |F|/|x| = 2.619e-11, verify = 5.546e-11, indicator = 3.2337e-02, flow-map calls = 412, wrote cont_c002.fld
[nekecs] point 3: Re = 48.010197, |F|/|x| = 8.975e-14, verify = 1.791e-13, indicator = 8.8968e-03, flow-map calls = 607, wrote cont_c003.fld
[nekecs] *** stability change between Re = 48.010197 and 42.018159, interpolated Re_c = 46.077063, lambda (-0.018680142,0.77003582) (Hopf, omega = 0.770036)
[nekecs] point 4: Re = 42.018159, |F|/|x| = 2.576e-13, verify = 7.486e-13, indicator = -1.8680e-02, flow-map calls = 790, wrote cont_c004.fld
[nekecs] point 5: Re = 36.028528, |F|/|x| = 1.578e-13, verify = 3.810e-13, indicator = -5.1239e-02, flow-map calls = 976, wrote cont_c005.fld
```

The steady wake is unstable at Re = 60 (indicator 0.052 > 0) and stable at
Re = 42 and below: it loses stability at Re_c ≈ 46.1. The eigenvalue that
crosses is complex (ω = 0.77), so this is a Hopf bifurcation, where vortex
shedding (a periodic orbit) starts. The value is interpolated over 6 units
of Re; to pin it down, solve directly (section 5) at two Re values just
either side and interpolate between those. The last line of the run repeats
each stability change as `[nekecs] summary: ...`.

Plot a branch with numpy:

```bash
python3 -c "import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; d = np.loadtxt('cont_branch.dat'); plt.plot(d[:,1], d[:,12], 'o-'); plt.axhline(0, color='k'); plt.xlabel('Re'); plt.ylabel('growth rate'); plt.savefig('branch.png')"
```

or with `tools/plot_branch.py`, which makes the two-panel figure of the
README.

**If a run stops** (power cut, killed job): every point is already saved.
Copy the session file, set `Re` to the last value in the branch file, point
the initial condition at that point's `cont_cNNN.fld`, and run again.

## 8. Periodic orbits

A periodic orbit needs a guess for the state **and** the period, both from
a DNS of the saturated oscillation. `cases/cyl2d` does this at Re = 100:

1. `shed mesh.xml.gz shed.xml`: DNS from the equilibrium plus a small
   perturbation until the oscillation saturates. It prints the period
   (`SHED_ProbeX/Y` set the probe point) and writes `shed_shed.fld`.
2. `nekecs mesh.xml.gz po.xml`, where `po.xml` has the initial condition
   `shed_shed.fld` and

   ```xml
   <P> ECS_Period = 5.9841596 </P>  <!-- the DNS period: > 0 means "periodic orbit" -->
   <P> ECS_Steps  = 600       </P>  <!-- time steps per period -->
   ```

   nekecs then solves for the state and the period together; the time step
   becomes `ECS_Period / ECS_Steps` (the `TimeStep` in the file is not
   used). The log shows `T = ...` at each step and `St = f D / U` at the
   end. With `ECS_Eigs > 0` it computes Floquet multipliers μ: one is always
   1 (shifting along the orbit), and the orbit is unstable if any other has
   |μ| > 1.

Continuation of orbits works the same way as for equilibria
(`cases/cyl2d/cont_po.xml`): add `ECS_ContParam` and the `ECS_Cont*` lines
to the orbit's session file.

## 9. ECS parameters: what to change when

All parameters, with defaults, are listed at the top of
[`apps/nekecs.cpp`](../apps/nekecs.cpp). The ones you will touch:

| Parameter | Default | Change it when |
|---|---|---|
| `ECS_Steps` | 100 | equilibria: sets T = Steps × TimeStep (aim for T ≈ 1–5). Orbits: steps per period |
| `ECS_Period` | 0 | 0: equilibrium. > 0: periodic orbit, with this period as the guess |
| `ECS_Tol` | 1e-10 | target residual. 1e-9 is plenty for orbits |
| `ECS_MaxNewton` | 30 | Newton step limit |
| `ECS_MaxKrylov` | 60 | if the log shows `gmres k` = this value and `lin.res` far above 1e-3, raise it or raise T |
| `ECS_Delta0` | off | a bad guess makes Newton wander: set e.g. 0.01 (a first step of at most 1% of ‖u‖) |
| `ECS_Eigs` | 0 | number of eigenvalues / Floquet multipliers wanted after convergence |
| `ECS_EigMaxDim` | 80 | if the stability line says `NOT converged`, raise it (cost: one flow-map call each) |
| `ECS_ContParam` (SOLVERINFO) | none | parameter to continue in |
| `ECS_ContEnd`, `ECS_ContDmu` | required | end value; first step (sign gives the direction) |
| `ECS_ContDs`, `ECS_ContDsMax` | \|Dmu\|, 4 × Ds | step along the branch and its maximum, in parameter units |
| `ECS_ContSteps` | 50 | maximum number of points |
| `ECS_ProbeX`, `ECS_ProbeY` | 2, 0 | probe point for the branch file |

## 10. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `NekMesh` or a solver: `Unable to read file` | the file is not there: gmsh failed (look at its output) or you are in the wrong directory |
| `Level 0 assertion violation` / `Composite not found.` | the session uses a `C[n]` the mesh does not have: compare `BOUNDARYREGIONS` with `grep -o '<C ID="[0-9]*"' mesh.xml` |
| run gives `NaN` or `inf` | time step too large: halve `TimeStep` (for orbits: raise `ECS_Steps`) |
| very slow | `export OPENBLAS_NUM_THREADS=1` |
| Newton residual stalls or grows; `trust region collapsed` | the guess is too far from a solution: run the DNS longer, or use the `<session>_ecs_latest.fld` of a previous attempt; try `ECS_Delta0 = 0.01` |
| residual falls only by a factor 10–30 per Newton step | GMRES is not solving the Newton step (`gmres k` at the limit): raise `ECS_Steps` (larger T) or `ECS_MaxKrylov` |
| `stability: NOT converged` | raise `ECS_EigMaxDim`; the eigenvalues printed are still estimates |
| `nekecs: only SolverType=VelocityCorrectionScheme is supported` | keep the `SOLVERINFO` of section 3.2 |
| continuation: `corrector failed, halving ds` repeatedly | the step is too large for the branch: lower `ECS_ContDs`, `ECS_ContDsMax` |
| unsure whether a new setup is sound | `~/build/nekecs/test_flowmap mesh.xml session.xml` checks that the flow map is repeatable and that its finite-difference Jacobian behaves; all lines should start with `ok` |
| MPI runs give erratic Newton convergence | add `<I PROPERTY="GlobalSysSoln" VALUE="XxtStaticCond"/>` to `SOLVERINFO` (or run serially) |
