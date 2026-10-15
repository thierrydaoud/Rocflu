# Rocflu + ppiclF fixes on the `PseudoTurb` branch

The goal of this branch is to validate the subgrid pseudo-turbulence (PT) model. Before
the PT results can be trusted, the particle–gas coupling has to be correct: interpolating
the gas to the particles, projecting the particles back onto the mesh, the drag and heat
closures, particle collisions with each other and with walls, and conservation across MPI
ranks and periodic boundaries. This page lists what
was wrong, how each problem was fixed, and how the fix was checked. The PT model itself
(the "k_sg gets too high" issue) is not addressed yet; see [Still open](#still-open).

All changes are covered by the regression suite in [`tests/`](../tests/README.md):
35 cases, compared against stored goldens and checked against physics on 1, 2 and 4 ranks.

## Summary

| Area | Problem | Fix | Commit |
|---|---|---|---|
| Interpolation | Inverse distance (1/d³) had errors of up to ±18% of the change in density across one cell, even in a linear field, and was discontinuous at cell faces. A particle sheet hit by a planar shock rippled with the cell size (57 µm spread after 250 µs) | Weighted linear least squares, exact for linear fields; the sheet stays flat (0.07 µm) | 89f9661, 6d00ca4 |
| Projection | The kernel used one length for every direction, so the spread on stretched cells was wrong | Each direction is scaled by its own cell length, as in MFC | 414ce7f |
| Projection | Cells just outside the bin domain were dropped, so cloud-edge particles reached 8 cells instead of 27; φ was ~8% too high at the edge | Map those cells to particle bins | affad85 |
| Closures | Forces and heat transfer were skipped below 1 mm/s slip, and the old value was kept | Evaluate at any slip; skip only the quasi-steady drag at zero slip | f0ae123 |
| Quasi-steady drag | Gidaspow's C_D was capped below Re = 0.1 (100% error); Osnes had a spurious Knudsen reduction (25% error) | Aligned with MFC; QSFORCE renumbered | 94607fd |
| MPI | Gas state in virtual cells was stale at rank borders with particles, causing spurious gas motion in a quiescent cloud and a 1% difference between 1 and 2 ranks | Convert and re-exchange after φ is known; exchange φ and k_sg every stage | e3ed370 |
| Conservation | Fluxes on periodic and symmetry patches were not conservative | Exchange and transform gradients in virtual cells | 11dc640 |
| Robustness | Uninitialized reads gave NaN on 2+ ranks; a single particle gave a filter width of 0 | Initialize; filter width from cell–box overlap | b10c5f1, 48f64b7 |
| Build | Required `SPEC=1` | Builds without it; stops if species are requested | aa3c0b6 |
| Collisions | The contact stiffness used particle i's own mass, so an unequal pair felt unequal forces: an 8:1 pair at 20 m/s ended with its total momentum reversed | Stiffness from the pair's effective mass (MFC soft-sphere model) | 3afd23a, 78eb711 |
| Particle walls | Walls existed only through `filein.vtk`; a particle over the edge between two wall triangles was pushed by both (e = 0.60 instead of 0.70); a fast particle went through the wall and the run crashed | Walls from the fluid wall patches (`PARTWALL 1`), contacts computed like pair contacts, a fluid side that pushes back; `filein.vtk` removed | 2ff0cd8 |
| Ghost particles | Overflowing the ghost arrays wrote past their end | Guards stop the run; ghost arrays 30% smaller | f7a3980, 0bb2e4b |
| Neighbour search | Dense subbin tables: ~28 GB requested for 50,000 particles in one 2 cm bin | Compact lists: ~1.2 MB for the same bin (about 24,000 times less), byte-identical results | 82d4ca4 |

## Interpolation: weighted linear least squares

A tracer moves across cells in a gas whose density varies linearly in x. Each cell
centre holds the exact value, so a correct interpolation returns the straight line.

![Interpolation, old vs new](figures/fig_interpolation.png)

The old 1/d³ weighting pulled the value towards the nearest cell centre. The result was a
staircase with errors of up to ±18% of the density change across one cell, and it jumped
at every face. Weighted linear least squares (radius 1.8h) reproduces linear fields to
round-off, and its value is continuous as the particle crosses a face. Each cell also now
enters each bin's stencil only once; before, a two-bin periodic case counted some cells
twice.

### Interpolation in shock-driven flows

**A particle sheet ripples with the old interpolation.** A flat sheet of particles hit by a
planar shock should stay flat. In the LLNL sheet case (8,281 particles, dp = 150 µm,
Ms 1.27, cells 0.25 × 0.517 × 0.517 mm, run with the old code) it did not: particles near
cell faces lagged and particles near cell centres led, so the sheet became corrugated with
the cell size as wavelength. Its spread in x grew to 193 µm (1.3 diameters) by 0.5 ms.

The cause is the inverse-distance weighting. It weights each cell by its 3D distance from
the particle, so a particle near a y face, roughly equidistant from several cell centres,
gives more weight to the cells upstream and downstream of it in x than a particle at the
cell centre does. Behind the shock the gas varies in x, so particles at different places
in their cells saw different gas velocities and drag. Least squares is exact for a
variation in x whatever the particle's position in y and z.

A periodic slice of that case (4 × 4 cells across, 13 × 13 particles at 13 different
positions in their cells, the same forces) reproduces it, and shows that the new
interpolation removes it:

![Particle sheet, old vs new interpolation](figures/fig_sheet.png)

| Spread of the sheet in x (std) | 100 µs | 150 µs | 200 µs | 250 µs |
|---|---|---|---|---|
| Old, one-way (`FEEDBACK 0`) | 2.2 µm | 11.2 µm | 24.9 µm | 42.8 µm |
| Old, two-way (`FEEDBACK 1`) | 3.4 µm | 15.2 µm | 33.1 µm | 56.6 µm |
| New, one-way | 0.09 µm | 0.24 µm | 0.48 µm | 0.78 µm |
| New, two-way | 0.01 µm | 0.03 µm | 0.05 µm | 0.07 µm |

The ripple appears with one-way coupling too, where the gas stays uniform in y and z, so
it comes from the interpolation and not from the feedback. Panels (a) and (b) show each
particle's offset from the sheet's mean x against its place in the cell; the spread at a
given place comes from the place in z.

![Particle sheet, side view](figures/anim_sheet.gif)

**Shock profiles.** The gas seen by a line of stationary 1 µm tracers across a planar
shock (one-way), at two Mach numbers and two cell sizes:

![Planar shock, old vs new interpolation](figures/fig_shock_interp.png)

| Shock width 10–90% (density / pressure) | Cell values | Old | New |
|---|---|---|---|
| Ms 1.27, 1 mm cells | 3.69 / 3.29 mm | 3.86 / 3.67 mm | 3.49 / 3.31 mm |
| Ms 1.27, 0.5 mm cells | 1.94 / 1.91 mm | 1.99 / 1.99 mm | 1.92 / 1.92 mm |
| Ms 2.0, 1 mm cells | 3.77 / 2.86 mm | 3.68 / 3.00 mm | 3.32 / 2.81 mm |

"Cell values" is the straight line between neighbouring cell centres. Rocflu spreads the
shock over 3–4 cells, and that width halves with the cell size. Least squares follows the
cell profile and does not widen it (its support radius is 1.8 cells, but the fit is
linear). Neither method overshoots the range of the cell values. The old method has kinks
at every cell face, clearest in dp/dx (bottom row), so a particle crossing a face felt a
sudden change in its pressure-gradient force. Both methods converge to the same profile
when the cells are halved.

## Projection: kernel width per direction, and cloud-edge cells

A single particle is moved along x, and the panels show how its volume fraction φ is
shared among the cells (left) and the spread σ/Δl of that share (right). The top row uses
cube cells, the bottom row 2:1:1 cells.

![Projection, old vs new](figures/fig_projection.png)

- On 2:1:1 cells the old spread along the long side fell to 0.07 cell lengths: the
  particle's φ landed almost entirely in one cell. The kernel now scales each direction by
  the receiving cell's length, as in MFC, and the spread stays between 0.41 and 0.51 cell
  lengths in every direction.
- On cube cells the old curve has a notch at the cell centre. Cells whose centre lay just
  outside the bin domain were dropped, so a particle reached only 8–12 cells instead of 27.
  The same bug made φ about 8% too high at the edge of a particle cloud.

The kernel width is still about 0.46h, so it depends on the mesh. How sensitive k_sg is to
this width is the next physics question (see [Still open](#still-open)).

## Low-slip closures

![Low-slip closures, old vs new](figures/fig_closures.png)

- **(a) Drag.** When the slip fell below 1 mm/s, the old code skipped the force update and
  kept its stale value. The particle went on accelerating and overshot the gas, reaching
  1.2 times the gas velocity. The new code follows the Stokes solution.
- **(b) Heat.** A particle at rest in gas at rest exchanged no heat at all. Its temperature
  now relaxes as the lumped model with Nu = 2 predicts.
- **(c) Fox Nusselt number.** The old code gave Nu = 0 below 1 mm/s and floored Re_p at
  0.1. The new formula tends smoothly to 2.

## Quasi-steady drag, aligned with MFC

Rows 1–2: C_D from each code's Fortran routine at M → 0 and φ = 0, compared with the
published correlations. Row 3: a particle released at rest in gas at 1 m/s (Re₀ = 1),
compared with the ODE using the published correlation.

![Quasi-steady drag, old vs new](figures/fig_qs_drag.png)

| Model | Worst error before | After |
|---|---|---|
| Gidaspow (Wen–Yu) | 100% (Re floored at 0.1 inside C_D) | 8e-16 |
| Parmar et al. (2010) | 3% | 2e-10 |
| Osnes et al. (2023) / Loth et al. (2021) | 25% (M floored at 0.01, giving a spurious Knudsen reduction) | 1e-5 |

The remaining "after" errors come from evaluating at M = 1e-9 and from the code's own
M ≥ 1e-6 floor in Loth's terms, which MFC also has. The modified Parmar model was removed.

> **Input decks:** `QSFORCE` now follows MFC: 0 none, 1 Gidaspow, 2 Parmar, 3 Osnes.
> 4 stops with an error. Decks that used `QSFORCE 2` for Osnes must change to `3`.

## Particle collisions

Collisions use a soft-sphere model: overlapping particles push each other apart with a
spring and a dashpot. The spring stiffness is the softer of two values:

- **k₁ = m_eff (π² + ln² e) / (KSP·dt)²**. With this stiffness a contact lasts exactly
  KSP time steps and rebounds with restitution e = `EREST` (the linear model in MFC's
  `m_collisions`). dt is the current time step, so the contact time follows dt.
- **k₂ = (4/3) E\* √(R\*) √δ**, the Hertz stiffness at overlap δ, with E = 1 GPa and
  ν = 0.35 for every particle.

The dashpot is set from the chosen stiffness so that a linear contact rebounds with e.

### Contact stiffness from the pair's effective mass (3afd23a)

k₁ used particle i's own mass instead of the pair's effective mass m₁m₂/(m₁+m₂). With
unequal masses each particle of a pair used a different stiffness (and could pick a
different branch of min(k₁, k₂)), so the two contact forces were not equal and opposite.
An 8:1 pair at 2 m/s kept 33% of its momentum; at 20 m/s its total momentum reversed. For
equal masses k₁ was twice the intended value, and contacts lasted about 0.72·KSP steps
instead of KSP.

The test case `box/binary-collisions` (np 1, 2, 4) collides four pairs across a bin face,
a rank border and the periodic faces. Particles start with an initial velocity, which
`points.dat` now accepts as `material x y z dp u v w` (786ecb1).

![Binary collisions](figures/anim_binary_collisions.gif)

Each panel follows one pair around its centre of mass, drawn to scale. The line on each
particle turns with its spin. Below each panel, the normal relative velocity divided by
the approach speed rebounds to e = 0.697 for the head-on pairs (EREST = 0.7). Momentum is
conserved to 1e-5 and angular momentum about the centre of mass to 0.5%, including in
the oblique pair (panel 4), where friction sets both particles spinning.

### Ghost particles and the neighbour search (f7a3980, 0bb2e4b, 82d4ca4)

A particle near a bin face is sent to the neighbouring bins as a ghost, so that pairs
across bins, ranks and periodic faces are found. The ghost algorithm was compared with
Sam's ppiclF: ghosts carry the whole state vector, velocity included. A test build whose
ghosts dropped the velocity fails the new case `box/ghost-collisions` on 1, 2 and 4
ranks; the current code passes it. The subbin neighbour search was compared with a brute-
force search over all real and ghost particles (scratch build) and found every pair on
1, 2, 4 and 8 ranks.

Two robustness fixes came with this:
- Sending or receiving more ghosts than `PPICLF_LPART_GP` wrote past the end of the ghost
  arrays. The run now stops with a message. Ghosts no longer carry the projection values,
  which shrinks the ghost arrays by 30% (about 520 to 360 MB per rank at LPART = 80000).
- The subbin contents were dense tables of (subbins × particles) 4-byte integers. A 2 cm
  bin with 0.4 mm subbins has 140,608 subbins, so 50,000 particles requested about 28 GB
  per call. They are now compact lists: one start index per subbin plus three integers
  per particle, about 1.2 MB for the same bin (12 bytes more per ghost), about 24,000
  times less. The interpolation stencil search had the same kind of table over cells, which
  grew with the square of the cells per bin; it uses the same lists now. Results are
  byte-identical.

## Particle–wall collisions (2ff0cd8)

### Turning walls on

A slip-wall or no-slip-wall section of the `.bc` file with `PARTWALL 1` makes that patch
a particle wall:

```
# BC_SLIPW
PATCH     1 2
NAME      SlipWall
PARTWALL  1     ! particles collide with these patches (default 0)
#
```

Patches without the key, and periodic patches, are not particle walls, so periodic runs
are unchanged. At start-up every rank gets every face of these patches, with
quadrilaterals split into two triangles. Their normals come from Rocflu and point into
the fluid. At start-up the code prints how many triangles came from which patches, and
warns if `COLLISIONS` is 0 (the walls then have no effect).

`filein.vtk` is no longer read, and a run that finds one in its directory stops with a
message. Its triangles had no known fluid side (in the cone file, for example, half of
them point each way), it had to be kept consistent with the mesh by hand, and walls that
the gas sees are exactly the fluid mesh walls. A deck that used it needs `PARTWALL 1` on
its wall patches instead.

### What changed

- **Same force law as pairs.** The pair contact force is now one routine
  (`ppiclf_user_ContactForce`) used by pairs and walls. The wall is a partner of infinite
  mass and radius at rest: m_eff = m, R\* = r, same E, ν, EREST, KSP, friction and rolling
  torque. Pair results are byte-identical to before. The old wall branch had its own Hertz
  constant (√2 too soft for a wall of the same material), no ln² e term in k₁, no
  friction, and started the contact 0.05·dp before the particle touched the wall.
- **Closest point on each triangle.** The old search projected the particle onto each
  triangle's plane and accepted the triangle if the projection was inside it (with 10%
  slack). Near the edge between two triangles both accepted it, and the particle was
  pushed twice. The contact now uses the closest point on each triangle (face, edge or
  vertex). Repeated contacts are dropped, so a flat wall made of many triangles acts as
  one surface, and edges and corners between walls are handled.
- **A fluid side.** A patch wall pushes a particle back towards the fluid even when the
  particle's centre has crossed the wall. The old wall turned its normal towards whichever
  side the centre was on, so once the centre crossed, the wall pushed it out of the domain.
  This assumes that the solid behind a `PARTWALL` patch is at least `NEIGHBORWIDTH` thick
  (at least 4 diameters), which holds for the boundary of the fluid domain.
- **Search.** Wall triangles are filed into the collision subbins of each rank, so a
  particle only checks the triangles near it. The old code checked every wall for every
  particle.
- **Particles just outside the domain.** ppiclF cut its bins back to the fluid domain's
  extent in every direction. A particle whose centre was briefly outside the domain was in
  no bin: it was dropped on 2+ ranks, and on 1 rank the subbin index went below 0 and the
  run crashed. Bins are now only cut back in periodic directions.

### Results

The probe below fires particles (dp = 100 µm, COLLISIONS 1, EREST 0.7, KSP 10) at the
walls of a closed 8 mm box with gas at rest. The gas CFL sets dt ≈ 1.9e-7 s, so a contact
governed by k₁ lasts about 1.9 µs. The old code reads the walls from a `filein.vtk` of the
box (two triangles per face); the new code uses the box's slip walls with `PARTWALL 1`.

![Particle-wall probe, old vs new](figures/fig_walls.png)

- **(a)** The new code rebounds with e = 0.697 on faces, on the diagonal between two
  triangles, at a box edge and in a corner, the same value as a head-on pair. At 1 m/s
  the Hertz stiffness is the softer one and e = 0.724. The old code gave 0.60 on the
  diagonal and 0.67 in the corner, where triangles pushed twice.
- **(b)** At 100 m/s the centre reaches the wall plane (penetration 1.0 r). At 300 and
  1000 m/s the old wall let the particle through and the run crashed; the new wall pushes
  it back from 3.1 r and 10 r behind the wall, with the same e.

![Particle-wall animation, old vs new](figures/anim_walls.gif)

Each panel looks along the wall. The vertical axis is the distance of the particle's
centre from the wall in radii: the particle touches at 1 and its centre has crossed the
wall below 0. Old (orange) and new (blue) particles are drawn to scale side by side, and
the strip below shows the velocity away from the wall over the approach speed. One
frame is one time step.

**How deep particles go.** With the k₁ stiffness the deepest penetration is about
v·KSP·dt/π, whatever the material: about 3.7 r for a 100 µm particle at 300 m/s with
KSP 10 and dt = 1.9e-7 s. The new wall recovers such particles, but the contact is not
resolved in any physical sense. Keeping v·KSP·dt/π well below the radius needs a smaller
KSP or a smaller time step. Rocflu's time step does not take particle speeds into account.

The test case `box/walls` (np 1, 2, 4) runs eight such shots with COLLISIONS 2 (friction
and spin). Each particle that hits one wall must rebound with e within 0.035 of EREST,
every particle must end inside the box, and 2 and 4 ranks must match 1 rank. The old code
fails it: the 300 m/s particle leaves the domain and the run crashes. At an edge or a
corner each wall's friction also slows the motion towards the other walls (e ≈ 0.55 at
the edge, 0.47 in the corner with μ = 0.4), so there the test only checks the rebound.

## Notes for the discussion with the MFC team

- **The normal force turns attractive at the end of a contact.** The normal force is
  −k δ − η v_n. While the particles separate, v_n changes sign, and near the end of the
  contact the dashpot term exceeds the spring term. The animation of the binary
  collisions shows a pull of 1–2% of the rebound speed. MFC's `m_collisions` has the same
  form. The usual fix allows only repulsive normal forces (clip at 0, as in LIGGGHTS).
  Clipping removes some damping, so the actual e comes out above EREST unless η is
  recalibrated (linear spring–dashpot, Schwager & Pöschel 2007):

  | EREST | e with clipping |
  |---|---|
  | 0.5 | 0.550 |
  | 0.7 | 0.718 |
  | 0.9 | 0.902 |

  The friction and rolling torques use |F_n|, so they also act during the attractive
  phase. The code keeps the unclipped force for now, as in MFC.
- **min(k₁, k₂).** MFC uses only the linear spring–dashpot (k₁). ppiclF takes the softer
  of k₁ and the Hertz stiffness. In the Hertz branch e depends on the impact speed (0.725
  for an equal pair at 2 m/s).
- **Material constants.** E = 1 GPa and ν = 0.35 are fixed in the code for every particle
  and wall. Aluminium is about 70 GPa.
- **Resolution.** KSP·dt sets the contact time, so fast impacts penetrate by
  v·KSP·dt/π (see above). A time-step limit from the particle speeds would bound it.

## MPI and periodic boundaries

These fixes have no figure. Each is checked by a regression case that compares 1, 2 and
4 ranks and checks the conservation budgets.

- **Stale virtual cells at rank borders (e3ed370).** `PICL_TEMP_InitSolver` converted the
  conserved variables before virtual cells had φ, and did not exchange them again. A
  quiescent particle cloud set the gas in motion, and results differed by 1% between 1 and
  2 ranks. φ and k_sg are now exchanged before the virtual-cell state is computed in every
  stage.
- **Periodic and symmetry patches (11dc640).** Virtual cells used their own gradients, so
  the fluxes through these patches did not cancel. Gradients are now exchanged
  (`RFLU_MPI_ExchangeGradCells`) and transformed (`RFLU_RELP_TransformGradWrapper`).
- **Uninitialized reads (b10c5f1).** k_sg was missing in the `VELTEMP` inflow branch, which
  gave NaN on 2+ ranks. Also fixed: PICL arrays zeroed, NSCBC patch `dv` set at
  initialization, and an out-of-bounds read of `dv(4)`.

## Still open

Next physics step:
- **Sensitivity of k_sg to the projection width.** The width is about 0.46h and depends on
  the mesh. The options are a fixed physical width or deposit-and-diffuse
  (Capecelatro & Desjardins).

Known limitations:
- Walls that only the particles see (a screen or an obstacle that is not in the fluid
  mesh) are not supported since `filein.vtk` was removed.
- In the closed box of `box/walls` the gas seen by the particles differs between 1 and 2
  ranks by about 0.2% of a 4e-4 m/s flow (driven by φ with `FEEDBACK 0`). The particles
  agree; the gas difference is not explained yet.
- Ghost particles are created and sent in every RK stage, even when nothing uses them.
- The subbin start index in `ppiclf_user_SetYdot` is declared for every subbin even when
  no neighbour search runs. Subbins are `NEIGHBORWIDTH` (at least 4 dp) wide, so very small
  particles with `NEIGHBORWIDTH 0` give a huge count (about 40 GB of address space for
  1 µm tracers in a 4 cm tube; it is never touched).
- The quasi-steady fluctuation mean (`QSFLUCT`) only sees neighbours within
  `NEIGHBORWIDTH`, not the filter width it is meant to average over.
- Feedback lags one RK stage, because particles are projected before `ppiclf_user_SetYdot`.
  It conserves, but the coupling is first order in time.
- Restart from a PICL dump is not verified: the files hold ρ/(1−φ), and the particle VTU
  files that restarts read are single precision.
- `FEEDBACK 0` (one-way) still changes the gas through φ in the pressure.
- PICL supports only one region per MPI rank, and there is no check for this.

## Reproducing the figures

Each figure compares the commit before a fix with the commit after it. The run and plot
scripts (`sweeps.py`, `make_plots.py`, `qs_driver.F`) are kept locally under the
git-ignored `tests/_work/plots/`. They need a built tree of each old commit
(aa3c0b6, 6d00ca4, 1a954d7, eb59a6c) as well as the current one.

The collision figures use the same folder: `animate_collisions.py` animates a run of
`box/binary-collisions` written every 0.1 µs; `walls.py` sets up and runs the wall probe
(`WALLS=vtk` or `WALLS=patch`, `SHOTS=` a subset) with any built tree, and
`animate_walls.py` and `wall_summary.py` draw the animation and the summary figure from
the old (3afd23a) and new (2ff0cd8) probe runs. `shock_interp.py` (planar shock) and
`sheet.py` (particle sheet) run and plot the interpolation comparisons in shock-driven
flow, with aa3c0b6 as the old interpolation.
