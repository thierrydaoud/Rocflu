# Rocflu + ppiclF fixes on the `PseudoTurb` branch

The goal of this branch is to validate the subgrid pseudo-turbulence (PT) model. Before
the PT results can be trusted, the particle–gas coupling has to be correct: interpolating
the gas to the particles, projecting the particles back onto the mesh, the drag and heat
closures, and conservation across MPI ranks and periodic boundaries. This page lists what
was wrong, how each problem was fixed, and how the fix was checked. The PT model itself
(the "k_sg gets too high" issue) is not addressed yet; see [Still open](#still-open).

All changes are covered by the regression suite in [`tests/`](../tests/README.md):
26 cases, compared against stored goldens and checked against physics on 1, 2 and 4 ranks.

## Summary

| Area | Problem | Fix | Commit |
|---|---|---|---|
| Interpolation | Inverse distance (1/d³) had errors of up to ±18% of the change in density across one cell, even in a linear field, and was discontinuous at cell faces | Weighted linear least squares, exact for linear fields | 89f9661, 6d00ca4 |
| Projection | The kernel used one length for every direction, so the spread on stretched cells was wrong | Each direction is scaled by its own cell length, as in MFC | 414ce7f |
| Projection | Cells just outside the bin domain were dropped, so cloud-edge particles reached 8 cells instead of 27; φ was ~8% too high at the edge | Map those cells to particle bins | affad85 |
| Closures | Forces and heat transfer were skipped below 1 mm/s slip, and the old value was kept | Evaluate at any slip; skip only the quasi-steady drag at zero slip | f0ae123 |
| Quasi-steady drag | Gidaspow's C_D was capped below Re = 0.1 (100% error); Osnes had a spurious Knudsen reduction (25% error) | Aligned with MFC; QSFORCE renumbered | 94607fd |
| MPI | Gas state in virtual cells was stale at rank borders with particles, causing spurious gas motion in a quiescent cloud and a 1% difference between 1 and 2 ranks | Convert and re-exchange after φ is known; exchange φ and k_sg every stage | e3ed370 |
| Conservation | Fluxes on periodic and symmetry patches were not conservative | Exchange and transform gradients in virtual cells | 11dc640 |
| Robustness | Uninitialized reads gave NaN on 2+ ranks; a single particle gave a filter width of 0 | Initialize; filter width from cell–box overlap | b10c5f1, 48f64b7 |
| Build | Required `SPEC=1` | Builds without it; stops if species are requested | aa3c0b6 |

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
