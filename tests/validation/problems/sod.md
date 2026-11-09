# Sod Shock Tube

A Riemann problem for the 1D Euler equations of an ideal gas. The solution has all three wave types: a rarefaction moving left, then a contact discontinuity and a shock moving right. The exact solution is known, so the error of a numerical solution can be measured directly.

Reference: G. A. Sod, "A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws", *J. Comput. Phys.* 27:1–31 (1978).

## Problem

Equations: the 1D Euler equations for an ideal gas, with no viscosity, heat conduction or source terms.

| Quantity | Value |
|---|---|
| Ratio of specific heats γ | 1.4 |
| Domain | 0 ≤ x ≤ 1 |
| Initial jump at | x₀ = 0.5 |
| End time | t = 0.2 |

Initial state (density ρ, velocity u, pressure p):

| | ρ | u | p |
|---|---|---|---|
| Left, x < 0.5 | 1 | 0 | 1 |
| Right, x > 0.5 | 0.125 | 0 | 0.1 |

All quantities are dimensionless. To run in physical units, choose a pressure scale P and a density scale R; velocity then scales with √(P/R) and time with L√(R/P) for a domain of length L (see [Running in physical units](#running-in-physical-units)).

### Boundaries

No wave reaches either end by t = 0.2 (see [Wave positions](#wave-positions-at-t--02)), so the gas next to each boundary stays at its initial state. Any boundary condition that keeps a uniform state at rest unchanged is acceptable: a zero-gradient (transmissive) outflow or a slip wall at both ends.

### Multi-dimensional solvers

Run the problem as a slab: x as above, and a few cells (we use 8) in y and, in 3D, z, with periodic or slip-wall boundaries on those sides. The exact solution has no y or z dependence and no transverse velocity, which gives an extra check (see [What to measure](#what-to-measure)).

## Exact solution

The solution depends only on x/t. Between the waves there are two "star" states of equal pressure p\* and velocity u\*, separated by the contact:

| Quantity | Value |
|---|---|
| p\* | 0.30313 |
| u\* | 0.92745 |
| ρ\* left of the contact | 0.42632 |
| ρ\* right of the contact | 0.26557 |
| Shock speed S | 1.75216 |
| Shock Mach number (into the right state) | 1.656 |

These values come from solving the standard pressure equation (E. F. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Springer 2009, ch. 4). Two independent implementations (`exact/riemann.py` and a bisection solver) agree to all digits shown.

### Wave positions at t = 0.2

| Wave | Speed | Position |
|---|---|---|
| Rarefaction head | −1.18322 (= −a_L) | 0.26336 |
| Rarefaction tail | −0.07027 (= u\* − a\*_L) | 0.48595 |
| Contact | 0.92745 (= u\*) | 0.68549 |
| Shock | 1.75216 | 0.85043 |

Here a = √(γp/ρ) is the sound speed: a_L = 1.18322 in the left state, a\*_L = 0.99772 in the star state left of the contact.

Inside the rarefaction (head < x < tail), with ξ = (x − x₀)/t and u_L = 0:

- u = 2/(γ+1) · (a_L + ξ)
- a = a_L − (γ−1)/2 · u
- ρ = ρ_L (a/a_L)^(2/(γ−1))
- p = p_L (a/a_L)^(2γ/(γ−1))

`exact/riemann.py` evaluates the full solution at any x and t.

## Grids

Uniform cells, Δx = 1/n with n = 200, 400 and 800 cells in x. Three grids give two convergence rates.

## What to measure

Compare the cell averages, or cell-centre values, at t = 0.2 with the exact solution at the cell centres. Definitions and pass limits are in `metrics/criteria.md`.

| Measure | Definition |
|---|---|
| L1 error of ρ, u, p | Σ \|q_i − q_exact(x_i)\| Δx |
| Convergence rate | log₂(L1 at n / L1 at 2n) |
| Shock position | where ρ crosses (ρ\*_R + ρ_R)/2, interpolated linearly between cells; compare with 0.85043 |
| Conservation | total mass, momentum and energy at t, compared with the initial totals plus the boundary fluxes (zero for mass and energy here, p_L − p_R per unit area for momentum per unit time) |
| Transverse velocity (2D, 3D) | rms of v and w over all cells; zero exactly |
| Dimension consistency (2D, 3D) | difference of the L1 errors between a 2D and a 3D run on the same n |

## Running in physical units

The cases in this repository use P = 10⁵ Pa and R = 1 kg/m³ on a 1 m domain. The velocity scale is √(P/R) = 316.228 m/s and the time scale 1/316.228 s, so the end time is 0.2/316.228 = 6.3246e-4 s. With γ = 1.4 the specific heat at constant pressure only sets the temperature, not the flow.
