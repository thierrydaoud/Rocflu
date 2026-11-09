# Lax Shock Tube

A Riemann problem for the 1D Euler equations of an ideal gas, like [Sod](sod.md), but harder: the gas on the left is already moving, the pressure ratio is larger, and the contact separates a light gas from one almost four times denser. Schemes that pass Sod often show overshoots or a smeared contact here.

Reference: P. D. Lax, "Weak solutions of nonlinear hyperbolic equations and their numerical computation", *Comm. Pure Appl. Math.* 7:159–193 (1954).

## Problem

Equations: the 1D Euler equations for an ideal gas, with no viscosity, heat conduction or source terms.

| Quantity | Value |
|---|---|
| Ratio of specific heats γ | 1.4 |
| Domain | 0 ≤ x ≤ 1 |
| Initial jump at | x₀ = 0.5 |
| End time | t = 0.14 |

Initial state (density ρ, velocity u, pressure p):

| | ρ | u | p |
|---|---|---|---|
| Left, x < 0.5 | 0.445 | 0.698 | 3.528 |
| Right, x > 0.5 | 0.5 | 0 | 0.571 |

All quantities are dimensionless (see [Running in physical units](#running-in-physical-units)). The literature uses several end times (0.13, 0.14, 0.16) and domains; this page fixes its own.

### Boundaries

No wave reaches either end by t = 0.14 (see [Wave positions](#wave-positions-at-t--014)). The rarefaction head, at x = 0.131, comes closest.

- **Left:** the gas flows in at Mach 0.21, so the boundary must keep supplying the left state. Use an inflow that holds it (for example fixed velocity and temperature, pressure from inside) or a zero-gradient boundary. A wall is not acceptable here.
- **Right:** the gas is at rest; a zero-gradient outflow or a slip wall.

### Multi-dimensional solvers

As for Sod: a slab with a few cells (we use 8) in y and z, periodic or slip-wall sides. The exact solution has no transverse velocity.

## Exact solution

| Quantity | Value |
|---|---|
| p\* | 2.46610 |
| u\* | 1.52872 |
| ρ\* left of the contact | 0.34457 |
| ρ\* right of the contact | 1.30408 |
| Shock speed S | 2.47932 |
| Shock Mach number (into the right state) | 1.961 |

Computed as for Sod (Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., ch. 4), with the same two independent implementations agreeing to all digits shown.

### Wave positions at t = 0.14

| Wave | Speed | Position |
|---|---|---|
| Rarefaction head | −2.63357 (= u_L − a_L) | 0.13130 |
| Rarefaction tail | −1.63669 (= u\* − a\*_L) | 0.27086 |
| Contact | 1.52872 (= u\*) | 0.71402 |
| Shock | 2.47932 | 0.84710 |

Sound speeds a = √(γp/ρ): a_L = 3.33157, a\*_L = 3.16541.

Inside the rarefaction (head < x < tail), with ξ = (x − x₀)/t:

- u = 2/(γ+1) · (a_L + (γ−1)/2 · u_L + ξ)
- a = a_L − (γ−1)/2 · (u − u_L)
- ρ = ρ_L (a/a_L)^(2/(γ−1))
- p = p_L (a/a_L)^(2γ/(γ−1))

## Grids

Uniform cells, Δx = 1/n with n = 200, 400 and 800.

## What to measure

The same measures as [Sod](sod.md#what-to-measure), with the shock position where ρ crosses (ρ\*_R + ρ_R)/2 = 0.90204, compared with 0.84710.

For conservation, the boundary fluxes are not zero here. Per unit area and unit time, with E = p/(γ−1) + ρu²/2:

| | Left (in) | Right (out) | Net rate |
|---|---|---|---|
| Mass ρu | 0.310610 | 0 | 0.310610 |
| Momentum ρu² + p | 3.744806 | 0.571 | 3.173806 |
| Energy u(E + p) | 8.694569 | 0 | 8.694569 |

So the totals at t should equal the initial totals plus t times the net rate.

## Running in physical units

With P = 10⁵ Pa, R = 1 kg/m³ and a 1 m domain, the velocity scale is 316.228 m/s: the left gas moves at 220.727 m/s and the end time is 0.14/316.228 = 4.4272e-4 s.
