# Pass Criteria

How a run of a [problem](../problems/) is judged. Each criterion is one of four kinds:

- **Absolute:** must hold for any correct solver.
- **Theory:** follows from the known convergence behaviour of the scheme's order.
- **Reference:** compared with a simple scheme whose behaviour is understood, run at the same grid and time step.
- **Consistency:** the same problem set up two ways must agree.

Limits are fixed before a rung is run, and not loosened afterwards to make a scheme pass. Some quantities are reported without a limit, because there is no defensible threshold for them.

Notation: n cells in x, Δx = 1/n; L1(q) = Σ|q_i − q_exact(x_i)| Δx; the rate between n and 2n is log₂(L1_n / L1_2n), taken between the two finest grids.

## Rung 1: Sod and Lax

| Kind | Criterion | Limit | Why |
|---|---|---|---|
| Absolute | Run reaches the end time | — | |
| Absolute | Conservation: \|total − (initial + t · net boundary flux)\| / \|total\|, for mass, momentum and energy | < 1e-10 when both end cells hold their initial state to 1e-12; otherwise reported | Round-off only, for a conservative scheme |
| Absolute | Rank independence: largest relative difference in ρ, u, p between 1 and N parallel ranks | < 1e-10 | Round-off only |
| Absolute | Transverse velocity (2D, 3D): rms of v, w divided by the largest initial sound speed | < 1e-10 | The exact solution has none |
| Theory | L1(ρ) rate, first order in space | ≥ 0.5 | Contacts converge at p/(p+1) in L1 (Banks, Aslam & Rider 2008), the slowest wave |
| Theory | L1(ρ) rate, second order in space | ≥ 0.67 | Same, p = 2 |
| Theory | Shock position error, in cells: \|x_s − x_s,exact\| / Δx | ≤ 1.5 on every grid | A captured shock converges at first order, so its error stays a fixed number of cells |
| Reference | First-order L1(ρ) divided by the exact-Riemann Godunov scheme's L1(ρ) at the same n and time step | ≤ 1.5 | Godunov is the classical first-order upwind reference ([`exact/godunov.py`](../exact/godunov.py)) |
| Consistency | \|L1_2D − L1_3D\| / L1_3D, same n | < 0.01 | The slab is 1D; the dimension must not change the answer |

Reported, no limit:

- contact width (cells between 10% and 90% of the jump), overshoot at the shock and contact;
- number of time steps and run time;
- L1 of u and p and their rates.

### Notes

- **Time step for the reference.** Solvers define the CFL number differently (some sum the wave speeds over all faces of a cell, so a "CFL 0.5" in 3D is about 0.09 in the 1D sense). The Godunov reference is therefore run with the same number of time steps as the solver took, as equal steps (`godunov.py <problem> <n> --steps N`), not with the solver's CFL number. `godunov.py` stops if that makes its 1D CFL number exceed 1.
- **Conservation needs undisturbed ends.** The check uses the exact boundary fluxes, so it only holds while the end cells keep their initial state. At first order the rarefaction is smeared far ahead of its exact head and can reach the boundary even though the exact wave does not. On Lax at first order this happens on every grid used here: at n = 400 the first cell's pressure is off by 5e-6 to 1e-4 depending on the scheme, and the totals by 3e-8 to 2e-6, both falling fast with n; at second order the first cell is exact to round-off and the totals are 1e-15. Such runs are reported, not failed.
- **Why 1e-10.** Measured round-off for correct runs is 1e-16 to 1e-12 (double precision, up to a few thousand steps). The limit leaves two orders of margin and still catches the smallest real faults we have seen (1e-10 to 1e-9 from a rank border).
