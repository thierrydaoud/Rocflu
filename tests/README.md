# Rocflu + ppiclF regression tests

Small cases (a few hundred cells, 1–40 particles) that run in seconds each. Every
case is compared against a stored golden, and most cases also have physics checks
that judge whether the result is correct, not just unchanged.

## Running

Build the code first (`makeall.sh`). Then, from this directory:

```bash
module load gcc/14.2.0 openmpi/5.0.7
./run_tests.py list                        # show cases
./run_tests.py run                         # all cases, compare to goldens (~40 s)
./run_tests.py run -o shocktube/two-way    # cases whose name contains a pattern
./run_tests.py run --exe /path/to/tree     # test executables from another tree
./run_tests.py run --generate              # rewrite goldens from this build
```

Pre-processing (`rflumap`, `rflupart`, `rfluinit`) runs serially where you launch the
script. The solver runs through `srun`:
- From a login node, each case is its own short job (`-A bala1s -p hpg-dev`). Override
  with `--account`/`--partition` or `ROCPICL_TEST_ACCOUNT`/`ROCPICL_TEST_PARTITION`.
- Inside an allocation (`salloc`), cases run as job steps.
- `--local-ranks N` runs cases with at most N ranks with `mpirun --bind-to none` where
  you launch the script (e.g. `-j 6 --local-ranks 4` on a login node, at most 16 cores;
  the whole suite takes ~5 min). Without `--bind-to none`, concurrent cases all pin
  to the first cores and run ~10x slower.

A case whose solver run fails is set up again and retried once; `srun` on a shared
partition can fail when the Slurm controller is unresponsive. The retry and any Slurm
messages are shown next to the case.

Run directories go to `_work/runs/<case>` (git-ignored). They are deleted for cases
that pass, unless you pass `--keep`.

## Results

| Status | Meaning |
|---|---|
| `PASS` | Matches the golden, and the physics checks pass |
| `FAIL` | The run failed, the result differs from the golden, or a physics check failed |
| `XFAIL` | The case has a `known_issue`. Its physics checks fail as expected, and it still matches its golden |
| `XPASS` | The physics checks of a `known_issue` case now pass. Remove the marker |

Tolerances:
- Gas fields are double precision: relative 1e-10 of each field's maximum.
- Particle files (`par*.vtu`) are single precision: relative 2e-6.

## What a case is

Cases are defined in `rocpicl_tests/cases.py`. Each one specifies:
- a box grid made with `mesh_generator`;
- `templates/rectshktb.inp` with per-case overrides, written as `{("SECTION", "KEY"): value}`;
- a boundary layout:
  - `shocktube`: `templates/rectshktb.bc`, with inflow at x−, outflow at x+, and periodic y and z;
  - `periodic`: periodic on all six faces;
  - `closed`: slip walls on all six faces;
  - `partwall`: slip walls on all six faces that are also particle walls (`PARTWALL 1`);
- a particle list, written to `points.dat` (`AL x y z dp`, optionally followed by an
  initial velocity `u v w`; particles without one start at rest);
- a rank count.

Particles never start on a cell face, so round-off can't decide which cell or rank owns
them. A case's id changes whenever its definition changes. A stale golden is reported
rather than silently compared.

Particles are matched by their line number in `points.dat`, which is recovered from the
ppiclF tags. That's what lets `same_as` compare a 2-rank run with a 1-rank run particle
by particle.

## Physics checks (`rocpicl_tests/checks.py`)

| Check | What it requires |
|---|---|
| `particles_moved` | Something happened, so the other checks aren't vacuous |
| `mass`, `gas_energy` | Gas totals are conserved (closed or periodic box, no inflow) |
| `momentum`, `energy` | The gas loses what the particles gain (periodic box); energy includes particle kinetic and thermal energy |
| `particles_heat_exchange` | Particles at rest in gas at another temperature exchange heat |
| `quiescent` | Gas at rest stays at rest (fixed particles) |
| `linear_interp` | Particles read a linear gas density exactly (`INITFLOW IVAL1 = 1`, density `RVAL2*(1 + RVAL7*x + RVAL8*y + RVAL9*z)` at rest) |
| `projection_shape` | One particle's volume fraction is centred on it (within 0.06 cell) and spread 0.35-0.50 cell lengths in each direction |
| `collision_pairs` | Overlapping pairs pushed apart by collisions alone conserve each pair's momentum, separate, and all move like the first pair, whether the partner is real, a ghost from another rank, or a periodic image |
| `binary_collisions` | Colliding pairs (ids 2k−1, 2k) conserve momentum (1e-5) and angular momentum about their centre of mass (0.5%), and separate; head-on pairs rebound with the input `EREST` (within 1%) |
| `wall_collisions` | Every particle that reaches a wall of the box rebounds from it and ends inside the box, clear of the wall; one that hits a single wall rebounds with the input `EREST` (within 0.035) |
| `same_as` | Particles match the same case on a different number of ranks (all particle fields, or only those named by the case's `same_as_fields`) |

## Known issues tracked by the suite

None at the moment. Fixed and now guarded:
- non-conservative periodic patches (`box/gas-jump/periodic`, 11dc640);
- stale gas state in virtual cells at region borders with particles, which set the gas
  moving around fixed particles and made 1- and 2-rank runs differ by ~1%
  (`box/quiescent-cloud/np2`, `np4`, and the `same_as` check of
  `box/momentum-exchange/np2`, e3ed370);
- resting particles skipping heat transfer (`box/resting-heat-exchange`, f0ae123);
- inverse-distance interpolation, which was not exact for a linear field (errors up to
  ~20% of the change across a cell, depending on where the particle sat), replaced by
  weighted linear least squares (`box/linear-interp/np1`, `np2`, `np4`, 89f9661);
- projection kernel sized by V^(1/3), which hardly spread particles along the long side
  of elongated cells (414ce7f), and cells just outside the bin domain left out of the
  stencil of particles at the edge of a cloud, so one particle projected onto 8 cells
  instead of 27 (affad85) (`box/elongated-projection/np1`, `np2`).

- collision stiffness from particle i's own mass, so unequal pairs did not conserve
  momentum (an 8:1 pair at 20 m/s ended with its momentum reversed) and equal pairs
  were twice as stiff as `KSP` intends (`box/binary-collisions/np1`, `np2`, `np4`,
  3afd23a). With the old stiffness the case fails momentum on both unequal pairs and
  restitution on three pairs.
- particle walls (`box/walls/np1`, `np2`, `np4`, 2ff0cd8): before, walls came only from
  `filein.vtk`, a particle over the edge between two wall triangles was pushed twice, a
  fast particle was pushed through the wall, and a particle whose centre was outside the
  domain was dropped or crashed the run. The old code crashes on this case (the 300 m/s
  particle leaves the domain). `same_as` compares only the particle state (`y`): the gas
  in this closed box differs slightly between decompositions, so the gas seen by the
  particles does too (see `docs/PseudoTurb_fixes.md`, Still open).

Guards with no known bug behind them:
- ghost particles (`box/ghost-collisions/np1`, `np2`, `np4`). Ghosts must carry the
  partner's velocity, not only its position and properties: a build whose ghosts lose
  `y(4:LRS)` fails `collision_pairs` on every rank count (periodic pairs separate 10% too
  fast) and `same_as` on 2 and 4 ranks. ppiclF splits bins in z on 2 ranks and in y and z
  on 4 ranks; the x direction only uses periodic ghosts.

## Flow-file conventions

What the five arrays of `SOL_*/<case>.mixt.cva_*` hold depends on whether particles are
on (`RFLU_TimeStepping` converts real cells just before each dump):

- PICL off: conserved bulk state `rho, rho*u, rho*v, rho*w, rho*E`
- PICL on: gas-phase primitive state `rho_g, u, v, w, E`, with `rho_g = rho/(1-phi)`

Goldens store the arrays as written. Conservation totals convert to bulk quantities,
taking `phi` from the "Particle Volume Fraction" field in `PARAVIEW_*/<case>_<t>_<reg>.vtu`.

## Adding a case

1. Add it to `rocpicl_tests/cases.py`.
2. Run `./run_tests.py run --generate -o <name>`. The golden is written only if the
   physics checks pass, or if the case has a `known_issue`.
3. Commit the new `goldens/<name>/` directory.
