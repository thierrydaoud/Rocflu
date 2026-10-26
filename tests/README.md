# Rocflu Regression Tests

Small cases (a few hundred cells, 1–40 particles) that run in seconds each. Every
case is compared against a stored golden, and most cases also have physics checks
that judge whether the result is correct, not just unchanged.

This suite tests the coupled code end to end: grid, pre-processing, the Rocflu solver
and ppiclF together, on 1, 2 and 4 MPI ranks. It replaces ppiclF's own unit tests
(`libpicl/unit_tests`, `short_tests` and the `#ifdef TEST` code), which were removed
in 3c41db2.

There are two kinds of test:
- **Golden comparison.** The final gas and particle state must match the one stored
  in `goldens/<case>/golden.json`. This catches any change in results, intended or not.
- **Physics checks.** Conservation, exact interpolation, restitution and so on, tested
  against the physics rather than against an older run. They stay meaningful when a fix
  changes every golden on purpose.

Particle fields are read by name (`vp_x`, `rho_g`, `dp`, ...; see
`libpicl/user_files/ppiclf_user_OutputFields.f`). The template input
`templates/rectshktb.inp` turns on `PARTICLE_WRT`, `GAS_WRT` and `GRADGAS_WRT`, because the
goldens and checks use the particle temperature, spin, density and volume and the gas
values at the particles; keep them on. Gas gradients are written only when a model that is
on needs them (`ppiclf_user_Layout`), so cases without pressure-gradient, added-mass,
viscous-unsteady or rotational-collision models have no `dpdx_g` and similar fields. A
field that is in a golden but no longer written fails the comparison.

The 38 cases fall into three groups: [shock tube](#shock-tube-cases) (force models,
one- and two-way coupling), [box](#box-cases) (conservation, interpolation, projection,
collisions, walls) and [gas only](#gas-only-cases) (the flow solver without particles).

## Running

### 1. Build

Build the code in the tree you want to test: `./rocflu.sh build` (see
[Getting Started](../docs/GettingStarted.md)). The suite runs the executables `rflumap`,
`rflupart`, `rfluinit` and `rflump` from `build/bin` (the last configuration built), or
from the directory given with `--exe`. `./rocflu.sh test [options]` runs
`./run_tests.py run` with the modules loaded and passes the options on.

### 2. Run

From this directory:

```bash
module purge; module load gcc/14.2.0 openmpi/5.0.7
./run_tests.py list                        # show cases
./run_tests.py run                         # all cases, compare to goldens (~40 s)
./run_tests.py run -o shocktube/two-way    # cases whose name contains a pattern
./run_tests.py run -o walls np1            # several patterns: any of them
./run_tests.py run --exe /path/to/tree     # test executables from another tree
./run_tests.py run -j 21                   # cases run at the same time (default 8)
./run_tests.py run --keep                  # keep run directories of passing cases
```

Load the modules in the same shell. Without them every case fails with a `GLIBCXX`
error.

Pre-processing (`rflumap`, `rflupart`, `rfluinit`) runs serially where you launch the
script. The solver runs through `srun`:
- From a login node, each case is its own short job (`-A bala1s -p hpg-dev`). Override
  with `--account`/`--partition` or `ROCPICL_TEST_ACCOUNT`/`ROCPICL_TEST_PARTITION`.
  `-j 21` runs the whole suite in about 40 s when the partition is free.
- Inside an allocation (`salloc`), cases run as job steps.
- `--local-ranks N` runs cases with at most N ranks with `mpirun --bind-to none` where
  you launch the script (e.g. `-j 4 --local-ranks 4` on a login node, at most 16 cores;
  the whole suite takes ~5 min). Without `--bind-to none`, concurrent cases all pin
  to the first cores and run ~10x slower.

Each case's input file goes through the case validator
(`utilities/scripts/case_validator.py`, see docs/Running.md) before it runs; a case
it rejects fails with the validator's messages. The validator has its own tests,
which need no cluster: `python3 tests/test_case_validator.py`.

A case whose solver run fails is set up again and retried once; `srun` on a shared
partition can fail when the Slurm controller is unresponsive. The retry and any Slurm
messages are shown next to the case.

### 3. Read the Result

The script prints one line per case, a total (e.g. `38 PASS`) and a closing line
(`✔ All 38 cases passed in 29 s`, or `✗ 2 of 38 cases failed`). On a terminal the statuses
are coloured (PASS green, FAIL red, XFAIL yellow, XPASS cyan); `NO_COLOR=1` turns that off. See [Results](#results)
for what each status means.

Run directories go to `_work/runs/<case>` (git-ignored). They are deleted for cases
that pass, unless you pass `--keep`. A failed case keeps its directory, with the
input files and one log per step: `mesh.log`, `rflumap.log`, `rflupart.log`,
`rfluinit.log`, `rflump.log`. A golden mismatch names the worst particle or cell and
its value against the golden.

### 4. After Changing the Code

- **Results should not change** (a clean-up, a speed-up): run the suite. It must pass
  as it is.
- **Results change on purpose** (a fix): the goldens of the affected cases have to be
  rewritten. Do it from a clean build of a committed tree, so each golden records which
  commit produced it:
  1. Commit the fix.
  2. Do a clean build of that commit, ideally in a separate worktree, so the binaries
     cannot come from uncommitted edits.
  3. `./run_tests.py run --generate --exe <tree>`. A golden is written only if the
     case's physics checks pass. `--force` writes it anyway; use it only for a case
     marked `known_issue`.
  4. Run the comparison twice (`./run_tests.py run --exe <tree>`). Both must pass.
  5. Commit `goldens/` separately from the fix.

  Never regenerate goldens to make a failing physics check go away: the check is
  telling you the new result is wrong.

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

## What a Case Is

Cases are defined in `rocpicl_tests/cases.py`. Each one specifies:
- a box grid made with `mesh_generator`;
- `templates/rectshktb.inp` with per-case overrides, written as `{("SECTION", "KEY"): value}`;
- a boundary layout:
  - `shocktube`: `templates/rectshktb.bc`, with inflow at x−, outflow at x+, and periodic y and z;
  - `periodic`: periodic on all six faces;
  - `closed`: slip walls on all six faces;
  - `partwall`: slip walls on all six faces that are also particle walls (`PARTWALL 1`);
- a particle list `(x, y, z, dp)` or `(x, y, z, dp, u, v, w)` (without a velocity the
  particle starts at rest), written to `points.dat` as `x y z vx vy vz dp`; density and
  specific heat come from `DENSITY` and `SPECIFICHEAT` in `# PICL`;
- a rank count.

Particles never start on a cell face, so round-off can't decide which cell or rank owns
them. A case's id changes whenever its definition changes. A stale golden is reported
rather than silently compared.

Particles are matched by their line number in `points.dat`, which is recovered from the
ppiclF tags. That's what lets `same_as` compare a 2-rank run with a 1-rank run particle
by particle.

## Physics Checks (`rocpicl_tests/checks.py`)

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
| `fluctuations_active` | Every particle feels a nonzero quasi-steady force fluctuation at the end (`Ffluct_*`, written with `MODEL_WRT 1`), so the noise is really tested |
| `wall_collisions` | Every particle that reaches a wall of the box rebounds from it and ends inside the box, clear of the wall; one that hits a single wall rebounds with the input `EREST` (within 0.035) |
| `same_as` | Particles match the same case on a different number of ranks (all particle fields, or only those named by the case's `same_as_fields`) |

The state a case is judged by (goldens, physics checks) comes from its dumps in
`restart_data/` (`rocpicl_tests/state.py`), in the form of the old-style flow files: with
particles `rho_g, u, v, w, E` and `phi`, without the conserved state.

The template switches the old-style output off (`OLD_OUTPUT_WRT 0`, as for a user's run).
Five cases switch it on (`OLD_OUTPUT_CASES` in `cases.py`: 2 regions and 1, NSCBC outflow
with history, random streams, gas only); for them the harness also checks that each `SOL_*`
flow file matches the fluid dump of the same time (1e-12 relative), that each `par*.vtu`
matches the particle dump of the same time field by field (single precision), and that
rfluinit's step 0 dump exists.

**Post program**: every case also runs `rflupost -s all`; the cases with the old-style
output compare `post/` with it at every step: each fluid field of each region with `PARAVIEW_*` of the
same time (1e-12 relative to the field's largest value; the template turns every
`*_WRT` field on), each particle field and position with the `par*.vtu` of the same step
(single precision). Without particles the old `Velocity` and `Energy` are the momentum
and the energy per volume, and are compared as such.

**Restart test** (cases with `restart_test=True`): after the run, the harness copies the run
directory to `<case>_restart`, sets `START_STEP` to the first dump after step 0, runs rflump
again to the same end time and requires every final dump (fluid of each region, particles)
to be **byte-identical** to the uninterrupted run's. It runs on `box/momentum-exchange/np2`
(two-way coupling), `shocktube/two-way/viscous-unsteady/np2` (history of past steps, NSCBC
outflow), `box/qs-fluctuations/np2` (random streams) and `box/binary-collisions/np2`
(collisions).

## The Test Cases

Every case writes two solution dumps and is compared at the last one. "np*N*" in a name
is the number of MPI ranks. Where a case exists on several rank counts, the 2- and
4-rank versions also run `same_as` against the 1-rank one (except where noted), so
they test that the result does not depend on the decomposition.

### Shock Tube Cases

A 40 × 4 × 4 mm tube of 1 mm cells, periodic in y and z, with a Mach 1.27 shock starting at
x = −10 mm and running into particles at rest (100 µm, quasi-steady drag, Osnes).
They run to 20 µs.

| Case | Ranks | Setup | Checks | What it guards |
|---|---|---|---|---|
| `shocktube/one-way/1-particle` | 1, 2 | One particle, `FEEDBACK 0` | `particles_moved`, `same_as` | ppiclF bins around a single particle (this once gave a filter width of 0) |
| `shocktube/one-way/4-particles` | 1, 2 | Four particles, `FEEDBACK 0` | `particles_moved`, `same_as` | Interpolation and drag with several particles |
| `shocktube/two-way/4-particles` | 1, 2 | Four particles, `FEEDBACK 1` | `particles_moved`, `same_as` | Feedback to the gas, across a rank border |
| `shocktube/two-way/pressure-gradient` | 2 | + `PGFORCE 1` | `particles_moved` | Pressure-gradient force |
| `shocktube/two-way/viscous-unsteady` | 2 | + `VISCOUSUNSTEADY 1` | `particles_moved` | History force |
| `shocktube/two-way/heat-transfer` | 2 | + `HEATTRANSFER 3` | `particles_moved` | Heat-transfer closure |
| `shocktube/two-way/qs-fluctuations` | 2 | + `QSFLUCT 2` | `particles_moved` | Quasi-steady force fluctuations run; the four particles have no neighbours, so the noise amplitude is 0 (see `box/qs-fluctuations`) |
| `shocktube/two-way/qs-parmar` | 2 | `QSFORCE 2` | `particles_moved` | Parmar drag |
| `shocktube/two-way/qs-gidaspow` | 2 | `QSFORCE 1` | `particles_moved` | Gidaspow drag |

The force-model cases have no physics check beyond `particles_moved`. Their golden is
what guards them: any change in the model's result shows up as a golden mismatch.

### Box Cases

An 8 mm cube of 1 mm cells (8 × 8 × 8), either periodic on all faces or closed by slip
walls. The gas is uniform (1.2 kg/m³, 101 325 Pa) unless stated otherwise.

| Case | Ranks | Setup | Checks | What it guards |
|---|---|---|---|---|
| `box/momentum-exchange` | 1, 2 | Periodic; gas at 50 m/s through 40 particles (20 µm) at rest, two-way | `mass`, `momentum`, `energy`, `same_as` | The gas loses exactly the momentum and energy the particles gain |
| `box/qs-fluctuations` | 1, 2 | Periodic; gas at 50 m/s through 40 particles (20 µm) with a ±1 m/s velocity spread, `QSFLUCT 2`, `MODEL_WRT 1`, `NEIGHBORWIDTH 3e-3` | `particles_moved`, `fluctuations_active`, `same_as` | Each particle draws from its own random stream: 1 and 2 ranks agree (96b28a4; one generator per rank failed this) |
| `box/quiescent-cloud` | 2, 4 | Periodic; 30 fixed particles (300 µm, uneven) in gas at rest, `VELOCITYMODEL 0`, `PGFORCE 1` | `quiescent`, `mass` | Stale virtual cells at rank borders set the gas moving (e3ed370) |
| `box/resting-heat-exchange` | 1 | Periodic; particles at rest at 296.4 K in gas at rest at 294.2 K, `HEATTRANSFER 2` | `particles_heat_exchange`, `mass`, `energy` | Heat transfer at zero slip used to be skipped (f0ae123) |
| `box/linear-interp` | 1, 2, 4 | Closed; gas at rest with density linear in x, y, z; 40 tracers (1 µm) | `linear_interp`, `mass` | Interpolation is exact for a linear field (89f9661) |
| `box/elongated-projection` | 1, 2 | Closed, 2:1:1 cells; one fixed particle | `projection_shape`, `mass` | The projection kernel spreads ~0.44 cell lengths in every direction (414ce7f, affad85) |
| `box/ghost-collisions` | 1, 2, 4 | Periodic; 9 overlapping pairs at rest, collisions as the only force | `collision_pairs`, `same_as` | Pairs across bin faces, rank borders, periodic faces, edges and corners all behave like a pair inside one bin (ghost particles) |
| `box/binary-collisions` | 1, 2, 4 | Periodic; 4 approaching pairs (equal and 8:1 mass, head-on and oblique, 1–10 m/s), `COLLISIONS 2` | `binary_collisions`, `same_as` | Momentum and angular momentum conserved, rebound with `EREST` (3afd23a) |
| `box/walls` | 1, 2, 4 | Closed, all walls `PARTWALL 1`; 7 particles hitting faces, a quad diagonal, an edge and a corner at 1–300 m/s | `wall_collisions`, `same_as` (particle state only) | Particle–wall contacts (2ff0cd8) |

All box cases with particles also check `particles_moved` where the particles are
meant to move.

### Gas-Only Cases

No particles (`USED 0`): a pressure and density jump at x = 4.3 mm in the box, to test
the flow solver's own conservation.

| Case | Ranks | Setup | Checks | What it guards |
|---|---|---|---|---|
| `box/gas-jump/closed` | 2 | Slip walls, two regions | `mass`, `gas_energy` | Conservation across a region border |
| `box/gas-jump/periodic` | 1, 2 | Periodic on all faces | `mass`, `gas_energy` | Conservative fluxes across periodic faces (11dc640) |
| `box/gas-jump/symmetry` | 1, 2 | Symmetry on all faces | `mass`, `gas_energy` | Conservative fluxes across symmetry faces (11dc640) |

## Known Issues Tracked by the Suite

None at the moment. Fixed and now guarded:
- non-conservative periodic and symmetry patches (`box/gas-jump/periodic`,
  `box/gas-jump/symmetry`, 11dc640);
- stale gas state in virtual cells at region borders with particles, which set the gas
  moving around fixed particles (1e-4 m/s) and made the total gas energy of 2- and 4-rank
  runs differ from 1 rank by 5e-8
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
  particles does too (see `docs/Improvements.md`, Still open).

Guards with no known bug behind them:
- ghost particles (`box/ghost-collisions/np1`, `np2`, `np4`). Ghosts must carry the
  partner's velocity, not only its position and properties: a build whose ghosts lose
  `y(4:LRS)` fails `collision_pairs` on every rank count (periodic pairs separate 10% too
  fast) and `same_as` on 2 and 4 ranks. ppiclF splits bins in z on 2 ranks and in y and z
  on 4 ranks; the x direction only uses periodic ghosts.

## Flow-File Conventions

What the five arrays of `SOL_*/<case>.mixt.cva_*` hold depends on whether particles are
on (`RFLU_TimeStepping` converts real cells just before each dump):

- PICL off: conserved bulk state `rho, rho*u, rho*v, rho*w, rho*E`
- PICL on: gas-phase primitive state `rho_g, u, v, w, E`, with `rho_g = rho/(1-phi)`

The fluid dumps (`restart_data/fluid_<step>_<region>.dat`) hold the state as it is in memory:
`(1-phi)*rho, ...` with particles; rfluinit's step 0 dump the pure gas state. Layouts of both
dump files: `modflu/RFLU_ModRestartData.F90`; readers: `rocpicl_tests/outputs.py`.

Goldens store the arrays in the old-style form above (made from the dumps). Conservation
totals convert to bulk quantities with `phi` from the fluid dump (at step 0, the particle
volume projected at the start of the run over the cell volume).

## Adding a Case

1. Add it to `rocpicl_tests/cases.py`.
2. Run `./run_tests.py run --generate -o <name>`. The golden is written only if the
   physics checks pass, or if the case has a `known_issue`.
3. Commit the new `goldens/<name>/` directory.
