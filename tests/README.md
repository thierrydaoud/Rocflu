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
- a particle list, written to `points.dat`;
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
| `same_as` | Particles match the same case on a different number of ranks |

## Known issues tracked by the suite

- `box/momentum-exchange/np2`: particles differ from the 1-rank run by up to ~1% in the
  interpolated pressure gradient and ~0.04% in transverse velocities (the shock tube
  agrees to 2e-6). Conservation itself passes.
- `box/quiescent-cloud`: gas starts moving (~5e-5 m/s) around fixed particles, and bulk
  mass drifts ~3e-8.

Fixed and now guarded: non-conservative periodic patches (`box/gas-jump/periodic`,
commit 11dc640).

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
