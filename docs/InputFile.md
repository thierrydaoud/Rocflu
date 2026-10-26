# The case input file (`<case>.inp`), `<case>.bc` and `points.dat`

A Rocflu case is described by `<case>.inp` (this page), `<case>.bc` (boundary
conditions), the grid (`<case>.hyb.asc`) and, with particles, `points.dat`. Before a run,
[`./rocflu.sh validate`](#the-case-validator) checks the input file, the `.bc` file and `points.dat`;
`./rocflu.sh run` does this first and stops before building anything if it finds a
problem.

The template used by the regression tests, [`tests/templates/rectshktb.inp`](../tests/templates/rectshktb.inp),
is a complete, working example.

## Contents

- [How Rocflu reads the file](#how-rocflu-reads-the-file)
- [Gas-side sections](#gas-side-sections)
- [`# PICL`: the particles](#picl-the-particles)
  - [Models that need other models](#models-that-need-other-models)
- [`<case>.bc`: the boundary conditions](#casebc-the-boundary-conditions)
- [`points.dat`: the particles' initial state](#pointsdat-the-particles-initial-state)
- [Output and restarts](#output-and-restarts)
- [The case validator](#the-case-validator)
  - [Removed keys](#removed-keys)

## How Rocflu reads the file

```
# PICL
USED             1 ! comment
QSFORCE          3
#
```

- A **section** starts at a line that is exactly `# NAME` (one space, no other text; `#PICL`
  or `# PICL extra` is not recognised) and ends at the next line that starts with `#`.
  Sections Rocflu does not know are skipped.
- A **key** must start in the first column: a line with a space before the key is
  ignored. The value follows the key; anything after it (usually a `!` comment) is
  ignored. Lines starting with `!` are comments.
- A key that is **not given** keeps its default (below). A key given twice: the first
  one counts.
- Rocflu matches a line to the **first key not yet read that the line starts with**. In
  `# PICL` this matters for `QSFLUCT`, `QSFLUCTFILT` and `QSFLUCTFILTADAPT`: write them in
  that order, or a `QSFLUCTFILT` line is read as `QSFLUCT`.
- Flags are whole numbers; Rocflu rounds other values silently (`2.6` becomes `3`).
- Rocflu silently ignores a key it does not know, so a typo (`QSFLUC 1`) turns a model off.

The validator reports all of these (with the line number), so they need not be
remembered.

## Gas-side sections

Only the keys in the test template are listed; the readers (`libfloflu/Read*Section.F90`)
accept more.

| Section | Key | Values | Meaning |
|---|---|---|---|
| `# FORMATS` | `GRID`, `SOLUTION` | 0 ASCII, 1 binary, 2 HDF | Rocflu's own grid files; `SOLUTION`: the old-style solution files (`OLD_OUTPUT_WRT 1` only) |
| | `GRIDSRC` | 0 CENTAUR ASCII, 1 VGRIDNS, 2 MESH3D | grid read by the pre-processing (`<case>.hyb.asc` for 0) |
| `# FLOWMODEL` | `MODEL` | 0 Euler, 1 Navier-Stokes | |
| | `MOVEGRID` | 0 static, 1 moving | |
| `# NUMERICS` | `CFL` | > 0 | CFL number |
| | `DISCR` | 0 central, 1 Roe, 2 MAPS, 3 HLLC, 4 AUSM+, 5 AUSM+up | flux scheme |
| | `ORDER` | 1, 2 | order of accuracy in space |
| | `ENTROPY` | | entropy correction (Roe) |
| | `DIMENS` | 1, 2, 3 | dimension of the computation |
| | `AXIFLAG` | 0, 1 | axisymmetric; **not available with particles** |
| | `SIGMA`, `KP`, `KU`, `MINF` | | AUSM+up parameters (defaults 1.0, 0.25, 0.75, 0.1) |
| | `RECONST` | 0 none, 1 WENO simple, 2 WENO XYZ, 10 Barth-Jespersen, 11 Venkatakrishnan | gradient reconstruction (limiter) |
| | `CRECONSTC`, `CRECONSTF` | 0, 1 | constrained cell / face gradients |
| | `CRECONSTCW`, `CRECONSTFW` | 0, 1, or a weight | their weights |
| `# TIMESTEP` | `FLOWTYPE` | 0 steady, 1 unsteady | |
| | `SOLVERTYPE` | 1 dissipative, 3 non-dissipative | |
| | `MAXTIME`, `RUNTIME` | s | end time; run time of this job |
| | `WRITIME`, `PRNTIME`, `FOMWRITIME` | s | intervals for solution files, printing, FOM files |
| | `TIMESTEP` | s | largest time step |
| | `DTMINLIMIT` | s | smallest time step, overrides the CFL limit |
| | `STARTTIME` | s | |
| | `WALLTIME`, `SAFEWRTIME` | s | wall time; write a restart this long before it |
| | `RKSCHEME` | 1 RK4 (default), 2 RK3 | **must be 2 with particles** |
| | `START_STEP` | 0 (default), -1, N | 0: a new run from rfluinit's step 0, replacing the results of an earlier run in the directory; -1: restart from the latest dump in `restart_data/`; N: restart from the dump of step N. Dumps newer than the start step are removed. A restart needs the same number of ranks (regions) and the same particle physics as the run that wrote the dump, and continues it exactly: it gives the same results, to the last bit, as a run that never stopped |
| `# INITFLOW` | `FLAG` | 1 from scratch, 2 from file, 3 hard-coded, 4-6 combinations | initial gas state |
| | `RVAL1`…, `IVAL1`… | | parameters of the initial state; with `FLAG 3` their meaning depends on the case name (`utilities/init/RFLU_InitFlowHardCode.F90`; for `rectshktb`: shock position, post- and pre-shock state) |
| `# REFERENCE` | `GAMMA`, `CP` | | ratio of specific heats; specific heat at constant pressure [J/(kg K)] |
| `# MIXTURE` | `GASMODEL` | 1 calorically perfect gas; 3-7 mixtures (7 JWL) | 3-7 need species transport, **disabled for now** (rejected) |
| `# SPECIES` | `USED` | 0 | species transport, **disabled for now**: `USED 1` is rejected. Species materials are given in `# MATERIAL` sections (not used by the particles) |
| `# VISCMODEL` | `MODEL` | 0 Sutherland, 1 fixed | |
| | `VISCOSITY` | Pa s | reference viscosity; **required with particles** |
| | `REFTEMP`, `SUTHCOEF` | K | Sutherland's law |
| `# POST` | `DENSITY_WRT`, `VELOCITY_WRT`, `PRESSURE_WRT`, `PHI_WRT` | 0, 1 (default 1) | fluid fields the post program writes: gas density, velocity, pressure, particle volume fraction (with particles) |
| | `TEMPERATURE_WRT`, `SOUNDSPEED_WRT`, `ENERGY_WRT`, `GRADRHO_WRT`, `PT_WRT` | 0, 1 (default 0) | also temperature, speed of sound, total energy per mass, gas-density gradient, pseudo-turbulence fields `DivPhiRsg` and `FCell` (with `PSEUDOTURB`) |
| | `OLD_OUTPUT_WRT` | 0 (default), 1 | also the old-style output (`SOL_*`, `PARAVIEW_*`, `par*.vtu`, `bin*.vtu`, patch coefficients); needed with `MOVEGRID 1` or `SOLVERTYPE 3`, which write into `SOL_*` |
| | other keys | | settings of the old-style ParaView output, as in the template |
| `# MISC` | | | miscellaneous settings, as in the template |

## `# PICL`: the particles

Defaults are those of `libflu/RFLU_InitGlobal.F90`. A key that belongs to a model that is
off (e.g. `KSP` without collisions) is not read and may stay in the file with any value.

| Key | Values | Default | Meaning |
|---|---|---|---|
| `USED` | 0, 1 | 0 | particles on; needs a build with ppiclF (not `--no-picl`) |
| `VELOCITYMODEL` | 0, 1 | 1 | 1: particles move under the forces. 0: they keep their position, velocity, temperature and spin; the forces are still computed, so with `FEEDBACK 1` the gas still feels them (a fixed bed) |
| `QSFORCE` | 0-3 | 0 | quasi-steady drag: 0 none, 1 Gidaspow, 2 Parmar, 3 Osnes |
| `QSFLUCT` | 0-2 | 0 | fluctuation of the quasi-steady force: 0 none, 1 Lattanzi, 2 Osnes |
| `QSFLUCTFILT` | 0, 1 | 0 | filter for the fluctuation: 0 box, 1 Gaussian |
| `QSFLUCTFILTADAPT` | 0, 1 | 0 | adaptive filter width |
| `AMFORCE` | 0-2 | 0 | added mass: 0 none, 1 Parmar, 2 Briney (binary; **not validated, rejected**) |
| `PGFORCE` | 0, 1 | 0 | pressure-gradient (undisturbed-flow) force |
| `VISCOUSUNSTEADY` | 0-4 | 0 | viscous-unsteady (history) force: 1 trapezoidal for uniform dt, 2 trapezoidal for variable dt, 3 modified trapezoidal, 4 Hinsberg |
| `VUSTEPS` | ≥ 1 | 50 | time steps of history kept for `VISCOUSUNSTEADY` |
| `COLLISIONS` | 0-4 | 0 | 1 normal contact force; 2 + tangential force and collision torque (particles spin); 3 + rolling torque and hydrodynamic torques; 4 + Saffman and Magnus lift |
| `KSP` | > 0 | (none) | collisions: a contact lasts `KSP` time steps (sets the spring stiffness); needed with `COLLISIONS ≥ 1` |
| `EREST` | 0 < e ≤ 1 | (none) | collisions: coefficient of restitution; needed with `COLLISIONS ≥ 1` |
| `HEATTRANSFER` | 0-5 | 0 | quasi-steady heat transfer (Nusselt number): 1 Stokes, 2 Ranz-Marshall, 3 Gunn, 4 Fox, 5 Sun |
| `DU HEATTRANSFER` | 0, 1 | 0 | diffusive-unsteady heat transfer (computed with the viscous-unsteady history) |
| `FEEDBACK` | 0, 1 | 0 | two-way coupling: the particles' forces and heat act on the gas |
| `PSEUDOTURB` | 0, 1 | 0 | pseudo-turbulence model: Reynolds stress, kinetic energy and heat flux the particles induce in the gas |
| `TEMP` | > 0 | 0 | initial temperature of every particle [K]; **required** |
| `DENSITY` | > 0 | 0 | density of every particle [kg/m³]; **required** |
| `SPECIFICHEAT` | > 0 | 0 | specific heat of every particle [J/(kg K)]; **required** |
| `NEIGHBORWIDTH` | ≥ 0 | 0 | neighbour search distance; the code uses max(4 × largest diameter, `NEIGHBORWIDTH`) |
| `PERIODICX`, `PERIODICY`, `PERIODICZ` | 0, 1 | 0 | particles periodic in x, y, z; should match the periodic patches of the fluid in `<case>.bc` (not checked yet) |
| `SUBBIN` | 0, 1 | 0 | neighbour-search method: 0 Zwick, 1 Avery sub-bins |
| `PARTICLE_WRT` | 0, 1 | 0 | particle files also hold temperature, spin, density, volume |
| `GAS_WRT` | 0, 1 | 0 | … and the gas density, velocity, sound speed and temperature at each particle |
| `GRADGAS_WRT` | 0, 1 | 0 | … and the gas gradients at each particle |
| `MODEL_WRT` | 0, 1 | 0 | … and the states and diagnostics of the models that are on |

The filter width of the particle projection is not an input: it is twice the largest cell
length around the particles. `./rocflu.sh run` builds only the particle data the models
that are on need (`libpicl/user_files/ppiclf_slots.txt`, see
[Getting started](GettingStarted.md#particle-slots)).

### Models that need other models

With `USED 1`, these combinations are rejected:

| Setting | Needs | Why |
|---|---|---|
| `QSFLUCT ≥ 1` | `QSFORCE ≥ 1` | the fluctuation is added to the quasi-steady force |
| `PSEUDOTURB 1` | `QSFORCE ≥ 1` | the model uses the drag |
| `PSEUDOTURB 1` | `QSFLUCT 2` | |
| `PSEUDOTURB 1` | `HEATTRANSFER ≥ 1` | for now: its heat flux divides by the Nusselt number. To be lifted when the pseudo-turbulence model gets its own switch for the heat flux |
| `DU HEATTRANSFER 1` | `VISCOUSUNSTEADY ≥ 1` | they share the history |
| `COLLISIONS ≥ 1` | `VELOCITYMODEL 1` | particles that keep their position cannot collide |
| `AMFORCE 2` | | not validated yet |

Other sections, with `USED 1`: `RKSCHEME 2` in `# TIMESTEP`, `VISCOSITY` in
`# VISCMODEL`, no `AXIFLAG 1` in `# NUMERICS`.

## `<case>.bc`: the boundary conditions

Each grid patch gets its boundary condition from one section of `<case>.bc`:

```
# BC_OUTFLOW
PATCH      1 3          ! first and last patch: always two numbers (PATCH 3 3 for one patch)
NAME       Outflow
TYPE       0
KIND       1
ORDER      1
#                       ! every section ends with a line '#'

# BC_SLIPW
PATCH      4 7
NAME       Walls
PARTWALL   1
#

# END
```

How Rocflu reads it (`modflu/RFLU_ModReadBcInputFile.F90`, `ReadPatchSection` in
`libfloflu/ReadInputUtil.F90`), and why each rule below exists:

- A header is the whole line `# NAME`, with one space. An unknown header is skipped
  silently, and its patches get no boundary condition.
- A section ends at the next line that starts with `#`, and Rocflu reads that line as
  part of the section. If the next header follows directly, it is used up as the end,
  and that whole section is skipped.
- `PATCH` is read as two whole numbers. With one number (`PATCH 3`, or `PATCH 3 ! comment`)
  Rocflu stops with `Fortran runtime error: End of file`. A comment after the two numbers
  is fine.
- A key is taken from a line that *starts with* it. A key the section does not know, or
  an indented one, is ignored; `PRESSURE` would be read as `PRESS`.
- Values out of range are replaced by defaults without a word (`KIND 5` becomes 0,
  `ORDER 3` becomes 1).

### Sections and their keys

Every section takes `PATCH` and `NAME` (both required) and `DISTRIB` (0: the values below;
1: a value per face from the file that `rfluinit` writes, so the data values are not needed).

| Section | Keys (whole-number values) | Required |
|---|---|---|
| `BC_OUTFLOW` | `TYPE` (0 supersonic, 1 subsonic, 2 mixed), `PRESS`, `KIND` (0, 1), `REFLECT` (0, 1), `NSCBCK`, `ORDER` (1, 2), `THRUSTFLAG` (0, 1) | `TYPE`; `PRESS` unless `TYPE 0` |
| `BC_FARF` | `CORR` (0, 1), `MACH`, `ATTACK`, `SLIP`, `PRESS`, `TEMP`, `KIND`, `ORDER`, `THRUSTFLAG` | `CORR`, `MACH`, `ATTACK`, `SLIP`, `PRESS`, `TEMP` |
| `BC_INFLOW_TOTANG` | `TYPE` (0, 1, 2), `FIXED` (0, 1), `PTOT`, `TTOT`, `BETAH`, `BETAV`, `MACH`, `KIND`, `REFLECT`, `ORDER`, `THRUSTFLAG` | `TYPE`, `PTOT`, `TTOT`, `BETAH`, `BETAV`; `MACH` unless `TYPE 1` |
| `BC_INFLOW_VELTEMP` | `TYPE` (0, 1, 2), `VELX`, `VELY`, `VELZ`, `TEMP`, `PRESS`, `KIND`, `REFLECT`, `ORDER`, `THRUSTFLAG` | `TYPE`, `VELX`, `VELY`, `VELZ`, `TEMP`; `PRESS` unless `TYPE 1` |
| `BC_INJECT` | `MFRATE`, `TEMP`, `CRECONST` (0, 1), `KIND`, `ORDER`, `THRUSTFLAG` | `MFRATE`, `TEMP` |
| `BC_NOSLIP_HFLUX` | `HFLUX`, `CRECONST`, `KIND`, `ORDER`, `THRUSTFLAG`, `PARTWALL` (0, 1) | `HFLUX` |
| `BC_NOSLIP_TEMP` | `TEMP`, `CRECONST`, `KIND`, `ORDER`, `THRUSTFLAG`, `PARTWALL` | `TEMP` |
| `BC_SLIPW` | `CRECONST`, `KIND`, `ORDER`, `THRUSTFLAG`, `PARTWALL` | |
| `BC_PERIODIC` | `RELPATCH`, `AXIS` (1, 2, 3), `ANGLE` | all three |
| `BC_SYMMETRY`, `BC_VIRTUAL` | | |

`KIND 1` is NSCBC (non-reflecting characteristic BC; 0 is the simple one). `ORDER` is the
order of the boundary reconstruction; without it the patch uses `ORDERBF` from
`# NUMERICS` (default 1). `PARTWALL 1` makes the patch a wall for the particles.
`THRUSTFLAG` only matters with `# FORCES` on.

Rules across sections:

- every patch of the grid in exactly one section. Rocflu drops the grid patches that have
  no faces and renumbers the others, so the patch numbers are the ones after that;
- `BC_PERIODIC` takes one patch per section, its `RELPATCH` is a `BC_PERIODIC` patch
  that points back at it, the pair has the same `AXIS`, and the two `ANGLE`s add up to 0;
- boundary gradients exist only at first order (`RFLU_CreateWtsBF2C`). An NSCBC patch
  (`KIND 1`) needs them, and with `MODEL 1` (Navier-Stokes) every patch needs them, so
  those patches must have order 1 (`ORDER 1`, or `ORDERBF 1` when `ORDER` is not given).
  Otherwise `rflump` stops in `RFLU_CreateWtsBF2C` with "reached default statement".

### Removed sections and keys

| In `<case>.bc` | Instead |
|---|---|
| `# BC_INFLOW` | `# BC_INFLOW_TOTANG` (the same section) |
| `# BC_NOSLIP` | `# BC_NOSLIP_HFLUX` (the same section) |
| `MVPATCH`, `SMGRID` | nothing: no longer used |
| `COUPLED` | nothing: only for the GENX coupling, which is gone |
| `MOVEDIR` | nothing: only for moving grids |
| `STATS`, `BFLAG`, `FILE` | nothing: nothing reads them |

## `points.dat`: the particles' initial state

`points.dat` (this name, in the case directory) lists the particles:

```
! optional comment lines, at the top only
3
1.0e-3  2.0e-3  3.0e-3    0.0  0.0  0.0    1.0e-4
4.0e-3  5.0e-3  6.0e-3   10.0  0.0  0.0    1.0e-4
...
```

The first line after the comments is the number of particles; then one line per
particle: position `x y z` [m], initial velocity `vx vy vz` [m/s] (0 0 0 for a particle at
rest) and diameter `dp` [m], all seven on every line. Density, specific heat and initial
temperature are the same for every particle (`DENSITY`, `SPECIFICHEAT`, `TEMP` in
`# PICL`).

Older forms are no longer read: a material name in the first column (removed with
`# MATERIAL` for particles), and `x y z dp [u v w]` (diameter before the velocity).
The validator and `rflump` both stop on them.

## Output and restarts

`rflump` writes its results as raw dumps in `restart_data/`, at step 0 (the particles),
every `WRITIME` and the end of the run:

- `fluid_<step>_<region>.dat`: the gas state of a region as it is in memory (real cells,
  double precision), with particles also the particle volume fraction and the gas-density
  gradient, with `PSEUDOTURB` the pseudo-turbulence fields, and the face states of the
  non-reflecting (NSCBC) boundaries. `rfluinit` writes step 0.
- `particles_<step>.dat`: every particle (state, properties, IDs, random-stream state,
  viscous-unsteady history), the fields projected to the fluid cells, the bins and the
  list of output fields.

Layouts: `modflu/RFLU_ModRestartData.F90`. A file is written as `.tmp` and renamed when it
is complete. The same files are what a restart reads and what the post program converts
(`./rocflu.sh post`, see [Running a simulation](Running.md#results-rocflush-post)), so results
can be looked at during a run or after it stopped.

**Starting a run** (`START_STEP` in `# TIMESTEP`): `0` (default) is a new run from
rfluinit's step 0; it removes what an earlier run in the directory left (dumps after step
0, the old-style output). `-1` restarts from the latest step with every dump, `N` from
step N; dumps after the start step are removed. A restart needs the same number of ranks
and the same particle physics as the run that wrote the dump, and continues it exactly:
the results are the same, to the last bit, as a run that never stopped (checked by the
suite's restart test).

The old-style output (`SOL_*`, `PARAVIEW_*`, `par*.vtu`, `bin*.vtu`, patch coefficients)
is written only with `OLD_OUTPUT_WRT 1` in `# POST`; the suite keeps it on in a few cases
to check that it agrees with the dumps and the post program. Without particles its
`Velocity` and `Energy` fields hold the momentum and the energy per volume; the post
program writes velocity and energy per mass.

## The case validator

```bash
./rocflu.sh validate /path/to/case/<case>.inp
python3 utilities/scripts/case_validator.py /path/to/case/<case>.inp      # the same
```

[`utilities/scripts/case_validator.py`](../utilities/scripts/case_validator.py) reads
the input file the way Rocflu does and prints every problem at once, with the file and
line, then exits with status 1 (0 when the case is fine). `./rocflu.sh run` runs it
before building, and the regression suite runs it on every case. It checks:

- `# PICL`: unknown keys (with the closest known key as a hint), [removed keys](#removed-keys),
  keys after a space, keys given twice, keys Rocflu would read as another key;
- flags that are not whole numbers or outside the values in the table above;
  `TEMP`, `DENSITY`, `SPECIFICHEAT` > 0; `NEIGHBORWIDTH` ≥ 0; with collisions `KSP` > 0
  and 0 < `EREST` ≤ 1; with viscous-unsteady `VUSTEPS` ≥ 1;
- the [combinations](#models-that-need-other-models) above, and the gas-side settings
  particles need;
- `# SPECIES USED 1` and `GASMODEL` 3-7: species transport is disabled for now;
- the case directory: `points.dat` is present and in the [format](#pointsdat-the-particles-initial-state)
  above (the count matches the lines, `dp` > 0); no `filein.vtk` (particle walls come from
  `PARTWALL 1` in the `.bc` file); no `par00002.vtu` from an earlier run (it would start a
  particle restart, which is not available yet);
- the `.bc` file: every rule of [`<case>.bc`](#casebc-the-boundary-conditions) above
  (sections, keys, values, required keys, every grid patch once, periodic pairs, order 1
  where boundary gradients are needed). The patch count comes from the grid
  (`<case>.hyb.asc`, `GRIDSRC 0`); without it, patches 1 to the highest one used must
  all be there;
- every particle inside the fluid domain, taken as the box around the grid's vertices
  (`<case>.hyb.asc`, `GRIDSRC 0`; about a second for a million vertices). In a domain
  that is not a box, a particle can be inside the box but outside the grid.

Without particles (`USED 0` or no `# PICL`), only the keys of `# PICL` are checked.

The validator mirrors `rocpicl/PICL_ReadPiclSection.F90` (keys and their order),
`libflu/RFLU_InitGlobal.F90` (defaults) and `modflu/RFLU_ModReadBcInputFile.F90` (the
`.bc` sections and keys): change it with them. Its tests need no cluster:
`python3 tests/test_case_validator.py`.

### Removed keys

| Key | Instead |
|---|---|
| `STATIONARY` | `VELOCITYMODEL` (`STATIONARY 1` is `VELOCITYMODEL 0`) |
| `FILTERWIDTH` | nothing: the filter width is set from the grid |
| `REACTIVE` | nothing: the particle burn-rate model was removed |
| `PICLDEBUG` | nothing |
| `ANGULARPERIODIC`, `ANGLE`, `XANGLE`, `RIN`, `ROUT` | nothing: angular periodicity is not available |
| `OUTPARTICLE`, `OUTGAS`, `OUTGRADGAS`, `OUTMODEL` | renamed `PARTICLE_WRT`, `GAS_WRT`, `GRADGAS_WRT`, `MODEL_WRT` |
| `PRNTENDTIME` (`# TIMESTEP`) | nothing: the end of every run is written |
| `# MATERIAL` for particles | `DENSITY` and `SPECIFICHEAT` in `# PICL` (`# MATERIAL` is still read for species) |
