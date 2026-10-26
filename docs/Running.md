# Running a simulation

`./rocflu.sh run` checks a case, builds the code for it, pre-processes it and submits the
job; `./rocflu.sh post` converts the results for ParaView. How the code is built is in
[Getting started](GettingStarted.md), the input file and `points.dat` in [Input file](InputFile.md).

## Contents

- [`rocflu.sh run`](#rocflush-run)
- [Validate a case: `rocflu.sh validate`](#validate-a-case-rocflush-validate)
- [Running the steps by hand](#running-the-steps-by-hand)
- [Post-process: `rocflu.sh post`](#post-process-rocflush-post)

## `rocflu.sh run`

```bash
./rocflu.sh run /path/to/case/<case>.inp -N 1 -n 4 -p hpg-default -w 02:00:00 -@ you@ufl.edu
```

The case directory must hold `<case>.inp`, `<case>.bc`, the grid (`<case>.hyb.asc`) and,
with particles, `points.dat`. `run` then:

0. **Checks the case** ([`rocflu.sh validate`](#validate-a-case-rocflush-validate)). An
   error stops here, before anything is built.
1. **Builds the code for this case.** It reads the PICL and FLOWMODEL sections of the
   input and builds only the particle slots the models that are on need (see
   [Particle slots](GettingStarted.md#particle-slots)). The build is named after those
   slot groups, e.g. `build/install/release-case-qsfluct/bin`; every case with the same
   physics shares it, so a second case usually compiles nothing.
2. **Pre-processes the case here**, on the login node: `rflumap` (into N x n regions),
   `rflupart`, `rfluinit`, one process each, with their output in `rflumap.log`,
   `rflupart.log` and `rfluinit.log` in the case directory. A mistake in the input, the
   boundary conditions or the grid stops here, before anything is queued. On a terminal
   a status line shows the step each tool is in and, region by region, how many are done
   and about how long is left; Ctrl-C stops the tool.

   Steps already done are skipped. After each step that succeeds, the keys of what it
   depends on go into `<dir>/.rocflu_pre.stamp` (`utilities/scripts/preprocess_key.py`):
   rflumap and rflupart depend on the number of regions, the grid (checksum), the `.bc`
   file and the input file; rfluinit on the partition and the input file. The input file
   counts without comments and without the settings that cannot change these outputs:
   `# POST`, `# PROBE`, `CFL` and the run-length keys of `# TIMESTEP` (`MAXTIME`,
   `RUNTIME`, `WRITIME`, `PRNTIME`, `FOMWRITIME`, `TIMESTEP`, `DTMINLIMIT`, `WALLTIME`,
   `SAFEWRTIME`, `START_STEP`) for both, and also `# INITFLOW`, `# PICL`, `# REFERENCE`,
   `# MIXTURE` and `STARTTIME` for the partition. So a longer run skips everything, a new
   initial condition reruns only `rfluinit`, and new numerics rerun all three. A step is
   also redone when one of its files is missing. The executables do not count: after a
   change to the pre-processing code, use `--repartition`. A run stopped halfway leaves no
   stamp for what it was redoing.
3. **Writes `<case>.sbatch`** in the case directory and **submits it**: the job runs
   `rflump`, with `srun` on N nodes x n ranks per node, then (unless `--no-post`)
   `rflupost` on one task if the simulation succeeded. The solver's output goes to
   `<case>.<job id>.out`: a header box from the job script (start time, partition,
   account, QOS, walltime, nodes and ranks, job ID, build and commit, modules), the Setup box,
   one progress line per `PRNTIME` (percent done, t, dt, step, time per step, ETA) and a
   footer box (total time, peak memory against the Slurm request, end time). `rflump -v 2`
   adds the detailed messages of every stage.

| Option | Meaning | Default |
|---|---|---|
| `-N, --nodes N` | nodes | 1 |
| `-n, --tasks-per-node N` | MPI ranks **per node**; the case is split into N x n regions | 1 |
| `-p, --partition P` | Slurm partition | `hpg-default` |
| `-a, --account A` | Slurm account | `bala1s` |
| `-q, --qos Q` | `primary`: the account's QOS (`bala1s`); `burst`: its burst QOS (`bala1s-b`) | `primary` |
| `-w, --walltime T` | wall time `HH:MM:SS` | `01:00:00` |
| `--mem M` | memory per node, e.g. `64gb` | 2 GB per rank (`--mem-per-cpu=2gb`) |
| `-@, --email A` | send all Slurm notifications (start, end, failure) to A | none |
| `-t pre_process` | build and pre-process only, no job | |
| `-t simulation` | build and submit only; the case was pre-processed before with the same N x n | |
| `-e`, `--engine` | `batch`: submit the job script with `sbatch`; `interactive`: run it right away in this shell, on a node you have already allocated (`srun --pty bash` or `salloc`; refused elsewhere), output also in `<case>.interactive.out` | `batch` |
| `--repartition` | pre-process again even when the stamp says it is up to date | skip what is done |
| `--dry-run` | build, pre-process and write the job script, but do not submit or run it | |
| `--no-post` | do not run `rflupost` after the simulation. By default, after a simulation that succeeded, the job runs `rflupost` (one task), so `post/` is ready when the job ends | post on |
| `--debug` | debug build of the case | |
| `-j, --jobs N` | compile with N processes | 16 |

`run` does not change `build/bin`, which stays the full build `rocflu.sh test` uses.
Pre-processing runs on the login node, which is fine for test-size grids; very large
grids (a million cells and more) can take long there (to be improved).

## Validate a case: `rocflu.sh validate`

```bash
./rocflu.sh validate /path/to/case/<case>.inp
```

`utilities/scripts/case_validator.py` reads the input file the way Rocflu does and
reports every problem at once, with the line it is on:

- keys in `# PICL` that Rocflu would ignore or misread: unknown keys (with the closest
  known key), removed keys (with what replaces them), a key after a space, a key given
  twice, and `QSFLUCTFILT` before `QSFLUCT` (Rocflu reads it as `QSFLUCT`);
- flags that are not whole numbers or out of range;
- particle models that do not work together, e.g. `QSFLUCT` without `QSFORCE`;
- with `USED 1`: `RKSCHEME 2` in `# TIMESTEP`, `VISCOSITY` in `# VISCMODEL`, no
  `AXIFLAG 1`, `DENSITY` and `SPECIFICHEAT` in `# PICL`, and no `filein.vtk` or
  `par00002.vtu` (which would start a particle restart) in the case directory;
- `points.dat`: `!` comment lines, the number of particles, then exactly that many
  lines of `x y z vx vy vz dp` (position, initial velocity, diameter) with `dp > 0`
  (lines in an older form, with a material name first or with `x y z dp` only, are
  reported as such), and every particle inside the fluid domain. The domain
  is the box around the grid's vertices, read from `<case>.hyb.asc` (CENTAUR ASCII,
  `GRIDSRC 0`; about a second for a million vertices). In a domain that is not a box
  (a cylinder, a wedge), a particle can be inside the box but outside the grid.

Settings of a model that is off (e.g. `KSP` without collisions) are ignored. `rocflu.sh
run` and the regression suite run the same checks. The validator's own tests:
`python3 tests/test_case_validator.py`.

## Running the steps by hand

With the modules loaded (`source ./rocflu.sh load`) and N MPI ranks, using a full build
(`build/bin`) or a case build (`build/install/<name>/bin`):

```bash
B=/path/to/repo/build/bin
$B/rflumap  -c <case> -m 1 -p N -r N -v 1
$B/rflupart -c <case> -v 1
$B/rfluinit -c <case> -v 1
srun -n N $B/rflump -c <case> -v 1      # or mpirun -n N
```

`rflump -v` sets how much the run log shows: `1` (default) only the run's progress,
warnings and errors; `2` the detailed messages of every stage (Rocflu and ppiclF);
`3` debugging detail as well. The pre-processing tools keep their own levels.

A case build stops at start-up if the input asks for a model whose slots it does not
have, with a message such as *Input needs ppiclF slots 'qsfluct' (QSFLUCT >= 1), not in
this build: use ./rocflu.sh run*.

The post program runs on one process, in the case directory:
`$B/rflupost -c <case> [-s new|all|last|N,M]`.

## Post-process: `rocflu.sh post`

```bash
./rocflu.sh post /path/to/case/<case>.inp               # the steps not yet in post/
./rocflu.sh post /path/to/case/<case>.inp -s all        # every step again
./rocflu.sh post /path/to/case/<case>.inp -e batch      # as a one-task Slurm job
```

`rflump` writes its results as raw dumps in `restart_data/` (the same files a restart
reads, see [Input file](InputFile.md#output-and-restarts)). `post` converts them into
ParaView files in `post/`, with the case's build: open `post/collection.pvd`, which lists
every step with its time. Each step is one dataset, `post/root/collection_<step>.vtm`,
with three blocks: the fluid of every region (`post/p<k>/<step>.vtu`, region k+1), the
particles (`post/particles/<step>.vtu`) and the ppiclF bins with the rank that owns each
(`post/bins/<step>.vtu`). All values are in double precision.

To show, hide and filter the fluid, the particles and the bins each on its own, open
`post/fluid.pvd`, `post/particles.pvd` and `post/bins.pvd` together instead (one series
each, same steps and times; ParaView lists them as separate objects in the Pipeline
Browser). In `collection.pvd` the three kinds are blocks of one object: View >
MultiBlock Inspector shows or hides them.

By default only the steps without a dataset in `post/` are converted, so `post` can run
while the simulation is still going and again later; `-s all` redoes every step, e.g.
after changing which fields to write (the `*_WRT` keys of `# POST` for the fluid,
`PARTICLE_WRT`, ... of `# PICL` for the particles).

| Option | Meaning | Default |
|---|---|---|
| `-s, --steps S` | `new`, `all`, `last`, or step numbers `N,M,...` | `new` |
| `-e, --engine E` | `here`: run on this node; `batch`: submit `<case>.post.sbatch` (one task) | `here` |
| `-p`, `-a`, `-q`, `-w`, `--mem` | partition, account, QOS (`primary` or `burst`), wall time, memory of the batch job | `hpg-default`, `bala1s`, `primary`, `00:30:00`, `8gb` |
| `--debug`, `-j N` | build options, as for `run` | |
