# Boundary Conditions (`<case>.bc`)

Each grid patch gets its boundary condition from one section of `<case>.bc`, in the case
directory next to `<case>.inp`. [`./rocflu.sh validate`](InputFile.md#the-case-validator)
checks the file against the rules on this page.

## Contents

- [The File](#the-file)
- [Sections and Their Keys](#sections-and-their-keys)
- [Removed Sections and Keys](#removed-sections-and-keys)

## The File

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
PARTWALL   1            ! new: the patch is a wall for particle-wall collisions
#

# END
```

How Rocflu reads the file (`modflu/RFLU_ModReadBcInputFile.F90`, `ReadPatchSection` in
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

## Sections and Their Keys

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
| `BC_SYMMETRY` | `PATCH`, `NAME` | `PATCH`, `NAME` |
| `BC_VIRTUAL` | `PATCH`, `NAME`. The two z patches of a `DIMENS 2` mesh: see [Running in 2D](Running.md#running-in-2d) | `PATCH`, `NAME` |

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

## Periodic Patches

A `BC_PERIODIC` pair is a translation when `ANGLE` is 0, else a rotation by `ANGLE`
degrees about the coordinate axis `AXIS` (1 x, 2 y, 3 z), through the origin. Angles
follow the right-hand rule, and `ANGLE` is the rotation that takes the related patch
(`RELPATCH`) onto the patch of the section. For a quarter cylinder about z with patch 3
at y = 0 and patch 4 at x = 0:

```
# BC_PERIODIC            # BC_PERIODIC
PATCH     3 3            PATCH     4 4
RELPATCH  4 4            RELPATCH  3 3
ANGLE     -90.0          ANGLE     90.0
AXIS      3              AXIS      3
```

With the signs the other way round, `rflupart` stops in `OCT_SEARCH` ("lies outside
convex hull") while matching the vertices of the two patches.

With particles on, the particles take their periodicity from these patches:

- a translation along x, y or z makes the particles periodic in that direction. The
  translation must be the length of the domain in that direction;
- a rotation is printed in the Setup box (e.g. `90 deg about z`). A particle that leaves
  through one plane re-enters through the other, rotated with every vector it carries
  (the `vec` slots of `ppiclf_slots.txt`). The gas interpolation and the projection do
  not reach across the planes yet;
- the two rotational patches must be flat, contain the axis, and bound a sector of at
  most 180 degrees; all rotational pairs must have the same axis and angle, and no
  translation may lie in the plane of the rotation. The run stops at start-up otherwise.

## Removed Sections and Keys

| In `<case>.bc` | Instead |
|---|---|
| `# BC_INFLOW` | `# BC_INFLOW_TOTANG` (the same section) |
| `# BC_NOSLIP` | `# BC_NOSLIP_HFLUX` (the same section) |
| `MVPATCH`, `SMGRID` | nothing: no longer used |
| `COUPLED` | nothing: only for the GENX coupling, which is gone |
| `MOVEDIR` | nothing: only for moving grids |
| `STATS`, `BFLAG`, `FILE` | nothing: nothing reads them |
