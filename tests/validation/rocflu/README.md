# Running Rung 1 With Rocflu

The only Rocflu-specific part of the validation tests. The problems, exact solutions and checks in `../problems`, `../exact` and `../metrics` do not depend on Rocflu.

## Files

| File | What it does |
|---|---|
| `make_cases.py` | Writes case directories from `tests/templates/rectshktb.*`. |
| `run_cases.sh` | Meshes, validates, pre-processes, runs and post-processes each listed case, then writes `result.txt`. |
| `to_columns.py` | Converts a run's last dump into the plain columns `metrics/compare.py` reads. |
| `units.py` | Scales between the dimensionless problems and Rocflu's SI input. |

## Steps

From this directory, with a built Rocflu (`./rocflu.sh build`; executables in `build/bin`):

```bash
python3 make_cases.py /blue/$USER/rung1                       # 96 cases: 2 problems x 4 schemes x 2 orders x 2D/3D x 3 grids
python3 make_cases.py /blue/$USER/rung1 --discr 3 --n 200     # or a subset
sbatch -n 8 -t 02:00:00 --wrap "bash $PWD/run_cases.sh $PWD/../../../build/bin /blue/$USER/rung1/cases.txt"
```

Then, for each scheme, order and dimension:

```bash
R=/blue/$USER/rung1
python3 ../metrics/compare.py sod --order 1 $R/sod_d3_o1_dim2_n{200,400,800}/result.txt
python3 ../metrics/compare.py sod --same dims $R/sod_d3_o1_dim2_n800/result.txt $R/sod_d3_o1_dim3_n800/result.txt
```

For rank independence, run the same case on two rank counts (two case roots, two `sbatch -n` values) and use `--same ranks`.

## The Rocflu setup

| Item | Setting |
|---|---|
| Grid | n × 8 cubic cells in x–y; 8 cells in z in 3D, 1 in 2D (`DIMENS 2`) |
| Initial state | `INITFLOW`: `RVAL1` jump position, `RVAL2–4` left ρ, u, p, `RVAL5–6` right ρ, p (right u is 0) |
| Left boundary | `BC_INFLOW_VELTEMP TYPE 1`: velocity and temperature of the left state, pressure from inside |
| Right boundary | `BC_OUTFLOW TYPE 0` (zero gradient), `KIND 0` |
| y, z | periodic; z `BC_VIRTUAL` in 2D |
| Numerics | `CFL 0.5`, RK3, `RECONST 2`, `CRECONSTC 0`, `CRECONSTF 0`, gas only (`PICL USED 0`) |
| Units | p in units of 10⁵ Pa, ρ of 1 kg/m³, 1 m domain: velocity scale 316.228 m/s |

Notes:

- **`KIND 1` (NSCBC) on the outflow** is acted on only by Roe (`DISCR 1`), where it grows from the outflow face and crashes the fine grids. HLLC and AUSM ignore `KIND`. The cases use `KIND 0`.
- **Rocflu's CFL number** sums the wave speeds over all faces of a cell, so `CFL 0.5` is about 0.13 in 2D and 0.09 in 3D in the 1D sense. On Sod at n = 800 Rocflu takes 2613 steps in 2D and 3823 in 3D, against 699 for `godunov.py` at CFL 0.5; `compare.py` runs Godunov at Rocflu's step count.
- **A run can stop on an error with exit status 0.** `to_columns.py` refuses a run whose `rflump.log` does not say "Finished".
