"""Write Rocflu case directories for rung 1 (Sod, Lax) from tests/templates/rectshktb.*

    python3 make_cases.py <root> [--problems sod,lax] [--discr 1,3,4,5] [--order 1,2] [--dim 2,3] [--n 200,400,800]

Writes <root>/<problem>_d<DISCR>_o<ORDER>_dim<DIMENS>_n<n>/ with rectshktb.inp, rectshktb.bc and
mesh.in (answers for tests/mesh_generator), and <root>/cases.txt listing them, for run_cases.sh.

Setup (see ../problems/): slab of n x 8 cubic cells in x-y (and 8 in z in 3D); y and z periodic, z
virtual in 2D; left BC_INFLOW_VELTEMP TYPE 1 holding the left state; right BC_OUTFLOW TYPE 0
(zero gradient) with KIND 0 (KIND 1, NSCBC, is acted on only by Roe and crashes on fine grids).
CFL 0.5, RK3, RECONST 2, CRECONSTC/F 0, gas only (PICL USED 0).
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "exact"))
from rocpicl_tests.case import TEMPLATES, bc_text, edit_inp  # noqa: E402
from riemann import PROBLEMS  # noqa: E402
import units  # noqa: E402

CROSS = 8  # cells across in y (and z in 3D)


def write_case(root, name, discr, order, dim, n):
    d = PROBLEMS[name]
    left, right = units.si(d["left"]), units.si(d["right"])
    t_end = d["t"] * units.T_SCALE
    if right[1] != 0:
        raise ValueError("rectshktb sets the right velocity to zero")
    case = "%s_d%d_o%d_dim%d_n%d" % (name, discr, order, dim, n)
    rd = os.path.join(root, case)
    os.makedirs(rd, exist_ok=True)
    over = {
        ("NUMERICS", "CFL"): 0.5, ("NUMERICS", "DISCR"): discr, ("NUMERICS", "ORDER"): order,
        ("NUMERICS", "DIMENS"): dim, ("NUMERICS", "RECONST"): 2,
        ("NUMERICS", "CRECONSTC"): 0, ("NUMERICS", "CRECONSTF"): 0,
        ("TIMESTEP", "MAXTIME"): "%.10E" % t_end, ("TIMESTEP", "RUNTIME"): "%.10E" % t_end,
        ("TIMESTEP", "WRITIME"): "%.10E" % t_end, ("TIMESTEP", "PRNTIME"): "%.10E" % (t_end / 10),
        ("TIMESTEP", "FOMWRITIME"): "1.0",
        ("INITFLOW", "RVAL1"): repr(float(d["x0"])),
        ("INITFLOW", "RVAL2"): repr(left[0]), ("INITFLOW", "RVAL3"): repr(left[1]), ("INITFLOW", "RVAL4"): repr(left[2]),
        ("INITFLOW", "RVAL5"): repr(right[0]), ("INITFLOW", "RVAL6"): repr(right[2]),
        ("PICL", "USED"): 0,
    }
    if dim == 2:
        over[("PICL", "PERIODICZ")] = 0
    with open(os.path.join(TEMPLATES, "rectshktb.inp")) as f:
        inp = edit_inp(f.read(), over)
    gas_r = units.CP * (d["gamma"] - 1) / d["gamma"]
    bc = bc_text("shocktube-2d" if dim == 2 else "shocktube")
    for key, value in (("PRESS", left[2]), ("TEMP", left[2] / (left[0] * gas_r)), ("VELX", left[1])):
        bc, k = re.subn(r"(?m)^%s\s+.*$" % key, "%-14s %r" % (key, value), bc, count=1)
        assert k == 1, key
    bc, k = re.subn(r"(?m)^KIND\s+1\s*$", "KIND       0", bc)
    assert k == 1, "outflow KIND"
    with open(os.path.join(rd, "rectshktb.inp"), "w") as f:
        f.write(inp)
    with open(os.path.join(rd, "rectshktb.bc"), "w") as f:
        f.write(bc)
    dx = 1.0 / n
    nz = CROSS if dim == 3 else 1
    answers = [1.0, CROSS * dx, nz * dx, 0.0, 0.0, 0.0, n + 1, CROSS + 1, nz + 1, 2]
    with open(os.path.join(rd, "mesh.in"), "w") as f:
        f.write("\n".join(repr(v) if isinstance(v, float) else str(v) for v in answers) + "\n")
    with open(os.path.join(rd, "meta.json"), "w") as f:
        json.dump(dict(problem=name, discr=discr, order=order, dim=dim, n=n), f)
    return case


def main():
    ints = lambda s: [int(v) for v in s.split(",")]
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root")
    ap.add_argument("--problems", default="sod,lax", type=lambda s: s.split(","))
    ap.add_argument("--discr", default=[1, 3, 4, 5], type=ints, help="1 Roe, 3 HLLC, 4 AUSM+, 5 AUSM+UP")
    ap.add_argument("--order", default=[1, 2], type=ints)
    ap.add_argument("--dim", default=[2, 3], type=ints)
    ap.add_argument("--n", default=[200, 400, 800], type=ints)
    a = ap.parse_args()
    cases = [write_case(a.root, p, dc, o, dm, n) for p in a.problems for dc in a.discr
             for o in a.order for dm in a.dim for n in a.n]
    with open(os.path.join(a.root, "cases.txt"), "w") as f:
        f.write("\n".join(cases) + "\n")
    print("%d cases in %s" % (len(cases), a.root))


if __name__ == "__main__":
    main()
