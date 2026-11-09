#!/usr/bin/env python3
"""Write a scaling case: the curtain deck (base/) on a Cartesian box.

The box has NX x NY x NZ cells of size DX and starts at x = -LX/4, so the
initial shock (x = -5 mm) is inside it. The run uses the deck's CFL time step
and stops at END: runs on the same cell size follow the same time steps, so a
common END gives them the same number of steps. No output is written before
the end. The grid answers go to mesh.in for tests/mesh_generator.
"""

import argparse
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from rocpicl_tests.case import edit_inp  # noqa: E402

CASENAME = "rectshktb"


def drop_probes(text):
    """Remove the PROBE section: its points block and the keys after it."""
    lines = text.split("\n")
    start = next(i for i, l in enumerate(lines) if l.split()[:2] == ["#", "PROBE"])
    ends = [i for i in range(start + 1, len(lines)) if lines[i].strip() == "#"]
    end = ends[1] if len(ends) > 1 and not any(
        lines[i].strip().startswith("#") for i in range(ends[0] + 1, ends[1])) else ends[0]
    return "\n".join(lines[:start] + lines[end + 1:])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("rundir")
    ap.add_argument("--cells", type=int, nargs=3, required=True, metavar=("NX", "NY", "NZ"))
    ap.add_argument("--dx", type=float, default=2.5e-4, help="cell size [m] (curtain: 0.25 mm)")
    ap.add_argument("--end", type=float, default=6.0e-6, help="simulated time [s] (~200 steps)")
    ap.add_argument("--print-time", type=float, default=6.0e-7, help="time between progress lines [s]")
    ap.add_argument("--picl", type=int, default=0, choices=(0, 1))
    args = ap.parse_args()

    nx, ny, nz = args.cells
    lx, ly, lz = nx * args.dx, ny * args.dx, nz * args.dx
    end = args.end

    with open(os.path.join(HERE, "base", CASENAME + ".inp")) as f:
        inp = f.read()
    inp = drop_probes(inp)
    inp = edit_inp(inp, {
        ("TIMESTEP", "MAXTIME"): "%.10E" % end,
        ("TIMESTEP", "RUNTIME"): "%.10E" % end,
        ("TIMESTEP", "WRITIME"): "%.10E" % (10 * end),
        ("TIMESTEP", "FOMWRITIME"): "%.10E" % (10 * end),
        ("TIMESTEP", "PRNTIME"): "%.10E" % args.print_time,
        ("PICL", "USED"): str(args.picl),
    })

    os.makedirs(args.rundir, exist_ok=True)
    with open(os.path.join(args.rundir, CASENAME + ".inp"), "w") as f:
        f.write(inp)
    shutil.copy(os.path.join(HERE, "base", CASENAME + ".bc"), args.rundir)
    with open(os.path.join(args.rundir, "mesh.in"), "w") as f:
        f.write("\n".join(str(v) for v in (
            lx, ly, lz, -lx / 4, 0.0, 0.0, nx + 1, ny + 1, nz + 1, 2)) + "\n")
    print("%s: %d x %d x %d = %.3gM cells, to t = %.4g s" % (
        args.rundir, nx, ny, nz, nx * ny * nz / 1e6, end))


if __name__ == "__main__":
    main()
