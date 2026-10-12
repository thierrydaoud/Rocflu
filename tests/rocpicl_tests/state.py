"""Extract a comparable state (particles, gas fields, totals) from a run directory."""

import glob
import math
import os
import re

from .case import CASENAME
from .outputs import particle_ids, read_flow_ascii, read_particle_vtu

# Particle fields kept in goldens: full state vector plus the interpolated
# gas quantities. Other rprop entries are scratch and not always set.
PARTICLE_FIELDS = (
    ["x", "y", "z"]
    + ["y%02d" % j for j in range(4, 8)]           # velocity, temperature
    + ["rpro%02d" % j for j in (1, 2, 3, 4, 5)]    # rhop, rhof, dp, volp, phip
    + ["rpro%02d" % j for j in range(6, 13)]       # ux uy uz cs dpdx dpdy dpdz
    + ["rpro22"]                                   # superparticle loading
)
CV_NAMES = ["rho", "rhou", "rhov", "rhow", "rhoE"]


def _real_cells(rundir, region):
    path = os.path.join(rundir, "GRIDDATA", "%s.dim_%05d" % (CASENAME, region))
    with open(path) as f:
        text = f.read()
    return int(re.search(r"# Cells\s*\n\s*(\d+)", text).group(1))


def _dumps(rundir):
    """Sorted list of (time, SOL dir)."""
    sols = glob.glob(os.path.join(rundir, "SOL_*"))
    return sorted((float(os.path.basename(p)[4:]), p) for p in sols)


def _particle_files(rundir):
    return sorted(glob.glob(os.path.join(rundir, "par[0-9]*.vtu")))


def snapshot(rundir, case, which=-1):
    """State at dump index `which` (0 = initial, -1 = last).

    Gas and particle dumps are written together, so the n-th SOL_* directory
    and the n-th par*.vtu file share a time.
    """
    dumps = _dumps(rundir)
    pfiles = _particle_files(rundir)
    _, soldir = dumps[which]

    gas = {}
    t_gas = None
    for path in sorted(glob.glob(os.path.join(soldir, "%s.mixt.cva_*" % CASENAME))):
        region = int(os.path.basename(path).split("_")[1])
        # Full-precision time from the file; the directory name has 6 digits
        t_gas, arrays = read_flow_ascii(path)
        n = _real_cells(rundir, region)
        gas[region] = {name: arr[:n] for name, arr in zip(CV_NAMES, arrays)}

    particles = {}
    t_par = None
    if case.particles and pfiles:
        t_par, f = read_particle_vtu(pfiles[which])
        ids = particle_ids(f, len(case.particles), case.nranks)
        for k, pid in enumerate(ids):
            particles[pid] = {name: f[name][k] for name in PARTICLE_FIELDS if name in f}

    # Particle files store the time with 7 significant digits
    if t_par is not None and abs(t_par - t_gas) > 1e-5 * max(abs(t_gas), 1e-30):
        raise RuntimeError("Gas dump t=%g and particle dump t=%g do not match" % (t_gas, t_par))

    return {"time": t_gas, "gas": gas, "particles": particles,
            "totals": totals(gas, particles, case.grid.cell_volume)}


def totals(gas, particles, cell_volume):
    """Domain integrals of gas mass, momentum, energy and particle momentum."""
    out = {}
    for name in CV_NAMES:
        out["gas_" + name] = sum(sum(g[name]) for g in gas.values()) * cell_volume
    for c, v in zip("xyz", ("y04", "y05", "y06")):
        out["par_mom" + c] = sum(
            p["rpro22"] * p["rpro01"] * p["rpro04"] * p[v] for p in particles.values())
    out["par_ke"] = sum(
        0.5 * p["rpro22"] * p["rpro01"] * p["rpro04"]
        * (p["y04"] ** 2 + p["y05"] ** 2 + p["y06"] ** 2) for p in particles.values())
    for c, cv in zip("xyz", ("rhou", "rhov", "rhow")):
        out["mom" + c] = out["gas_" + cv] + out["par_mom" + c]
    return out


def all_finite(state):
    vals = [v for g in state["gas"].values() for arr in g.values() for v in arr]
    vals += [v for p in state["particles"].values() for v in p.values()]
    return all(math.isfinite(v) for v in vals)
