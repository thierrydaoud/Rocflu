"""Extract a comparable state (particles, gas fields, totals) from a run directory."""

import glob
import math
import os
import re

from .case import CASENAME
from .outputs import (particle_ids, read_flow_ascii, read_particle_vtu, read_vtu_cell_centroids,
                      read_vtu_cell_field)

# Particle fields kept in goldens: full state vector plus the interpolated
# gas quantities. Other rprop entries are scratch and not always set.
PARTICLE_FIELDS = (
    ["x", "y", "z"]
    + ["y%02d" % j for j in range(4, 11)]          # velocity, temperature, spin
    + ["rpro%02d" % j for j in (1, 2, 3, 4, 5)]    # rhop, rhof, dp, volp, phip
    + ["rpro%02d" % j for j in range(6, 13)]       # ux uy uz cs dpdx dpdy dpdz
    + ["rpro22"]                                   # superparticle loading
)
# What the five arrays of a Rocflu ASCII flow file hold depends on whether
# particles are on (RFLU_TimeStepping converts real cells just before writing):
#   PICL off: conserved bulk state      rho,  rho*u, rho*v, rho*w, rho*E
#   PICL on:  gas-phase primitive state rho_g, u, v, w, E   with rho_g = rho/(1-phi)
# Goldens store the arrays as written; totals are converted to bulk quantities.
FLOW_NAMES = ["q1", "q2", "q3", "q4", "q5"]


def _real_cells(rundir, region):
    path = os.path.join(rundir, "GRIDDATA", "%s.dim_%05d" % (CASENAME, region))
    with open(path) as f:
        text = f.read()
    return int(re.search(r"# Cells\s*\n\s*(\d+)", text).group(1))


def _dumps(rundir):
    """Sorted list of (time, SOL dir), skipping directories without flow files
    (the solver can create the directory of a final dump it never writes)."""
    sols = [p for p in glob.glob(os.path.join(rundir, "SOL_*"))
            if glob.glob(os.path.join(p, "%s.mixt.cva_*" % CASENAME))]
    return sorted((float(os.path.basename(p)[4:]), p) for p in sols)


def _picl_used(rundir):
    with open(os.path.join(rundir, CASENAME + ".inp")) as f:
        in_picl = False
        for line in f:
            s = line.strip()
            if s.startswith("#"):
                in_picl = s[1:].split() == ["PICL"]
            elif in_picl and s.split()[:1] == ["USED"]:
                return int(s.split()[1]) != 0
    return False


def _particle_cp(rundir):
    """Particle specific heat (SPHT in the # MATERIAL section)."""
    with open(os.path.join(rundir, CASENAME + ".inp")) as f:
        in_mat = False
        for line in f:
            s = line.strip()
            if s.startswith("#"):
                in_mat = s[1:].split() == ["MATERIAL"]
            elif in_mat and s.split()[:1] == ["SPHT"]:
                return float(s.split()[1])
    return 0.0


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

    picl = _picl_used(rundir)
    gas = {}
    centroids = {}  # only for checks that need geometry; not written to goldens
    t_gas = None
    for path in sorted(glob.glob(os.path.join(soldir, "%s.mixt.cva_*" % CASENAME))):
        region = int(os.path.basename(path).split("_")[1])
        # Full-precision time from the file; the directory name has 6 digits
        t_gas, arrays = read_flow_ascii(path)
        n = _real_cells(rundir, region)
        gas[region] = {name: arr[:n] for name, arr in zip(FLOW_NAMES, arrays)}
        if picl:
            # Gas volume fraction comes from the ParaView file of the same dump
            suffix = os.path.basename(soldir)[4:]
            vtu = os.path.join(rundir, "PARAVIEW_" + suffix,
                               "%s_%s_%05d.vtu" % (CASENAME, suffix, region))
            phi = read_vtu_cell_field(vtu, "Particle Volume Fraction")
            if len(phi) != n:
                raise RuntimeError("%s has %d cells, expected %d" % (vtu, len(phi), n))
            gas[region]["phi"] = phi
            if "projection_shape" in case.checks:
                centroids[region] = read_vtu_cell_centroids(vtu)

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

    snap = {"time": t_gas, "gas": gas, "particles": particles,
            "totals": totals(gas, particles, case.grid.cell_volume, picl,
                             _particle_cp(rundir))}
    if centroids:
        snap["_centroids"] = centroids
    return snap


def bulk_state(g, i, picl):
    """Bulk (mixture-volume) density, momentum and energy of cell i."""
    if picl:
        rho = (1.0 - g["phi"][i]) * g["q1"][i]
        return rho, rho * g["q2"][i], rho * g["q3"][i], rho * g["q4"][i], rho * g["q5"][i]
    return g["q1"][i], g["q2"][i], g["q3"][i], g["q4"][i], g["q5"][i]


def velocity(g, i, picl):
    if picl:
        return g["q2"][i], g["q3"][i], g["q4"][i]
    return g["q2"][i] / g["q1"][i], g["q3"][i] / g["q1"][i], g["q4"][i] / g["q1"][i]


def totals(gas, particles, cell_volume, picl, cp_particle=0.0):
    """Domain integrals of gas mass, momentum, energy and particle momentum."""
    sums = [0.0] * 5
    for g in gas.values():
        for i in range(len(g["q1"])):
            for k, val in enumerate(bulk_state(g, i, picl)):
                sums[k] += val
    out = {name: val * cell_volume for name, val in
           zip(("gas_rho", "gas_rhou", "gas_rhov", "gas_rhow", "gas_rhoE"), sums)}
    for c, v in zip("xyz", ("y04", "y05", "y06")):
        out["par_mom" + c] = sum(
            p["rpro22"] * p["rpro01"] * p["rpro04"] * p[v] for p in particles.values())
    out["par_ke"] = sum(
        0.5 * p["rpro22"] * p["rpro01"] * p["rpro04"]
        * (p["y04"] ** 2 + p["y05"] ** 2 + p["y06"] ** 2) for p in particles.values())
    # Particle thermal energy m*cp*T (temperature is y07)
    out["par_eth"] = sum(
        p["rpro22"] * p["rpro01"] * p["rpro04"] * cp_particle * p["y07"]
        for p in particles.values())
    for c, cv in zip("xyz", ("rhou", "rhov", "rhow")):
        out["mom" + c] = out["gas_" + cv] + out["par_mom" + c]
    return out


def all_finite(state):
    vals = [v for g in state["gas"].values() for arr in g.values() for v in arr]
    vals += [v for p in state["particles"].values() for v in p.values()]
    return all(math.isfinite(v) for v in vals)
