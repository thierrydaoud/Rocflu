"""Extract a comparable state (particles, gas fields, totals) from a run directory."""

import glob
import math
import os
import re

from . import outputs
from .case import CASENAME
from .outputs import particle_ids

# Particle fields kept in goldens: the particle state plus the gas quantities
# interpolated to the particle. Named as in ppiclf_user_OutputFields; the
# template input turns on PARTICLE_WRT, GAS_WRT and GRADGAS_WRT so they are written.
VEL = ("vp_x", "vp_y", "vp_z")
SPIN = ("omp_x", "omp_y", "omp_z")
PARTICLE_FIELDS = (
    ["x", "y", "z"]
    + list(VEL) + ["Tp"] + list(SPIN)
    + ["rhop", "rho_g", "dp", "volp", "phi_p"]
    + ["ug_x", "ug_y", "ug_z", "c_g", "dpdx_g", "dpdy_g", "dpdz_g"]
    # Quasi-steady force fluctuation: written only with QSFLUCT and MODEL_WRT 1
    + ["Ffluct_x", "Ffluct_y", "Ffluct_z"]
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
    """Particle specific heat (SPECIFICHEAT in the # PICL section)."""
    with open(os.path.join(rundir, CASENAME + ".inp")) as f:
        in_picl = False
        for line in f:
            s = line.strip()
            if s.startswith("#"):
                in_picl = s[1:].split() == ["PICL"]
            elif in_picl and s.split()[:1] == ["SPECIFICHEAT"]:
                return float(s.split()[1])
    return 0.0


def _particle_files(rundir):
    return sorted(glob.glob(os.path.join(rundir, "par[0-9]*.vtu")))


def snapshot(rundir, case, which=-1):
    """State at dump step index `which` (0 = initial, -1 = last), from the dumps in
    restart_data/ (the old-style SOL_*/par*.vtu are checked against them when a
    case writes them).

    The gas arrays are in the form of the old-style flow files (FLOW_NAMES): with
    particles rho_g, u, v, w, E; without, the conserved state. phi with particles:
    from the fluid dump, or at step 0 (rfluinit's dump) the volume projected at the
    start of the run (particle dump) over the cell volume. Particle fields: those of
    the particle dump's output-field table, with x, y, z."""
    fluid = outputs.fluid_dumps(rundir)
    steps = sorted({s for s, _ in fluid})
    step = steps[which]
    picl = _picl_used(rundir)
    particle_dumps = outputs.particle_dumps(rundir)
    pdump = outputs.read_particle_dump(particle_dumps[step]) if picl and step in particle_dumps \
        else None

    gas = {}
    centroids = {}  # only for checks that need geometry; not written to goldens
    t_gas = None
    for region in sorted(r for s, r in fluid if s == step):
        d = outputs.read_fluid_dump(fluid[(step, region)])
        t_gas = d["time"]
        cv = d["cv"]
        n = len(cv[0])
        if picl:
            if d["form"] == outputs.FORM_VOLUME_WEIGHTED:
                phi = list(d["phi"])
                rho_g = [cv[0][i] / (1.0 - phi[i]) for i in range(n)]
            else:
                rank = max(region, 1) - 1
                phi = [v / case.grid.cell_volume for v in pdump["projected"][rank][0]] \
                    if pdump else [0.0] * n
                rho_g = list(cv[0])
            arrays = [rho_g] + [[cv[k][i] / cv[0][i] for i in range(n)] for k in range(1, 5)]
            gas[region] = dict(zip(FLOW_NAMES, arrays))
            gas[region]["phi"] = phi
        else:
            gas[region] = dict(zip(FLOW_NAMES, [list(a) for a in cv]))
        if "projection_shape" in case.checks:
            centroids[region] = outputs.read_post_vtu(
                os.path.join(rundir, "post", "p%d" % max(region - 1, 0), "%08d.vtu" % step),
                centroids=True)

    particles = {}
    if case.particles and pdump:
        names = {name: (source, slot) for source, slot, name in pdump["fields"]}
        arrays = {"y": pdump["y"], "rprop": pdump["rprop"], "iprop": pdump["iprop"]}
        source_name = {outputs.OUT_Y: "y", outputs.OUT_RPROP: "rprop", outputs.OUT_TAG: "iprop"}
        f = {name: [p[slot - 1] for p in arrays[source_name[source]]]
             for name, (source, slot) in names.items()}
        f.update({c: [p[k] for p in pdump["y"]] for k, c in enumerate("xyz")})
        ids = particle_ids(f, len(case.particles), case.nranks)
        for k, pid in enumerate(ids):
            particles[pid] = {name: f[name][k] for name in PARTICLE_FIELDS if name in f}

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
    # Particle mass rhop*volp: the superparticle loading is always 1 (not written)
    for c, v in zip("xyz", VEL):
        out["par_mom" + c] = sum(
            p["rhop"] * p["volp"] * p[v] for p in particles.values())
    out["par_ke"] = sum(
        0.5 * p["rhop"] * p["volp"]
        * (p["vp_x"] ** 2 + p["vp_y"] ** 2 + p["vp_z"] ** 2) for p in particles.values())
    # Particle thermal energy m*cp*T (temperature is y07)
    out["par_eth"] = sum(
        p["rhop"] * p["volp"] * cp_particle * p["Tp"]
        for p in particles.values())
    for c, cv in zip("xyz", ("rhou", "rhov", "rhow")):
        out["mom" + c] = out["gas_" + cv] + out["par_mom" + c]
    return out


def all_finite(state):
    vals = [v for g in state["gas"].values() for arr in g.values() for v in arr]
    vals += [v for p in state["particles"].values() for v in p.values()]
    return all(math.isfinite(v) for v in vals)
