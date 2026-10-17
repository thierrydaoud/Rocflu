"""Golden comparison and physics checks.

Goldens say whether a result changed. Physics checks say whether it is
right, and stay meaningful when a change (e.g. a new interpolation scheme)
moves every golden on purpose.
"""

import glob
import math
import os

from . import outputs, state
from .case import CASENAME
from .state import SPIN, VEL

# Gas fields are double precision; particle files are single precision
GAS_RTOL = 1e-10
PAR_RTOL = 2e-6


def _close(a, b, rtol, scale):
    return abs(a - b) <= rtol * max(abs(b), scale)


# Spin components are compared on the scale of the largest one, so a
# component that is only round-off (a head-on hit) is not judged on its own


def compare_states(new, ref, label="golden", fields=None):
    """Return a list of failure messages (empty if new matches ref).
    fields: only particle fields starting with one of these prefixes."""
    errors = []
    if not _close(new["time"], ref["time"], 1e-9, 0.0):
        errors.append("time %.10g != %s %.10g" % (new["time"], label, ref["time"]))

    if sorted(map(str, new["particles"])) != sorted(map(str, ref["particles"])):
        errors.append("particle ids differ from %s" % label)
    else:
        refp = {str(k): v for k, v in ref["particles"].items()}
        names = sorted({f for p in refp.values() for f in p})
        if fields:
            names = [f for f in names if f.startswith(tuple(fields))]
        for f in names:
            if any(f not in p for p in new["particles"].values()):
                errors.append("particle field %s is in %s but not written" % (f, label))
                continue
            group = SPIN if f in SPIN else (f,)
            scale = max(abs(p.get(g, 0.0)) for p in refp.values() for g in group) or 1.0
            worst = None
            for pid, p in new["particles"].items():
                a, b = p[f], refp[str(pid)][f]
                err = abs(a - b) / scale
                if not _close(a, b, PAR_RTOL, scale) and (worst is None or err > worst[0]):
                    worst = (err, pid, a, b)
            if worst:
                errors.append("particle %s %s: %.9g vs %s %.9g (rel. to field max %.2e)"
                              % (worst[1], f, worst[2], label, worst[3], worst[0]))

    if set(map(str, new["gas"])) != set(map(str, ref["gas"])):
        errors.append("gas regions differ from %s" % label)
    else:
        refg = {str(k): v for k, v in ref["gas"].items()}
        for region, g in sorted(new["gas"].items()):
            for name, arr in g.items():
                b_arr = refg[str(region)][name]
                if len(arr) != len(b_arr):
                    errors.append("region %s %s: %d cells vs %d" % (region, name, len(arr), len(b_arr)))
                    continue
                scale = max(abs(v) for v in b_arr) or 1.0
                bad = [(abs(a - b) / scale, i, a, b) for i, (a, b) in enumerate(zip(arr, b_arr))
                       if not _close(a, b, GAS_RTOL, scale)]
                if bad:
                    e, i, a, b = max(bad)
                    errors.append("region %s %s: %d cells differ, worst cell %d: %.15g vs %.15g"
                                  % (region, name, len(bad), i + 1, a, b))
    return errors


def check_restart_data(rundir):
    """Each SOL_* dump must have a fluid dump (restart_data/) of the same time
    holding the same state, and rfluinit's step 0 dump must exist.

    The dump holds the state in memory; with particles that is volume-weighted
    ((1-phi)*rho, ...), while the SOL_* files hold rho_g = rho/(1-phi), u, v, w,
    E (see state.FLOW_NAMES). rfluinit's step 0 dump holds the pure gas state."""
    errors = []
    picl = state._picl_used(rundir)
    dumps = {}
    for (step, region), path in outputs.fluid_dumps(rundir).items():
        dumps.setdefault(region, []).append((step, path))
    if not dumps:
        return ["restart_data: no fluid dumps"]
    for region, items in sorted(dumps.items()):
        if 0 not in [step for step, _ in items]:
            errors.append("restart_data: no step 0 dump of region %d" % region)
        dumps[region] = [outputs.read_fluid_dump(path) for _, path in sorted(items)]

    for _, soldir in state._dumps(rundir):
        for path in glob.glob(os.path.join(soldir, "%s.mixt.cva_*" % CASENAME)):
            region = int(os.path.basename(path).split("_")[1])
            # The directory name has 6 digits; the file holds the full time
            time, sol = outputs.read_flow_ascii(path)
            match = [d for d in dumps.get(region, []) if _close(d["time"], time, 1e-12, 1e-30)]
            if not match:
                errors.append("restart_data: no dump of region %d at t = %.17g" % (region, time))
                continue
            d = match[-1]
            n = len(d["cv"][0])
            cv = d["cv"]
            if picl and d["form"] == outputs.FORM_VOLUME_WEIGHTED:
                expected = [[cv[0][i] / (1.0 - d["phi"][i]) for i in range(n)]]
            else:
                expected = [list(cv[0])]
            if picl:
                expected += [[cv[k][i] / cv[0][i] for i in range(n)] for k in range(1, len(cv))]
            else:
                expected += [list(cv[k]) for k in range(1, len(cv))]
            for k, (e, s) in enumerate(zip(expected, sol)):
                scale = max(abs(v) for v in e)
                bad = [i for i in range(n) if not _close(s[i], e[i], 1e-12, scale)]
                if bad:
                    i = bad[0]
                    errors.append("restart_data: step %d region %d var %d: %d cells differ "
                                  "from SOL_ (cell %d: %.17g vs %.17g)"
                                  % (d["step"], region, k + 1, len(bad), i + 1, e[i], s[i]))
    return errors


def compare_particles_only(new, ref, label, fields=None):
    """Decomposition-independence: particles must agree, gas regions differ."""
    return compare_states(dict(new, gas={}), dict(ref, gas={}), label, fields)


# ---------------------------------------------------------------------------
# Physics checks: fn(initial_state, final_state, case) -> list of failures
# ---------------------------------------------------------------------------

def check_quiescent(init, final, case):
    """Gas at rest must stay at rest (max |u| below 1e-6 m/s)."""
    from .state import velocity
    picl = all("phi" in g for g in final["gas"].values())
    worst = 0.0
    for g in final["gas"].values():
        for i in range(len(g["q1"])):
            u, v, w = velocity(g, i, picl)
            worst = max(worst, math.sqrt(u * u + v * v + w * w))
    return [] if worst < 1e-6 else ["gas velocity %.3e m/s in a quiescent case" % worst]


def _conserved(init, final, keys, rtol, what):
    errors = []
    for k in keys:
        a, b = init["totals"][k], final["totals"][k]
        scale = max(abs(a), abs(b))
        if abs(b - a) > rtol * scale:
            errors.append("%s %s: %.12g -> %.12g (rel. change %.2e)"
                          % (what, k, a, b, abs(b - a) / scale))
    return errors


# Gas conservation tolerance. Round-off is ~1e-15 on one region; across
# region borders energy drifts ~5e-12 (box/gas-jump/closed/np2), small
# but above round-off. 1e-10 still flags the periodic-patch drift (~4e-5).
GAS_CONSERVATION_RTOL = 1e-10


def check_mass(init, final, case):
    return _conserved(init, final, ["gas_rho"], GAS_CONSERVATION_RTOL, "mass not conserved:")


# Exchange checks compare what the gas lost with what the particles gained.
# Particle data is single precision, so allow 0.1% of the exchanged amount.
EXCHANGE_RTOL = 1e-3


def _exchange(init, final, gas_key, par_key, what, scale=None):
    d_gas = final["totals"][gas_key] - init["totals"][gas_key]
    d_par = final["totals"][par_key] - init["totals"][par_key]
    scale = abs(d_par) if scale is None else scale
    if scale == 0.0:
        return []
    if abs(d_gas + d_par) > EXCHANGE_RTOL * scale:
        return ["%s not conserved: gas change %.6e, particle change %.6e (imbalance %.2e of exchange)"
                % (what, d_gas, d_par, abs(d_gas + d_par) / scale)]
    return []


def check_momentum(init, final, case):
    """Gas + particle momentum in a periodic box (no walls, no inflow).

    Each component is judged against the magnitude of the particle momentum
    change vector, so near-zero transverse components don't amplify noise.
    """
    d = [final["totals"]["par_mom" + c] - init["totals"]["par_mom" + c] for c in "xyz"]
    scale = math.sqrt(sum(v * v for v in d))
    errors = []
    for c, cv in zip("xyz", ("rhou", "rhov", "rhow")):
        errors += _exchange(init, final, "gas_" + cv, "par_mom" + c, "momentum " + c, scale)
    return errors


def check_energy(init, final, case):
    """Gas total energy + particle kinetic and thermal energy."""
    for st in (init, final):
        st["totals"]["par_energy"] = st["totals"]["par_ke"] + st["totals"].get("par_eth", 0.0)
    return _exchange(init, final, "gas_rhoE", "par_energy", "energy")


def check_gas_energy(init, final, case):
    """Gas total energy in a closed or periodic box without particles."""
    return _conserved(init, final, ["gas_rhoE"], GAS_CONSERVATION_RTOL, "energy not conserved:")


def check_particles_heat_exchange(init, final, case):
    """Particles at rest in gas at a different temperature exchange heat."""
    t0 = {k: p["Tp"] for k, p in init["particles"].items()}
    changed = any(abs(p["Tp"] - t0[k]) > 1e-6 * abs(t0[k]) for k, p in final["particles"].items())
    return [] if changed else ["particle temperatures never changed; resting particles exchange no heat"]


def check_particles_moved(init, final, case):
    """Guards conservation checks against a run where nothing happened."""
    moved = any(abs(p["vp_x"]) > 0 for p in final["particles"].values())
    return [] if moved else ["particles never moved; the case tests nothing"]


# Particle data is single precision (relative resolution ~6e-8).
LINEAR_INTERP_RTOL = 1e-6


def check_linear_interp(init, final, case):
    """Gas density linear in x, y, z (INITFLOW IVAL1 = 1) must be interpolated
    exactly at every particle: rhof = RVAL2*(1 + RVAL7*x + RVAL8*y + RVAL9*z)."""
    rho0 = float(case.inp[("INITFLOW", "RVAL2")])
    a, b, c = (float(case.inp[("INITFLOW", k)]) for k in ("RVAL7", "RVAL8", "RVAL9"))
    errors = []
    for pid, p in sorted(final["particles"].items()):
        exact = rho0 * (1.0 + a * p["x"] + b * p["y"] + c * p["z"])
        if abs(p["rho_g"] - exact) > LINEAR_INTERP_RTOL * exact:
            errors.append("particle %s at (%.4g, %.4g, %.4g): rhof %.9g, exact %.9g (rel. error %.2e)"
                          % (pid, p["x"], p["y"], p["z"], p["rho_g"], exact,
                             abs(p["rho_g"] - exact) / exact))
    return errors


# Gaussian projection on hex cells (MFC kernel, exp(-2.305*sum (dx_l/L_l)^2)
# over the 3x3x3 neighbours): the deposit's centroid is within ~0.04 cell of
# the particle and its spread is ~0.41-0.46 cell in every direction. The old
# kernel, sized by V^(1/3), gave 0.08-0.17 along the long side of 2:1:1 cells
# and 0.57 across it.
PROJ_CENTROID_TOL = 0.06
PROJ_SPREAD_RANGE = (0.35, 0.50)


def check_projection_shape(init, final, case):
    """One particle's volume fraction, in cell lengths per direction: its
    centroid must sit on the particle and its spread match the cell size."""
    (_, p), = final["particles"].items()
    xp = (p["x"], p["y"], p["z"])
    g = case.grid
    cell = (g.lx / g.nx, g.ly / g.ny, g.lz / g.nz)
    w, xs = [], []
    for region, gas in final["gas"].items():
        for phi, xc in zip(gas["phi"], final["_centroids"][region]):
            w.append(phi)
            xs.append(xc)
    total = sum(w)
    errors = []
    for l, axis in enumerate("xyz"):
        c = sum(wi * x[l] for wi, x in zip(w, xs)) / total
        s = math.sqrt(sum(wi * (x[l] - c) ** 2 for wi, x in zip(w, xs)) / total)
        err, spread = (c - xp[l]) / cell[l], s / cell[l]
        if abs(err) > PROJ_CENTROID_TOL:
            errors.append("%s: volume-fraction centroid off the particle by %.3f cells"
                          % (axis, err))
        if not PROJ_SPREAD_RANGE[0] <= spread <= PROJ_SPREAD_RANGE[1]:
            errors.append("%s: volume-fraction spread %.3f cells, expected %.2f-%.2f"
                          % ((axis, spread) + PROJ_SPREAD_RANGE))
    return errors


# Particle velocities are single precision; positions are ~4e-3 m in single
# precision (resolution ~5e-10 m), against a pair separation of ~1e-4 m.
PAIR_VEL_RTOL = 1e-5
PAIR_POS_ATOL = 2e-9


def check_collision_pairs(init, final, case):
    """Particles come in pairs (ids 2k-1, 2k) with the same offset, at rest and
    overlapping. Collisions alone must push every pair apart exactly like the
    first pair, conserving the pair's momentum, whether the partner is real,
    a ghost from another rank, or a periodic image."""
    g = case.grid
    length = (g.lx, g.ly, g.lz)
    p = {int(k): v for k, v in final["particles"].items()}
    pairs = []
    for k in range(1, len(p) // 2 + 1):
        a, b = p[2 * k - 1], p[2 * k]
        d = []
        for l, axis in enumerate("xyz"):  # separation, unwrapped across periodic faces
            dl = b[axis] - a[axis]
            d.append(dl - length[l] * round(dl / length[l]))
        va = [a[c] for c in VEL]
        vb = [b[c] for c in VEL]
        pairs.append((k, d, va, vb, a["dp"]))

    errors = []
    vscale = max(abs(v) for _, _, va, vb, _ in pairs for v in va + vb)
    if vscale == 0.0:
        return ["no pair moved; collisions never acted"]
    _, d1, va1, vb1, _ = pairs[0]
    for k, d, va, vb, dp in pairs:
        mom = max(abs(x + y) for x, y in zip(va, vb))
        if mom > PAIR_VEL_RTOL * vscale:
            errors.append("pair %d: momentum not conserved, |va + vb| = %.3e m/s (%.2e of max |v|)"
                          % (k, mom, mom / vscale))
        dist = math.sqrt(sum(x * x for x in d))
        if dist <= dp or sum(x * (y - z) for x, y, z in zip(d, vb, va)) <= 0.0:
            errors.append("pair %d: still in contact or approaching (distance %.4g dp)"
                          % (k, dist / dp))
        if k == 1:
            continue
        dv = max(abs(x - y) for x, y in zip(va + vb, va1 + vb1))
        dd = max(abs(x - y) for x, y in zip(d, d1))
        if dv > PAIR_VEL_RTOL * vscale or dd > PAIR_POS_ATOL:
            errors.append("pair %d moved differently from pair 1: velocity differs by %.2e of "
                          "max |v|, separation by %.3e m" % (k, dv / vscale, dd))
    return errors


# Binary collisions (pairs of ids 2k-1, 2k). Momentum is exact to single
# precision. Angular momentum about the pair's centre of mass changes by
# ~0.1% in an oblique contact: the friction torque acts at r_i from each
# centre, not at the overlapped contact point. Head-on pairs must rebound
# with the input restitution coefficient; with KSP = 10 steps per contact
# the measured e is within 0.5% (0.697 for EREST 0.7).
COLL_MOM_RTOL = 1e-5
COLL_ANGMOM_RTOL = 5e-3
COLL_E_RTOL = 1e-2


def _min_image(d, length):
    return d - length * round(d / length)


def check_binary_collisions(init, final, case):
    """Each pair conserves momentum and angular momentum and separates; a
    head-on pair rebounds with the input coefficient of restitution."""
    g = case.grid
    length = (g.lx, g.ly, g.lz)
    erest = float(case.inp.get(("PICL", "EREST"), 0.7))
    p0 = {int(k): v for k, v in init["particles"].items()}
    p1 = {int(k): v for k, v in final["particles"].items()}
    errors = []
    for k in range(1, len(p0) // 2 + 1):
        a, b = 2 * k - 1, 2 * k
        m = {i: p0[i]["rhop"] * p0[i]["volp"] for i in (a, b)}
        inertia = {i: m[i] * p0[i]["dp"] ** 2 / 10.0 for i in (a, b)}

        def vel(p, i):
            return [p[i][c] for c in VEL]

        def sep(p):
            return [_min_image(p[b][ax] - p[a][ax], length[l]) for l, ax in enumerate("xyz")]

        def momentum(p):
            return [m[a] * va + m[b] * vb for va, vb in zip(vel(p, a), vel(p, b))]

        def ang_mom(p):
            # about the centre of mass: reduced-mass orbital term plus spins
            d = sep(p)
            dv = [vb - va for va, vb in zip(vel(p, a), vel(p, b))]
            mu = m[a] * m[b] / (m[a] + m[b])
            orb = [mu * (d[1] * dv[2] - d[2] * dv[1]), mu * (d[2] * dv[0] - d[0] * dv[2]),
                   mu * (d[0] * dv[1] - d[1] * dv[0])]
            return [o + sum(inertia[i] * p[i][SPIN[l]] for i in (a, b))
                    for l, o in enumerate(orb)]

        scale = sum(m[i] * sum(v * v for v in vel(p0, i)) ** 0.5 for i in (a, b))
        dmom = max(abs(x - y) for x, y in zip(momentum(p1), momentum(p0)))
        if dmom > COLL_MOM_RTOL * scale:
            errors.append("pair %d: momentum changed by %.2e of sum |m v|" % (k, dmom / scale))

        d0, d1 = sep(p0), sep(p1)
        dist0 = math.sqrt(sum(x * x for x in d0))
        n0 = [x / dist0 for x in d0]
        dv0 = [vb - va for va, vb in zip(vel(p0, a), vel(p0, b))]
        dv1 = [vb - va for va, vb in zip(vel(p1, a), vel(p1, b))]
        lscale = (m[a] * m[b] / (m[a] + m[b])) * dist0 * math.sqrt(sum(x * x for x in dv0))
        dl = max(abs(x - y) for x, y in zip(ang_mom(p1), ang_mom(p0)))
        if dl > COLL_ANGMOM_RTOL * lscale:
            errors.append("pair %d: angular momentum changed by %.2e of mu*d*|dv|" % (k, dl / lscale))

        dist1 = math.sqrt(sum(x * x for x in d1))
        dp_mean = 0.5 * (p0[a]["dp"] + p0[b]["dp"])
        vn1 = sum(x * y for x, y in zip(dv1, d1)) / dist1
        if dist1 <= dp_mean or vn1 <= 0.0:
            errors.append("pair %d: never collided or still in contact (distance %.3g dp, "
                          "normal speed %.3g m/s)" % (k, dist1 / dp_mean, vn1))
            continue
        vn0 = sum(x * y for x, y in zip(dv0, n0))
        head_on = math.sqrt(max(0.0, sum(x * x for x in dv0) - vn0 * vn0)) < 1e-6 * abs(vn0)
        if head_on:
            e = vn1 / -vn0
            if abs(e - erest) > COLL_E_RTOL * erest:
                errors.append("pair %d: restitution %.5f, input %.3f" % (k, e, erest))
    return errors


WALL_E_TOL = 0.035


def check_wall_collisions(init, final, case):
    """Every particle that reaches a wall of the box rebounds from it and ends
    inside the box, clear of the wall. A particle that hits a single wall
    rebounds with normal restitution within WALL_E_TOL of the input; at an
    edge or corner each wall's friction (COLLISIONS >= 2) also slows the
    motion towards the other walls, so only the rebound is checked there."""
    g = case.grid
    lo = (g.x0, g.y0, g.z0)
    hi = (g.x0 + g.lx, g.y0 + g.ly, g.z0 + g.lz)
    erest = float(case.inp.get(("PICL", "EREST"), 0.7))
    t_end = float(case.inp[("TIMESTEP", "MAXTIME")])
    p1 = final["particles"]
    errors = []
    for pid, p in sorted(init["particles"].items(), key=lambda kv: int(kv[0])):
        r = 0.5 * p["dp"]
        q = p1.get(pid)
        if q is None:
            errors.append("particle %s lost" % pid)
            continue
        hits = []
        for k, ax in enumerate("xyz"):
            v0 = p[VEL[k]]
            gap0 = (hi[k] - p[ax]) if v0 > 0 else (p[ax] - lo[k])
            if v0 != 0.0 and gap0 - r <= abs(v0) * t_end:
                hits.append(k)
        for k in hits:
            ax = "xyz"[k]
            v0 = p[VEL[k]]
            e = -q[VEL[k]] / v0
            if e <= 0.0 or (len(hits) == 1 and abs(e - erest) > WALL_E_TOL):
                errors.append("particle %s, %s wall: restitution %.4f, input %.3f"
                              % (pid, ax, e, erest))
            gap1 = (hi[k] - q[ax]) if v0 > 0 else (q[ax] - lo[k])
            if gap1 < (1.0 - 1e-3) * r:
                errors.append("particle %s ends %.3g radii from the %s wall (outside if < 0)"
                              % (pid, gap1 / r, ax))
    return errors


def check_fluctuations_active(init, final, case):
    """Every particle must feel a nonzero force fluctuation at the end; a zero
    noise amplitude (no neighbours, no velocity spread) would test nothing."""
    quiet = [pid for pid, p in final["particles"].items()
             if all(p.get("Ffluct_" + c, 0.0) == 0.0 for c in "xyz")]
    if quiet:
        return ["%d of %d particles have no force fluctuation (e.g. particle %s)"
                % (len(quiet), len(final["particles"]), quiet[0])]
    return []


PHYSICS = {
    "quiescent": check_quiescent,
    "mass": check_mass,
    "momentum": check_momentum,
    "energy": check_energy,
    "gas_energy": check_gas_energy,
    "particles_moved": check_particles_moved,
    "particles_heat_exchange": check_particles_heat_exchange,
    "linear_interp": check_linear_interp,
    "projection_shape": check_projection_shape,
    "collision_pairs": check_collision_pairs,
    "binary_collisions": check_binary_collisions,
    "wall_collisions": check_wall_collisions,
    "fluctuations_active": check_fluctuations_active,
}
