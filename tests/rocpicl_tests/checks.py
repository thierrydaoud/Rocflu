"""Golden comparison and physics checks.

Goldens say whether a result changed. Physics checks say whether it is
right, and stay meaningful when a change (e.g. a new interpolation scheme)
moves every golden on purpose.
"""

import math

# Gas fields are double precision; particle files are single precision
GAS_RTOL = 1e-10
PAR_RTOL = 2e-6


def _close(a, b, rtol, scale):
    return abs(a - b) <= rtol * max(abs(b), scale)


def compare_states(new, ref, label="golden"):
    """Return a list of failure messages (empty if new matches ref)."""
    errors = []
    if not _close(new["time"], ref["time"], 1e-9, 0.0):
        errors.append("time %.10g != %s %.10g" % (new["time"], label, ref["time"]))

    if sorted(map(str, new["particles"])) != sorted(map(str, ref["particles"])):
        errors.append("particle ids differ from %s" % label)
    else:
        refp = {str(k): v for k, v in ref["particles"].items()}
        fields = sorted({f for p in refp.values() for f in p})
        for f in fields:
            scale = max(abs(p[f]) for p in refp.values()) or 1.0
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


def compare_particles_only(new, ref, label):
    """Decomposition-independence: particles must agree, gas regions differ."""
    return compare_states(dict(new, gas={}), dict(ref, gas={}), label)


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
    t0 = {k: p["y07"] for k, p in init["particles"].items()}
    changed = any(abs(p["y07"] - t0[k]) > 1e-6 * abs(t0[k]) for k, p in final["particles"].items())
    return [] if changed else ["particle temperatures never changed; resting particles exchange no heat"]


def check_particles_moved(init, final, case):
    """Guards conservation checks against a run where nothing happened."""
    moved = any(abs(p["y04"]) > 0 for p in final["particles"].values())
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
        if abs(p["rpro02"] - exact) > LINEAR_INTERP_RTOL * exact:
            errors.append("particle %s at (%.4g, %.4g, %.4g): rhof %.9g, exact %.9g (rel. error %.2e)"
                          % (pid, p["x"], p["y"], p["z"], p["rpro02"], exact,
                             abs(p["rpro02"] - exact) / exact))
    return errors


PHYSICS = {
    "quiescent": check_quiescent,
    "mass": check_mass,
    "momentum": check_momentum,
    "energy": check_energy,
    "gas_energy": check_gas_energy,
    "particles_moved": check_particles_moved,
    "particles_heat_exchange": check_particles_heat_exchange,
    "linear_interp": check_linear_interp,
}
