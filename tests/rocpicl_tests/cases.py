"""The regression cases.

Particles never start on a cell face (faces sit on whole millimetres here),
so round-off cannot decide which cell or rank a particle starts in.
"""

from .case import Case, Grid

# ---------------------------------------------------------------------------
# Shared settings
# ---------------------------------------------------------------------------

# Shock tube: 40 x 4 x 4 cells of 1 mm, shock at x = -10 mm moving into
# particles at rest. Dumps at ~10 and ~20 us; compared at the last one.
SHOCKTUBE = {
    ("TIMESTEP", "MAXTIME"): "2.05E-05",
    ("TIMESTEP", "RUNTIME"): "2.05E-05",
    ("TIMESTEP", "WRITIME"): "1.0E-05",
    ("TIMESTEP", "PRNTIME"): "1.0E-05",
}
ONE_WAY = {("PICL", "FEEDBACK"): 0}
TWO_WAY = {("PICL", "FEEDBACK"): 1}

ONE_PARTICLE = [(-0.0093, 0.00113, 0.00137, 1.0e-4)]
FOUR_PARTICLES = [
    (-0.0093, 0.00113, 0.00137, 1.0e-4),
    (-0.0071, 0.00287, 0.00211, 1.0e-4),
    (-0.0052, 0.00169, 0.00263, 1.0e-4),
    (-0.0037, 0.00241, 0.00089, 1.0e-4),
]

# Periodic 8 mm cube, 8^3 cells, uniform gas (x < RVAL1 = whole box gets
# the "post-shock" state RVAL2-4)
BOX = Grid(lx=0.008, ly=0.008, lz=0.008, x0=0.0, y0=0.0, z0=0.0, nx=8, ny=8, nz=8)


def uniform_gas(u):
    return {
        ("INITFLOW", "RVAL1"): "1.0",
        ("INITFLOW", "RVAL2"): "1.2",
        ("INITFLOW", "RVAL3"): str(u),
        ("INITFLOW", "RVAL4"): "101325.0",
        ("INITFLOW", "RVAL5"): "1.2",
        ("INITFLOW", "RVAL6"): "101325.0",
        ("PICL", "PERIODICX"): 1,
        ("PICL", "PERIODICY"): 1,
        ("PICL", "PERIODICZ"): 1,
        ("TIMESTEP", "MAXTIME"): "1.05E-04",
        ("TIMESTEP", "RUNTIME"): "1.05E-04",
        ("TIMESTEP", "WRITIME"): "5.0E-05",
        ("TIMESTEP", "PRNTIME"): "5.0E-05",
    }


def cloud(n, dp, lo=0.0021, hi=0.0059, seed=12345):
    """n particles in [lo, hi]^3, off cell faces (whole mm), reproducible."""
    pts = []
    state = seed
    while len(pts) < n:
        xyz = []
        for _ in range(3):
            state = (1103515245 * state + 12345) % 2 ** 31
            xyz.append(lo + (hi - lo) * state / 2 ** 31)
        if all(0.1 < (v * 1e3) % 1.0 < 0.9 for v in xyz):
            pts.append(tuple(round(v, 7) for v in xyz) + (dp,))
    return pts


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

CASES = []


def add(**kw):
    CASES.append(Case(**kw))


# -- Shock tube, one-way coupling --------------------------------------------
for n in (1, 2):
    add(name="shocktube/one-way/1-particle/np%d" % n,
        description="Shock over one particle; ppiclF bins around a single particle",
        nranks=n, inp={**SHOCKTUBE, **ONE_WAY}, particles=ONE_PARTICLE,
        checks=["particles_moved"],
        same_as=None if n == 1 else "shocktube/one-way/1-particle/np1")
    add(name="shocktube/one-way/4-particles/np%d" % n,
        description="Shock over four particles, quasi-steady drag only",
        nranks=n, inp={**SHOCKTUBE, **ONE_WAY}, particles=FOUR_PARTICLES,
        checks=["particles_moved"],
        same_as=None if n == 1 else "shocktube/one-way/4-particles/np1")

# -- Shock tube, two-way coupling --------------------------------------------
for n in (1, 2):
    add(name="shocktube/two-way/4-particles/np%d" % n,
        description="Shock over four particles with feedback to the gas",
        nranks=n, inp={**SHOCKTUBE, **TWO_WAY}, particles=FOUR_PARTICLES,
        checks=["particles_moved"],
        same_as=None if n == 1 else "shocktube/two-way/4-particles/np1")

# -- Force models, one at a time on top of two-way QS (2 ranks) ---------------
FORCES = {
    "added-mass": {("PICL", "AMFORCE"): 2},
    "pressure-gradient": {("PICL", "PGFORCE"): 1},
    "viscous-unsteady": {("PICL", "VISCOUSUNSTEADY"): 1},
    "heat-transfer": {("PICL", "HEATTRANSFER"): 3},
    "qs-fluctuations": {("PICL", "QSFLUCT"): 2},
}
for force, opts in FORCES.items():
    add(name="shocktube/two-way/%s/np2" % force,
        description="Two-way QS drag plus %s" % force.replace("-", " "),
        nranks=2, inp={**SHOCKTUBE, **TWO_WAY, **opts}, particles=FOUR_PARTICLES,
        checks=["particles_moved"])

# -- Periodic box: gas moving through particles at rest -----------------------
for n in (1, 2):
    add(name="box/momentum-exchange/np%d" % n,
        description="Uniform gas at 50 m/s through 40 particles at rest in a periodic "
                    "box; gas + particle momentum and energy must be conserved",
        nranks=n, grid=BOX, bc="periodic",
        inp={**uniform_gas(50.0), **TWO_WAY}, particles=cloud(40, 2.0e-5),
        checks=["particles_moved", "mass", "momentum", "energy"],
        same_as=None if n == 1 else "box/momentum-exchange/np1",
        known_issue="Gas takes up ~83% of the particle momentum gain, independent of dt; "
                    "periodic BCs are not conservative (see box/gas-jump/periodic); "
                    "np1 vs np2 interpolated gradients differ by ~1%")

# -- Periodic box: fixed particles in gas at rest -----------------------------
add(name="box/quiescent-cloud/np2",
    description="Fixed, unevenly spaced particles in gas at rest at uniform pressure; "
                "the gas must stay at rest",
    nranks=2, grid=BOX, bc="periodic",
    inp={**uniform_gas(0.0), **TWO_WAY, ("PICL", "STATIONARY"): 1,
         ("PICL", "PGFORCE"): 1},
    particles=cloud(30, 3.0e-4, lo=0.0021, hi=0.0045),
    checks=["quiescent", "mass"],
    known_issue="Gas starts moving (~4e-5 m/s): cv is initialised without the particle "
                "volume fraction, but pressure is computed with 1/(1-phi)")

# -- Gas only: conservation of the flow solver itself --------------------------
JUMP = {("PICL", "USED"): 0, ("INITFLOW", "RVAL1"): "0.0043", ("INITFLOW", "RVAL3"): "0.0",
        ("INITFLOW", "RVAL5"): "0.97", ("INITFLOW", "RVAL6"): "82700.0"}
add(name="box/gas-jump/closed/np2",
    description="Pressure/density jump in a slip-wall box across two regions; "
                "mass and energy must be conserved to round-off",
    nranks=2, grid=BOX, bc="closed", inp={**uniform_gas(0.0), **JUMP},
    checks=["mass", "gas_energy"])
add(name="box/gas-jump/periodic/np1",
    description="Same jump in a fully periodic box",
    nranks=1, grid=BOX, bc="periodic", inp={**uniform_gas(0.0), **JUMP},
    checks=["mass", "gas_energy"],
    known_issue="Rocflu periodic patches are not conservative: mass and energy drift ~4e-5 "
                "in 100 us with no particles (closed box: <1e-11)")


def by_name():
    return {c.name: c for c in CASES}
