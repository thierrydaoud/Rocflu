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
    # Quasi-steady models other than the default Osnes (QSFORCE 3, MFC numbering)
    "qs-parmar": {("PICL", "QSFORCE"): 2},
    "qs-gidaspow": {("PICL", "QSFORCE"): 1},
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
        same_as=None if n == 1 else "box/momentum-exchange/np1")

# -- Periodic box: fixed particles in gas at rest -----------------------------
# Virtual cells at region borders used to keep the gas-phase state after the
# start-up conversion to phi_g-weighted variables, so the first stage leaked
# mass and set the gas moving (~5e-5 m/s), on 2+ ranks only.
for n in (2, 4):
    add(name="box/quiescent-cloud/np%d" % n,
        description="Fixed, unevenly spaced particles in gas at rest at uniform "
                    "pressure; the gas must stay at rest",
        nranks=n, grid=BOX, bc="periodic",
        inp={**uniform_gas(0.0), **TWO_WAY, ("PICL", "STATIONARY"): 1,
             ("PICL", "PGFORCE"): 1},
        particles=cloud(30, 3.0e-4, lo=0.0021, hi=0.0045),
        checks=["quiescent", "mass"])

# -- Periodic box: particles at rest exchanging heat with gas at rest ----------
# Particles start at 296.4 K (template TEMP), the gas at 101325 Pa and 1.2 kg/m^3
# (294.2 K). Before 2026-10 particles slipping < 1 mm/s skipped every force,
# heat transfer included, so their temperature never changed.
add(name="box/resting-heat-exchange/np1",
    description="Particles at rest in gas at rest at a different temperature; heat "
                "must flow and total energy (gas + particle thermal) be conserved",
    nranks=1, grid=BOX, bc="periodic",
    inp={**uniform_gas(0.0), **TWO_WAY, ("PICL", "HEATTRANSFER"): 2},
    particles=cloud(40, 2.0e-5),  # small: thermal time ~ dp^2, so T changes measurably
    checks=["particles_heat_exchange", "mass", "energy"])

# -- Closed box: interpolation of a linear gas field --------------------------
# Gas at rest at uniform pressure with density linear in x, y, z (a steady
# state). Particles at rest anywhere at least half a cell from the walls
# (inside the hull of cell centres) must see the exact density. The old
# inverse-distance (1/d^3) interpolation was off by up to ~15% of the
# density change across a cell, depending on where the particle sat in it.
LINEAR_RHO = {("INITFLOW", "IVAL1"): 1, ("INITFLOW", "RVAL7"): "25.0",
              ("INITFLOW", "RVAL8"): "-15.0", ("INITFLOW", "RVAL9"): "10.0",
              ("PICL", "PERIODICX"): 0, ("PICL", "PERIODICY"): 0,
              ("PICL", "PERIODICZ"): 0, ("PICL", "STATIONARY"): 1}
for n in (1, 2, 4):
    add(name="box/linear-interp/np%d" % n,
        description="Gas density linear in x, y, z at rest in a slip-wall box; "
                    "particles at rest must interpolate it exactly",
        nranks=n, grid=BOX, bc="closed",
        inp={**uniform_gas(0.0), **ONE_WAY, **LINEAR_RHO},
        particles=cloud(40, 1.0e-6, lo=0.0006, hi=0.0074),
        checks=["linear_interp", "mass"])

# -- Closed box of 2:1:1 cells: shape of the projected volume fraction --------
# One particle at rest, at least two cells from every wall. The old kernel
# sized each cell's Gaussian by V^(1/3), so on elongated cells the particle
# hardly spread along the long side (0.17 of a cell here) and spread too far
# across it (0.57). Scaling by each direction's cell length, as MFC does,
# gives ~0.44 everywhere.
ELONGATED = Grid(lx=0.016, ly=0.008, lz=0.008, x0=0.0, y0=0.0, z0=0.0, nx=8, ny=8, nz=8)
for n in (1, 2):
    add(name="box/elongated-projection/np%d" % n,
        description="One particle at rest in 2:1:1 cells; its volume fraction must be "
                    "centred on it and spread ~0.44 cell lengths in every direction",
        nranks=n, grid=ELONGATED, bc="closed",
        inp={**uniform_gas(0.0), **ONE_WAY, ("PICL", "PERIODICX"): 0,
             ("PICL", "PERIODICY"): 0, ("PICL", "PERIODICZ"): 0,
             ("PICL", "STATIONARY"): 1},
        particles=[(0.0066, 0.0044, 0.0037, 1.0e-4)],
        checks=["projection_shape", "mass"])

# -- Periodic box: particle collisions through ghost particles ----------------
# Pairs of equal particles start at rest, overlapping by 5% of a diameter, with
# every force off except collisions. The contact spring pushes each pair apart;
# its damping uses the partner's velocity. Every pair has the same offset, so
# all pairs must move identically, whether the partner is a real particle, a
# ghost on another rank, or a periodic image. On 1 rank the bin covers the whole
# box: the reference pair and the pairs at the 4 mm planes are real-real, and the
# pairs across the periodic faces use periodic ghosts. On 2 and 4 ranks the bin
# faces sit at 2, 4 or 6 mm, and the pairs across them use ghosts from another rank.
GHOST_DP = 1.0e-4
GHOST_OFFSET = (85.0e-6, 35.0e-6, 25.0e-6)  # |d| = 95.3 um < dp
GHOST_PAIR_CENTRES = [  # mm
    (3.0, 5.0, 3.0),     # reference: away from every bin face
    (4.0, 1.5, 5.5),     # across x = 4 mm
    (1.5, 4.0, 6.5),     # across y = 4 mm
    (6.5, 6.5, 4.0),     # across z = 4 mm
    (4.0, 4.0, 4.0),     # across all three: edge and corner ghosts
    (0.0, 6.5, 1.5),     # across the periodic x faces
    (5.5, 0.0, 6.5),     # across the periodic y faces
    (1.5, 5.5, 0.0),     # across the periodic z faces
    (0.0, 0.0, 0.0),     # across the periodic corner
]


def collision_pairs():
    pts = []
    for c in GHOST_PAIR_CENTRES:
        for sign in (-0.5, 0.5):
            pts.append(tuple(round((c[l] * 1e-3 + sign * GHOST_OFFSET[l]) % BOX.lx, 9)
                             for l in range(3)) + (GHOST_DP,))
    return pts


COLLISIONS_ONLY = {("PICL", "QSFORCE"): 0, ("PICL", "COLLISIONS"): 1,
                   ("TIMESTEP", "MAXTIME"): "2.05E-05", ("TIMESTEP", "RUNTIME"): "2.05E-05",
                   ("TIMESTEP", "WRITIME"): "1.0E-05", ("TIMESTEP", "PRNTIME"): "1.0E-05"}
for n in (1, 2, 4):
    add(name="box/ghost-collisions/np%d" % n,
        description="Overlapping pairs pushed apart by collisions only, across bin "
                    "faces, rank borders and periodic faces; every pair must move "
                    "like the reference pair and conserve its momentum",
        nranks=n, grid=BOX, bc="periodic",
        inp={**uniform_gas(0.0), **ONE_WAY, **COLLISIONS_ONLY},
        particles=collision_pairs(),
        checks=["particles_moved", "collision_pairs"],
        same_as=None if n == 1 else "box/ghost-collisions/np1")

# -- Periodic box: binary collisions with known outcomes ---------------------
# Four pairs approach each other (initial velocities in points.dat) with
# collisions as the only force (COLLISIONS 2: normal, tangential friction and
# spin). Each pair collides across a different boundary: a z bin face (2+
# ranks), a y bin face (4 ranks), the periodic x faces and the periodic y
# faces. Before 2026-10 the spring stiffness used particle i's own mass, so
# unequal pairs did not conserve momentum (an 8:1 pair at 20 m/s ended with
# its total momentum reversed).


def collision_pair(centre_mm, axis, dp1, dp2, speed, offset=0.0):
    """Two particles approaching along `axis` at +/-speed, 10 us from contact,
    with an impact offset (fraction of the mean diameter) along the next axis."""
    r = 0.5 * (dp1 + dp2)
    half = 0.5 * (r + 20e-6 * speed)
    side = (axis + 1) % 3
    pts = []
    for sign, dp in ((-1, dp1), (1, dp2)):
        xyz = [c * 1e-3 for c in centre_mm]
        xyz[axis] += sign * half
        xyz[side] += sign * 0.5 * offset * r
        vel = [0.0, 0.0, 0.0]
        vel[axis] = -sign * speed
        pts.append(tuple(round(v % BOX.lx, 9) for v in xyz) + (dp,) + tuple(vel))
    return pts


BINARY_COLLISIONS = (
    collision_pair((1.5, 1.5, 4.0), 2, 1.0e-4, 1.0e-4, 10.0)           # equal, z bin face
    + collision_pair((6.5, 4.0, 1.5), 1, 1.0e-4, 5.0e-5, 10.0)         # 8:1, y bin face
    + collision_pair((0.0, 6.5, 6.5), 0, 1.0e-4, 5.0e-5, 1.0)          # 8:1, periodic x
    + collision_pair((3.5, 0.0, 5.5), 1, 1.0e-4, 1.0e-4, 1.0, 0.5))    # oblique, periodic y
for n in (1, 2, 4):
    add(name="box/binary-collisions/np%d" % n,
        description="Four colliding pairs (equal and 8:1, head-on and oblique) across bin "
                    "faces and periodic faces; momentum, angular momentum and restitution",
        nranks=n, grid=BOX, bc="periodic",
        inp={**uniform_gas(0.0), **ONE_WAY, **COLLISIONS_ONLY, ("PICL", "COLLISIONS"): 2},
        particles=BINARY_COLLISIONS,
        checks=["particles_moved", "binary_collisions"],
        same_as=None if n == 1 else "box/binary-collisions/np1")

# -- Closed box: particle-wall collisions --------------------------------------
# The six slip walls are particle walls (PARTWALL 1 in the .bc file). Each
# particle starts 1 us from touching a wall and has collisions as its only
# force: faces at 1-300 m/s, a hit right on the diagonal that splits a quad
# face into two triangles, a box edge, a box corner and a 45-degree hit.
# Fast hits push the centre past the wall (soft contact lasting KSP time
# steps); the wall must push it back. Before 2026-10 particle walls existed
# only through filein.vtk, a particle on a shared triangle edge was pushed
# by both triangles, and one whose centre crossed a wall was pushed through.
WALL_T_HIT = 1.0e-6
WALL_SHOTS = [  # start point (mm, None = set from the speed), direction, speed
    ((None, 2.3, 5.6), (-1, 0, 0), 1.0),     # x- face, Hertz regime
    ((None, 6.3, 1.4), (-1, 0, 0), 10.0),    # x- face
    ((2.4, None, 6.6), (0, -1, 0), 100.0),   # y- face, centre crosses the wall
    ((4.6, 5.4, None), (0, 0, -1), 300.0),   # z- face on a quad diagonal
    ((None, 1.6, None), (1, 0, 1), 10.0),    # x+/z+ box edge
    ((None, None, None), (1, -1, -1), 10.0), # x+/y-/z- box corner
    ((5.6, 2.6, None), (1, 0, 1), 10.0),     # z+ face at 45 degrees
]


def wall_shots(dp=1.0e-4):
    """Particles touching the wall(s) they move towards after WALL_T_HIT."""
    r = 0.5 * dp
    pts = []
    for start, d, v in WALL_SHOTS:
        n = sum(c * c for c in d) ** 0.5
        u = [c / n for c in d]
        # a free coordinate is given; a wall coordinate is set from the speed
        walls = [k for k in range(3) if start[k] is None]
        xyz = []
        for k in range(3):
            if k in walls:
                gap = r + abs(u[k]) * v * WALL_T_HIT
                xyz.append(gap if d[k] < 0 else BOX.lx - gap)
            else:
                xyz.append(start[k] * 1e-3)
        pts.append(tuple(round(c, 9) for c in xyz) + (dp,) + tuple(c * v for c in u))
    return pts


for n in (1, 2, 4):
    add(name="box/walls/np%d" % n,
        description="Particles hitting the slip walls of a closed box (PARTWALL 1): faces, "
                    "a quad diagonal, an edge and a corner, 1-300 m/s; each rebounds with "
                    "the input restitution and ends inside the box",
        nranks=n, grid=BOX, bc="partwall",
        inp={**uniform_gas(0.0), **ONE_WAY, **COLLISIONS_ONLY, ("PICL", "COLLISIONS"): 2,
             ("PICL", "PERIODICX"): 0, ("PICL", "PERIODICY"): 0, ("PICL", "PERIODICZ"): 0,
             ("TIMESTEP", "MAXTIME"): "8.05E-06", ("TIMESTEP", "RUNTIME"): "8.05E-06",
             ("TIMESTEP", "WRITIME"): "4.0E-06", ("TIMESTEP", "PRNTIME"): "4.0E-06"},
        particles=wall_shots(),
        checks=["particles_moved", "wall_collisions"],
        same_as=None if n == 1 else "box/walls/np1",
        # The gas in this closed box differs slightly between decompositions
        # (~0.2% of a 4e-4 m/s flow driven by the particles' volume fraction
        # with FEEDBACK 0), so the gas seen by the particles (rprop) does too
        same_as_fields=("y",))

# -- Gas only: conservation of the flow solver itself --------------------------
JUMP = {("PICL", "USED"): 0, ("INITFLOW", "RVAL1"): "0.0043", ("INITFLOW", "RVAL3"): "0.0",
        ("INITFLOW", "RVAL5"): "0.97", ("INITFLOW", "RVAL6"): "82700.0"}
add(name="box/gas-jump/closed/np2",
    description="Pressure/density jump in a slip-wall box across two regions; "
                "mass and energy must be conserved to round-off",
    nranks=2, grid=BOX, bc="closed", inp={**uniform_gas(0.0), **JUMP},
    checks=["mass", "gas_energy"])
# Periodic virtual cells used to compute their own gradients, so the two sides
# of a periodic face disagreed (mass/energy drift ~4e-5 in 100 us; 11dc640)
for n in (1, 2):
    add(name="box/gas-jump/periodic/np%d" % n,
        description="Same jump in a fully periodic box; second-order fluxes across "
                    "periodic faces must conserve",
        nranks=n, grid=BOX, bc="periodic", inp={**uniform_gas(0.0), **JUMP},
        checks=["mass", "gas_energy"])


def by_name():
    return {c.name: c for c in CASES}
