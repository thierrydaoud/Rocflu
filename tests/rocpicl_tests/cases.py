"""The regression cases.

Particles never start on a cell face (faces sit on whole millimetres here),
so round-off cannot decide which cell or rank a particle starts in.
"""

import math
import re

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

# -- Shock tube, pressure-gradient force alone --------------------------------
# The template's shock state, given here so the check can read it. At 40 us the
# shock is 8-10 cells past the particles (non-corner cell columns) and the
# start-up contact from x = -10 mm has not reached them.
SHOCK_STATE = {("INITFLOW", "RVAL2"): "1.4190765531700", ("INITFLOW", "RVAL3"): "138.9205016483263",
               ("INITFLOW", "RVAL4"): "141869.8040079374623", ("INITFLOW", "RVAL5"): "0.9700779339532",
               ("INITFLOW", "RVAL6"): "82700.0", ("PICL", "DENSITY"): "1880.0"}
add(name="shocktube/one-way/pg-impulse/np2",
    description="Pressure-gradient force alone on two particles the shock passes",
    nranks=2, inp={**SHOCKTUBE, **ONE_WAY, **SHOCK_STATE, ("PICL", "QSFORCE"): 0,
                   ("PICL", "PGFORCE"): 1, ("TIMESTEP", "MAXTIME"): "4.05E-05",
                   ("TIMESTEP", "RUNTIME"): "4.05E-05"},
    particles=[(-0.0030, 0.00163, 0.00137, 1.0e-4), (-0.0010, 0.00241, 0.00263, 1.0e-4)],
    checks=["pg_impulse"])

# -- Shock tube, two-way coupling --------------------------------------------
for n in (1, 2):
    add(name="shocktube/two-way/4-particles/np%d" % n,
        description="Shock over four particles with feedback to the gas",
        nranks=n, inp={**SHOCKTUBE, **TWO_WAY}, particles=FOUR_PARTICLES,
        checks=["particles_moved"],
        same_as=None if n == 1 else "shocktube/two-way/4-particles/np1")

# -- Force models, one at a time on top of two-way QS (2 ranks) ---------------
# Added mass (AMFORCE 2, Briney) is left out until it is validated; the case
# validator rejects it.
FORCES = {
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
        checks=["particles_moved"],
        # history of past steps and the NSCBC outflow must survive a restart
        restart_test=(force == "viscous-unsteady"))

# -- Shock tube, added mass in a dense slab ----------------------------------
# 288 particles of 0.4 mm fill x in [-5, -3] mm (phi_p ~0.3 in 32 cells). Added
# mass reads the gas Du/Dt, whose pressure part must be -grad(p)/rho_g; taken
# from the residual alone it is -grad(p)/(phi_g rho_g), 1.4x larger here.
def slab(n, dp, x=(-0.0049, -0.0031), yz=(0.0001, 0.0039), seed=2024):
    """n particles in x[0..1] by yz[0..1]^2, off cell faces (whole mm), reproducible."""
    pts = []
    state = seed
    while len(pts) < n:
        xyz = []
        for lo, hi in (x, yz, yz):
            state = (1103515245 * state + 12345) % 2 ** 31
            xyz.append(lo + (hi - lo) * state / 2 ** 31)
        if all(0.1 < (v * 1e3) % 1.0 < 0.9 for v in xyz):
            pts.append(tuple(round(v, 7) for v in xyz) + (dp,))
    return pts


add(name="shocktube/two-way/added-mass-dense/np2",
    description="Shock through a dense slab (phi_p ~0.3) with PGFORCE 1 and AMFORCE 1 "
                "(Parmar); added mass sees the undisturbed Du/Dt",
    nranks=2, inp={**SHOCKTUBE, **TWO_WAY, ("PICL", "PGFORCE"): 1, ("PICL", "AMFORCE"): 1},
    particles=slab(288, 4.0e-4), checks=["particles_moved"])

# -- Shock tube with slip side walls: constrained reconstruction --------------
# 8 x 8 cells across, so nearly every cell is within the three layers next to a
# wall where CRECONSTC 1 refits gradients. It used to zero the gradients of
# variables without a Dirichlet face there (e7daa9b): 6-9% in the gas at the
# shock and 5-21% in vp_x against CRECONSTC 0. vp_y, vp_z feel the wall-normal
# velocity constraint, so only x and vp_x must match.
SLIPWALL_TUBE = Grid(ly=0.008, lz=0.008, ny=8, nz=8)
for cr in (0, 1):
    add(name="shocktube/slip-walls/creconst%d/np2" % cr,
        description="Shock over four particles between slip walls, PGFORCE 1, "
                    "CRECONSTC %d; the shock must stay planar" % cr,
        nranks=2, grid=SLIPWALL_TUBE, bc="shocktube-walls",
        inp={**SHOCKTUBE, **TWO_WAY, ("PICL", "PGFORCE"): 1,
             ("NUMERICS", "CRECONSTC"): cr},
        particles=FOUR_PARTICLES, checks=["particles_moved", "planar"],
        same_as=None if cr == 0 else "shocktube/slip-walls/creconst0/np2",
        same_as_fields=("x", "vp_x"))

# -- Periodic box: gas moving through particles at rest -----------------------
for n in (1, 2):
    add(name="box/momentum-exchange/np%d" % n,
        description="Uniform gas at 50 m/s through 40 particles at rest in a periodic "
                    "box; gas + particle momentum and energy must be conserved",
        nranks=n, grid=BOX, bc="periodic",
        inp={**uniform_gas(50.0), **TWO_WAY}, particles=cloud(40, 2.0e-5),
        checks=["particles_moved", "mass", "momentum", "energy"],
        same_as=None if n == 1 else "box/momentum-exchange/np1",
        restart_test=(n == 2))

# -- Periodic box: Sun heat transfer in the dilute blend ----------------------
# Gas at 5 m/s (Re_p ~ 150) through 40 hot particles; cell volume fractions
# stay below 0.1 (max ~0.075), so Sun's model blends to Ranz-Marshall.
add(name="box/sun-dilute/np1",
    description="Gas at 5 m/s through 40 hot particles with Sun's heat transfer "
                "below eps_s = 0.1 (Ranz-Marshall blend); energy must be conserved",
    nranks=1, grid=BOX, bc="periodic",
    inp={**uniform_gas(5.0), **TWO_WAY, ("PICL", "HEATTRANSFER"): 5,
         ("PICL", "TEMP"): "394.17"},
    particles=cloud(40, 4.5e-4),
    checks=["particles_heat_exchange", "mass", "momentum", "energy"])

# -- Periodic box: quasi-steady force fluctuations (QSFLUCT 2, Osnes) ---------
# The noise amplitude grows with the granular temperature, the velocity spread
# among a particle's neighbours, so the particles start with a spread of up to
# 1 m/s per component and the neighbour width (3 mm) covers the cloud. Each
# particle draws from its own random stream, so 1 and 2 ranks must agree; with
# one generator per rank they did not.
def with_velocity_spread(particles, spread=1.0, seed=7):
    """Add reproducible velocities in [-spread, spread] to (x, y, z, dp) particles."""
    out = []
    state = seed
    for p in particles:
        uvw = []
        for _ in range(3):
            state = (1103515245 * state + 12345) % 2 ** 31
            uvw.append(spread * (2.0 * state / 2 ** 31 - 1.0))
        out.append(tuple(p[:4]) + tuple(uvw))
    return out


for n in (1, 2):
    add(name="box/qs-fluctuations/np%d" % n,
        description="Uniform gas at 50 m/s through 40 particles with a velocity spread; "
                    "the force fluctuation must be on and not depend on the rank count",
        nranks=n, grid=BOX, bc="periodic",
        inp={**uniform_gas(50.0), **TWO_WAY, ("PICL", "QSFLUCT"): 2,
             ("PICL", "MODEL_WRT"): 1, ("PICL", "NEIGHBORWIDTH"): "3.0e-3"},
        particles=with_velocity_spread(cloud(40, 2.0e-5)),
        checks=["particles_moved", "fluctuations_active"],
        same_as=None if n == 1 else "box/qs-fluctuations/np1",
        restart_test=(n == 2))

# -- Periodic box: fixed particles in gas at rest -----------------------------
# Virtual cells at region borders used to keep the gas-phase state after the
# start-up conversion to phi_g-weighted variables, so the first stage leaked
# mass and set the gas moving (~5e-5 m/s), on 2+ ranks only.
for n in (2, 4):
    add(name="box/quiescent-cloud/np%d" % n,
        description="Fixed, unevenly spaced particles in gas at rest at uniform "
                    "pressure; the gas must stay at rest",
        nranks=n, grid=BOX, bc="periodic",
        inp={**uniform_gas(0.0), **TWO_WAY, ("PICL", "VELOCITYMODEL"): 0,
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
              ("PICL", "VELOCITYMODEL"): 0}
for n in (1, 2, 4):
    add(name="box/linear-interp/np%d" % n,
        description="Gas density linear in x, y, z at rest in a slip-wall box; "
                    "particles at rest must interpolate it exactly",
        nranks=n, grid=BOX, bc="closed",
        inp={**uniform_gas(0.0), **ONE_WAY, **LINEAR_RHO},
        particles=cloud(40, 1.0e-6, lo=0.0006, hi=0.0074),
        checks=["linear_interp", "linear_grad_rhog", "mass"])

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
        inp={**uniform_gas(0.0), **ONE_WAY,
             ("PICL", "VELOCITYMODEL"): 0},
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
        same_as=None if n == 1 else "box/binary-collisions/np1",
        restart_test=(n == 2))

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


def wall_shots(dp=1.0e-4, shots=WALL_SHOTS):
    """Particles touching the wall(s) they move towards after WALL_T_HIT."""
    r = 0.5 * dp
    pts = []
    for start, d, v in shots:
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
             ("TIMESTEP", "MAXTIME"): "8.05E-06", ("TIMESTEP", "RUNTIME"): "8.05E-06",
             ("TIMESTEP", "WRITIME"): "4.0E-06", ("TIMESTEP", "PRNTIME"): "4.0E-06"},
        particles=wall_shots(),
        checks=["particles_moved", "wall_collisions"],
        same_as=None if n == 1 else "box/walls/np1",
        # The gas in this closed box differs slightly between decompositions
        # (~0.2% of a 4e-4 m/s flow driven by the particles' volume fraction
        # with FEEDBACK 0), so the gas seen by the particles (rprop) does too.
        # Particle state only (position y, velocity, temperature, spin)
        same_as_fields=("y", "vp_", "Tp", "omp_"))

# -- Gas only: conservation of the flow solver itself --------------------------
JUMP = {("PICL", "USED"): 0, ("INITFLOW", "RVAL1"): "0.0043", ("INITFLOW", "RVAL3"): "0.0",
        ("INITFLOW", "RVAL5"): "0.97", ("INITFLOW", "RVAL6"): "82700.0"}
add(name="box/gas-jump/closed/np2",
    description="Pressure/density jump in a slip-wall box across two regions; "
                "mass and energy must be conserved to round-off",
    nranks=2, grid=BOX, bc="closed", inp={**uniform_gas(0.0), **JUMP},
    checks=["mass", "gas_energy"])
add(name="box/gas-jump/closed-creconst0/np2",
    description="Same jump in the slip-wall box with CRECONSTC 0 (unconstrained "
                "cell gradients)",
    nranks=2, grid=BOX, bc="closed",
    inp={**uniform_gas(0.0), **JUMP, ("NUMERICS", "CRECONSTC"): 0},
    checks=["mass", "gas_energy"])
# Periodic virtual cells used to compute their own gradients, so the two sides
# of a periodic face disagreed (mass/energy drift ~4e-5 in 100 us; 11dc640)
for n in (1, 2):
    add(name="box/gas-jump/periodic/np%d" % n,
        description="Same jump in a fully periodic box; second-order fluxes across "
                    "periodic faces must conserve",
        nranks=n, grid=BOX, bc="periodic", inp={**uniform_gas(0.0), **JUMP},
        checks=["mass", "gas_energy"])
# Symmetry virtual cells had the same fault (drift ~3e-5 in 100 us; 11dc640)
for n in (1, 2):
    add(name="box/gas-jump/symmetry/np%d" % n,
        description="Same jump with symmetry patches on all six faces; second-order "
                    "fluxes across symmetry faces must conserve",
        nranks=n, grid=BOX, bc="symmetry", inp={**uniform_gas(0.0), **JUMP},
        checks=["mass", "gas_energy"])


# -- Periodic box: uniform lattice with the torque on -------------------------
# One particle per cell at its centre: phi and the forces are the same in all
# cells and stay so while the particles move together. The torque reads the
# gradient of the feedback force, which reads the virtual cells.
TORQUE_LATTICE = [(0.0005 + i * 1.0e-3, 0.0005 + j * 1.0e-3, 0.0005 + k * 1.0e-3, 7.26e-4)
                  for i in range(8) for j in range(8) for k in range(8)]
for n in (1, 2, 4):
    add(name="box/lattice-torque/np%d" % n,
        description="Uniform lattice with COLLISIONS 2, whose torque reads the gradient "
                    "of the feedback force; the gas must stay uniform from cell to cell",
        nranks=n, grid=BOX, bc="periodic",
        inp={**uniform_gas(50.0), **TWO_WAY,
             ("TIMESTEP", "MAXTIME"): "2.0E-06", ("TIMESTEP", "RUNTIME"): "2.0E-06",
             ("TIMESTEP", "WRITIME"): "1.0E-06", ("TIMESTEP", "PRNTIME"): "1.0E-06",
             ("PICL", "QSFORCE"): 3, ("PICL", "COLLISIONS"): 2,
             ("PICL", "NEIGHBORWIDTH"): "3.0e-3"},
        particles=TORQUE_LATTICE,
        checks=["uniform_gas_cells"])


# -- 2D: DIMENS 2, one cell thick in z, BC_VIRTUAL z patches -------------------
# Particles on the mid-plane. Each case mirrors a 3D one above in x and y; in_plane
# requires that z stays inactive.
TWO_D = {("NUMERICS", "DIMENS"): 2}
TWO_D_TUBE = Grid(lx=0.04, ly=0.004, lz=0.001, x0=-0.02, y0=0.0, z0=0.0, nx=40, ny=4, nz=1)
TWO_D_BOX = Grid(lx=0.008, ly=0.008, lz=0.001, x0=0.0, y0=0.0, z0=0.0, nx=8, ny=8, nz=1)
TWO_D_THIN = Grid(lx=0.008, ly=0.008, lz=0.0001, x0=0.0, y0=0.0, z0=0.0, nx=8, ny=8, nz=1)
TWO_D_ELONGATED = Grid(lx=0.016, ly=0.008, lz=0.001, x0=0.0, y0=0.0, z0=0.0, nx=8, ny=8, nz=1)
MID_PLANE = 0.0005
_3D = {c.name: c for c in CASES}


def on_plane(particles):
    return [p[:2] + (MID_PLANE,) + p[3:] for p in particles]


for n in (1, 2, 4):
    add(name="2d/shocktube/two-way/np%d" % n,
        description="Shock over four particles with feedback to the gas, in 2D",
        nranks=n, grid=TWO_D_TUBE, bc="shocktube-2d",
        inp={**SHOCKTUBE, **TWO_WAY, **TWO_D}, particles=on_plane(FOUR_PARTICLES),
        checks=["particles_moved", "in_plane"],
        same_as=None if n == 1 else "2d/shocktube/two-way/np1")
for n in (1, 2):
    add(name="2d/momentum-exchange/np%d" % n,
        description="Gas at 50 m/s through 40 particles at rest in a box periodic in x "
                    "and y; momentum and energy conserved in 2D",
        nranks=n, grid=TWO_D_BOX, bc="periodic-2d",
        inp={**uniform_gas(50.0), **TWO_WAY, **TWO_D}, particles=on_plane(cloud(40, 2.0e-5)),
        checks=["particles_moved", "mass", "momentum", "energy", "in_plane"],
        same_as=None if n == 1 else "2d/momentum-exchange/np1",
        restart_test=(n == 2))
    add(name="2d/linear-interp/np%d" % n,
        description="Gas density linear in x and y; particles at rest must interpolate "
                    "it exactly in 2D",
        nranks=n, grid=TWO_D_BOX, bc="closed-2d",
        inp={**uniform_gas(0.0), **ONE_WAY, **LINEAR_RHO, ("INITFLOW", "RVAL9"): "0.0",
             **TWO_D},
        particles=on_plane(cloud(40, 1.0e-6, lo=0.0006, hi=0.0074)),
        checks=["linear_interp", "mass", "in_plane"])
    add(name="2d/elongated-projection/np%d" % n,
        description="One particle at rest in 2:1 cells; its volume fraction must be "
                    "centred on it and spread ~0.44 cell lengths in x and y",
        nranks=n, grid=TWO_D_ELONGATED, bc="closed-2d",
        inp={**_3D["box/elongated-projection/np1"].inp, **TWO_D},
        particles=[(0.0066, 0.0044, MID_PLANE, 1.0e-4)],
        checks=["projection_shape", "mass", "in_plane"])

# Same pairs idea as box/ghost-collisions, in the plane: across the bin faces at
# x = 4 and y = 4 mm, the periodic faces and the periodic corner
TWO_D_GHOST_OFFSET = (85.0e-6, 35.0e-6)
TWO_D_GHOST_CENTRES = [(3.0, 5.0), (4.0, 1.5), (1.5, 4.0), (4.0, 4.0),
                       (0.0, 6.5), (5.5, 0.0), (0.0, 0.0)]


def collision_pairs_2d():
    pts = []
    for c in TWO_D_GHOST_CENTRES:
        for sign in (-0.5, 0.5):
            pts.append(tuple(round((c[l] * 1e-3 + sign * TWO_D_GHOST_OFFSET[l]) % BOX.lx, 9)
                             for l in range(2)) + (MID_PLANE, GHOST_DP))
    return pts


TWO_D_BINARY = on_plane(
    collision_pair((4.0, 1.5, 0.5), 0, 1.0e-4, 1.0e-4, 10.0)           # equal, x bin face
    + collision_pair((6.5, 4.0, 0.5), 1, 1.0e-4, 5.0e-5, 10.0)         # 8:1, y bin face
    + collision_pair((0.0, 6.5, 0.5), 0, 1.0e-4, 5.0e-5, 1.0)          # 8:1, periodic x
    + collision_pair((2.5, 5.5, 0.5), 0, 1.0e-4, 1.0e-4, 1.0, 0.5))    # oblique, spins about z
TWO_D_WALL_SHOTS = [
    ((None, 2.3, 0.5), (-1, 0, 0), 1.0),     # x- face, Hertz regime
    ((None, 6.3, 0.5), (-1, 0, 0), 10.0),    # x- face
    ((2.4, None, 0.5), (0, -1, 0), 100.0),   # y- face, centre crosses the wall
    ((None, 5.5, 0.5), (1, 0, 0), 300.0),    # x+ face on a quad diagonal
    ((None, None, 0.5), (1, -1, 0), 10.0),   # x+/y- box corner
    ((5.6, None, 0.5), (1, 1, 0), 10.0),     # y+ face at 45 degrees
]
for n in (1, 2, 4):
    add(name="2d/ghost-collisions/np%d" % n,
        description="Overlapping in-plane pairs pushed apart by collisions only, across "
                    "bin faces, rank borders and periodic faces",
        nranks=n, grid=TWO_D_BOX, bc="periodic-2d",
        inp={**uniform_gas(0.0), **ONE_WAY, **COLLISIONS_ONLY, **TWO_D},
        particles=collision_pairs_2d(),
        checks=["particles_moved", "collision_pairs", "in_plane"],
        same_as=None if n == 1 else "2d/ghost-collisions/np1")
    add(name="2d/binary-collisions/np%d" % n,
        description="Approaching in-plane pairs (equal and 8:1 mass, head-on and oblique); "
                    "momentum and angular momentum conserved, spin only about z",
        nranks=n, grid=TWO_D_BOX, bc="periodic-2d",
        inp={**_3D["box/binary-collisions/np1"].inp, **TWO_D}, particles=TWO_D_BINARY,
        checks=["particles_moved", "binary_collisions", "in_plane"],
        same_as=None if n == 1 else "2d/binary-collisions/np1")
    add(name="2d/walls/np%d" % n,
        description="Particles hitting the x and y slip walls (PARTWALL 1) in 2D: faces, "
                    "a quad diagonal and a corner; each rebounds and ends inside the box",
        nranks=n, grid=TWO_D_BOX, bc="partwall-2d",
        inp={**_3D["box/walls/np1"].inp, **TWO_D},
        particles=wall_shots(shots=TWO_D_WALL_SHOTS),
        checks=["particles_moved", "wall_collisions", "in_plane"])

# Gas only: the same jump as box/gas-jump, at two mesh thicknesses
add(name="2d/gas-jump/lz1mm/np1",
    description="Pressure/density jump in a 2D slip-wall box, 1 mm thick",
    nranks=1, grid=TWO_D_BOX, bc="closed-2d", inp={**uniform_gas(0.0), **JUMP, **TWO_D},
    checks=["mass", "gas_energy"])
add(name="2d/gas-jump/lz1mm/np2",
    description="Same on two regions; conservation across a region border in 2D",
    nranks=2, grid=TWO_D_BOX, bc="closed-2d", inp={**uniform_gas(0.0), **JUMP, **TWO_D},
    checks=["mass", "gas_energy"])
add(name="2d/gas-jump/lz100um/np1",
    description="Same jump 0.1 mm thick; the gas must not depend on the thickness",
    nranks=1, grid=TWO_D_THIN, bc="closed-2d", inp={**uniform_gas(0.0), **JUMP, **TWO_D},
    checks=["mass", "gas_energy"], same_gas_as="2d/gas-jump/lz1mm/np1")


# Cases that also write the old-style output (OLD_OUTPUT_WRT 1): the harness checks
# that it agrees with the dumps and with the post program (2 regions and 1, NSCBC
# outflow with history, random streams, gas only)
OLD_OUTPUT_CASES = (
    "box/momentum-exchange/np2",
    "shocktube/two-way/viscous-unsteady/np2",
    "shocktube/one-way/1-particle/np1",
    "box/qs-fluctuations/np2",
    "box/gas-jump/periodic/np2",
    "2d/momentum-exchange/np2",
)
for _case in CASES:
    if _case.name in OLD_OUTPUT_CASES:
        _case.inp[("POST", "OLD_OUTPUT_WRT")] = 1


def by_name():
    return {c.name: c for c in CASES}


# -- Periodicity: Newton's cradles and sectors of a cylinder -------------------
# Sector meshes come from utilities/mesh/sector_mesh.py (bc "sector": its .bc,
# the two planes a rotational periodic pair). Gas at rest, no gas forces:
# particles move by collisions only, or not at all, so the reference is exact.
from .case import SectorGrid  # noqa: E402

COLLISIONS_ELASTIC = {("PICL", "QSFORCE"): 0, ("PICL", "COLLISIONS"): 1,
                      ("PICL", "KSP"): 10, ("PICL", "EREST"): "1.0",
                      ("TIMESTEP", "MAXTIME"): "3.0E-05", ("TIMESTEP", "RUNTIME"): "3.0E-05",
                      ("TIMESTEP", "WRITIME"): "1.5E-05", ("TIMESTEP", "PRNTIME"): "1.5E-05"}
NO_FORCES = {("PICL", "QSFORCE"): 0, ("PICL", "AMFORCE"): 0, ("PICL", "PGFORCE"): 0,
             ("PICL", "COLLISIONS"): 0,
             ("TIMESTEP", "MAXTIME"): "2.0E-05", ("TIMESTEP", "RUNTIME"): "2.0E-05",
             ("TIMESTEP", "WRITIME"): "1.0E-05", ("TIMESTEP", "PRNTIME"): "1.0E-05"}


def cradle(axis_pos, normal, along, dp, n=5, gap=1.0e-3, speed=10.0, z=None):
    """A resting chain of n touching particles centred on axis_pos along
    `along` (unit vector), and one particle hitting it from the + side at
    `speed`. Returns (x, y, z, dp, u, v, w) tuples."""
    pitch = dp * (1.0 + gap)
    out = []
    for k in range(n + 1):
        s = (k - (n - 1) / 2.0) * pitch + 0.37 * gap * dp
        if k == n:
            s += 0.5 * dp
        p = [axis_pos[i] + s * along[i] for i in range(3)]
        if z is not None:
            p[2] = z
        v = [-speed * along[i] if k == n else 0.0 for i in range(3)]
        out.append(tuple(p) + (dp,) + tuple(v))
    return out


def _rot_z_point(p, deg):
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return (c * p[0] - s * p[1], s * p[0] + c * p[1]) + tuple(p[2:3])


def rotated(particles, deg):
    out = []
    for p in particles:
        x = _rot_z_point(p[:3], deg)
        v = _rot_z_point(p[4:7], deg)
        out.append(x + (p[3],) + v)
    return out


def wrap_quarter(particles):
    """Particles below y = 0 moved to their periodic image near x = 0."""
    return [rotated([p], 90.0)[0] if p[1] < 0.0 else p for p in particles]


QUARTER = SectorGrid(angle=90.0, r_out=0.02, ntheta=48, core=(0.002, 0.004))
_chain = cradle((0.015, 0.0, 0.0), None, (0.0, 1.0, 0.0), 1.0e-4, z=5.0e-5)
add(name="sector/cradle/np4",
    description="Newton's cradle across a rotational periodic plane (y = 0) of a quarter "
                "cylinder, and its twin rotated 45 degrees into the sector: the momentum "
                "must cross the plane (rotated ghosts) exactly as in the twin",
    nranks=4, grid=QUARTER, bc="sector",
    inp={**uniform_gas(0.0), **ONE_WAY, **TWO_D, **COLLISIONS_ELASTIC},
    particles=wrap_quarter(_chain) + rotated(_chain, 45.0),
    checks=["particles_moved", "cradle_twin", "in_plane"],
    twin=("rotate", 45.0, 90.0))

_lchain = cradle((0.0, 0.00237, 0.00263), None, (1.0, 0.0, 0.0), 1.0e-4)
_lchain_in = [((p[0] + 0.008) if p[0] < 0.0 else p[0],) + p[1:] for p in _lchain]
_lshift = (0.004, 0.003, 0.0)
add(name="box/cradle-linear/np2",
    description="Newton's cradle across the periodic face x = 0 of a box, and its twin "
                "shifted into the box: the momentum must cross the face exactly as in "
                "the twin",
    nranks=2, grid=BOX, bc="periodic",
    inp={**uniform_gas(0.0), **ONE_WAY, **COLLISIONS_ELASTIC},
    particles=_lchain_in + [(p[0] + _lshift[0], p[1] + _lshift[1], p[2] + _lshift[2]) + p[3:]
                            for p in _lchain],
    checks=["particles_moved", "cradle_twin"],
    twin=("shift", _lshift, (0.008, 0.008, 0.008)))


def _polar_particle(r, deg, v_r, v_t, dp, z):
    t = math.radians(deg)
    c, s = math.cos(t), math.sin(t)
    return (r * c, r * s, z, dp, v_r * c - v_t * s, v_r * s + v_t * c, 0.0)


WEDGE5 = SectorGrid(angle=5.0, theta0=-2.5, r_in=0.02, r_out=0.03, ntheta=12)
add(name="2d/sector-wedge-crossing/np2",
    description="Particles without forces crossing the planes of a 5 degree wedge "
                "several times (500 m/s across it): straight paths brought back by the "
                "rotation, to round-off",
    nranks=2, grid=WEDGE5, bc="sector",
    inp={**uniform_gas(0.0), **ONE_WAY, **TWO_D, **NO_FORCES},
    particles=[_polar_particle(0.025, 0.3, 40.0, 500.0, 2.0e-5, 5.0e-5),
               _polar_particle(0.0252, -1.1, -30.0, -500.0, 2.0e-5, 5.0e-5),
               _polar_particle(0.0248, 1.9, 0.0, 120.0, 2.0e-5, 5.0e-5)],
    checks=["particles_moved", "ballistic_periodic", "in_plane"])

QUARTER3D = SectorGrid(angle=90.0, r_out=0.01, ntheta=16, core=(0.002, 0.004),
                       z=(0.0, 0.001), nz=4)
add(name="sector/z-periodic-crossing/np4",
    description="3D quarter cylinder, z periodic: particles without forces crossing the "
                "rotational planes, the z faces, and both planes near the axis",
    nranks=4, grid=QUARTER3D, bc="sector",
    inp={**uniform_gas(0.0), **ONE_WAY, **NO_FORCES},
    particles=[(0.006, 0.0007, 0.00041, 2.0e-5, 20.0, -120.0, 37.0),
               (0.0005, 0.0006, 0.00013, 2.0e-5, -55.0, -61.0, -29.0),
               (0.0031, 0.0042, 0.00087, 2.0e-5, 13.0, 9.0, 61.0)],
    checks=["particles_moved", "ballistic_periodic"])

# All particles near one z face, not spanning z: the bin box in z must cover
# the whole period (BinMinLen 0.25 mm here, the period 1 mm)
QUARTER3D_THIN = SectorGrid(angle=90.0, r_out=0.01, ntheta=16, core=(0.002, 0.004),
                            z=(0.0, 0.001), nz=8)
add(name="sector/z-cluster-near-face/np2",
    description="3D quarter cylinder, z periodic: particles without forces clustered "
                "near the face z = 1 mm, crossing it to the face z = 0",
    nranks=2, grid=QUARTER3D_THIN, bc="sector",
    inp={**uniform_gas(0.0), **ONE_WAY, **NO_FORCES},
    particles=[_polar_particle(0.006, deg, 5.0, 0.0, 2.0e-5, z)[:6] + (20.0,)
               for deg, z in ((40.0, 0.00090), (45.0, 0.00092), (50.0, 0.00095))],
    checks=["particles_moved", "ballistic_periodic"])

# Quarter cylinder against the full cylinder of four exact copies: radial
# blast (case cylds: r < RVAL1 at the post-shock state), particles frozen,
# one-way coupling; the gas at the particles and their volume fraction must
# match the unrotated copy in the full cylinder
QUARTER_SMALL = SectorGrid(angle=90.0, r_out=0.012, ntheta=32, core=(0.002, 0.004))
FULL_SMALL = SectorGrid(angle=90.0, r_out=0.012, ntheta=32, core=(0.002, 0.004), full=True)
BLAST = {("INITFLOW", "RVAL1"): "0.006", ("INITFLOW", "RVAL3"): "0.0",
         ("PICL", "VELOCITYMODEL"): 0, ("PICL", "PGFORCE"): 1,
         ("TIMESTEP", "MAXTIME"): "4.0E-06", ("TIMESTEP", "RUNTIME"): "4.0E-06",
         ("TIMESTEP", "WRITIME"): "2.0E-06", ("TIMESTEP", "PRNTIME"): "2.0E-06"}


def _quarter_cloud(n, dp, seed=31):
    state = seed
    out = []
    while len(out) < n:
        vals = []
        for _ in range(2):
            state = (1103515245 * state + 12345) % 2 ** 31
            vals.append(state / 2 ** 31)
        r, t = 0.0045 + 0.0045 * vals[0], 90.0 * vals[1]
        out.append(_polar_particle(r, t, 0.0, 0.0, dp, 5.0e-5)[:4])
    return out


_qcloud = _quarter_cloud(24, 5.0e-5)
_fcloud = [q for p in _qcloud for q in
           [p] + [rotated([p + (0.0, 0.0, 0.0)], 90.0 * k)[0][:4] for k in (1, 2, 3)]]
add(name="sector/quarter-vs-full/full/np8",
    description="Full cylinder of four copies of the quarter mesh, radial blast, frozen "
                "particles in four rotated copies (reference of the quarter case)",
    nranks=8, grid=FULL_SMALL, bc="sector", casename="cylds",
    inp={**ONE_WAY, **TWO_D, **BLAST}, particles=_fcloud,
    checks=["in_plane"])
add(name="sector/quarter-vs-full/quarter/np2",
    description="Quarter cylinder with rotational periodic planes against the full "
                "cylinder: gas interpolated at the particles and their volume fraction "
                "(cell images across the planes and around the axis)",
    nranks=2, grid=QUARTER_SMALL, bc="sector", casename="cylds",
    inp={**ONE_WAY, **TWO_D, **BLAST}, particles=_qcloud,
    checks=["in_plane"], quarter_of="sector/quarter-vs-full/full/np8")

# Stencils larger than the old fixed 125 cells: near 45 degrees and the
# refinement ring the cells' bounding boxes are large, and all 4 z layers
# are in the support (z = 0.29-0.30 mm): 152, 144 and 144 cells
QUARTER3D_RING = SectorGrid(angle=90.0, r_out=0.012, ntheta=16, core=(0.002, 0.004),
                            max_cell=4.0e-4, z=(0.0, 0.001), nz=4)
add(name="sector/stencil-bound/np4",
    description="3D quarter cylinder with a refinement ring, radial blast, frozen "
                "particles whose interpolation stencils hold more than 125 cells",
    nranks=4, grid=QUARTER3D_RING, bc="sector", casename="cylds",
    inp={**ONE_WAY, **BLAST},
    particles=[(0.003719, 0.003719, 0.00030, 5.0e-5),
               (0.003726, 0.003598, 0.00030, 5.0e-5),
               (0.003598, 0.003726, 0.00029, 5.0e-5)],
    checks=[])


# -- Particles outside the fluid domain: the run must stop -------------------
def _outside_stop(particle, patch, t=None):
    lines = [r"^PPICLF: particle %d \(line of points\.dat\) is outside the fluid domain$"
             % particle,
             r"^  beyond patch %d by " % patch,
             r"PPICLF: ERROR .* Particles outside the fluid domain:"]
    if t is not None:
        lines.append(r"^  t = +%s s, x y z = " % re.escape(t))
    return lines


STOP_TIME = {("TIMESTEP", "MAXTIME"): "3.0E-06", ("TIMESTEP", "RUNTIME"): "3.0E-06",
             ("TIMESTEP", "WRITIME"): "1.5E-06", ("TIMESTEP", "PRNTIME"): "1.5E-06"}
add(name="stop/outside-at-start/np2",
    description="A particle outside the outer arc of a quarter cylinder (inside its "
                "bounding box) at t = 0: the run stops at start-up naming it",
    nranks=2, grid=QUARTER, bc="sector",
    inp={**uniform_gas(0.0), **ONE_WAY, **TWO_D, **NO_FORCES, **STOP_TIME},
    particles=[_polar_particle(0.015, 45.0, 0.0, 0.0, 1.0e-4, 5.0e-5),
               _polar_particle(0.0201, 45.0, 0.0, 0.0, 1.0e-4, 5.0e-5)],
    expect_stop=_outside_stop(2, 3, t="0.000000E+00"))
add(name="stop/outflow-near-plane/np4",
    description="A particle leaving a quarter cylinder through the outflow 0.3 degrees "
                "from a periodic plane (rotated cell images): the run stops",
    nranks=4, grid=QUARTER, bc="sector",
    inp={**uniform_gas(0.0), **ONE_WAY, **TWO_D, **NO_FORCES, **STOP_TIME},
    particles=[_polar_particle(0.015, 45.0, 0.0, 0.0, 1.0e-4, 5.0e-5),
               _polar_particle(0.0195, 0.3, 400.0, 0.0, 1.0e-4, 5.0e-5)],
    expect_stop=_outside_stop(2, 3))
add(name="stop/wall-without-partwall/np1",
    description="A particle reaching a slip wall of a closed box without PARTWALL (no "
                "particle wall): the run stops when it crosses the wall",
    nranks=1, grid=BOX, bc="closed",
    inp={**uniform_gas(0.0), **ONE_WAY, **NO_FORCES, **STOP_TIME},
    particles=[(0.004, 0.004, 0.004, 1.0e-4),
               (0.004, 0.004, 0.0005, 1.0e-4, 0.0, 0.0, -400.0)],
    expect_stop=_outside_stop(2, 5))
