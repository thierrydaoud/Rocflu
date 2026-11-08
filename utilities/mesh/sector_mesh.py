#!/usr/bin/env python3
"""Hex meshes of a circular sector about the z axis, for Rocflu (CENTAUR ASCII
<case>.hyb.asc), with a matching <case>.bc template.

The sector spans theta0 to theta0 + angle (degrees, at most 180), or the whole
cylinder (--full: four exact 90 degree copies of a 90 degree sector at
theta0 = 0). All cells are hexahedra: the 2D quads are extruded in z with
--nz layers (one layer: a DIMENS 2 mesh).

Layouts:
  annulus  (--r-in R): polar cells from r_in > 0 to r_out; a wall at r_in.
           Any angle; the only choice for thin wedges.
  core     (--core A R_POLAR): a quadrilateral core at the axis (a square for
           90 degrees: the origin, A on each plane, and the point where the
           normals to the planes at A meet), a transition ring from its outer
           edges to the circle r = R_POLAR, then polar cells to r_out. The
           core has ntheta/2 cells per side. Refused when a cell angle falls
           below 30 degrees (thin sectors: use the annulus).

Polar cells (--grading):
  square   r grows geometrically so cells are about square (dr = r dtheta).
           With --max-cell L, a refinement ring triples the cells around
           (each cell's outer edge split in three, the cell into 4 quads)
           wherever the cells would grow past L, so they stay between about
           L/3 and L. Default.
  uniform  --nr rings of equal width, no refinement rings.

Patches, in this order: plane1 (theta0), plane2 (theta0 + angle), outer,
inner (annulus only), zlow, zhigh; with --full: outer, inner, zlow, zhigh.

The .bc template: the two planes as a rotational periodic pair (ANGLE takes
the related patch onto the section's patch: plane1 -angle, plane2 +angle),
outer as a non-reflecting outflow, inner as a slip wall, z patches as
BC_VIRTUAL for one layer, else a linear periodic pair (--z-bc periodic,
default) or slip walls (--z-bc slip). Edit it for your case.

--dp D reports the smallest cell from each radius outward, in particle
diameters (cells should be a few dp wherever particles go).

Examples:
  # 2D quarter cylinder to 0.1 m, square core 1.9 mm, polar from 3.8 mm,
  # 48 cells around, refinement rings keeping cells below 0.6 mm
  sector_mesh.py quarter --angle 90 --r-out 0.1 --ntheta 48 \\
      --core 0.0019 0.0038 --max-cell 6e-4 --z 0 1e-4 --dp 1e-4 --preview
  # the same in 3D, 10 layers over 1 mm, z periodic
  sector_mesh.py quarter3d --angle 90 --r-out 0.1 --ntheta 48 \\
      --core 0.0019 0.0038 --max-cell 6e-4 --z 0 1e-3 --nz 10
  # 5 degree wedge about x, from r = 0.02 m
  sector_mesh.py wedge --angle 5 --theta0 -2.5 --r-in 0.02 --r-out 0.1 \\
      --ntheta 24 --z 0 1e-4
"""
import argparse
import math
import sys

import numpy as np

MIN_ANGLE = 30.0


# ---------------------------------------------------------------------------
# 2D quads: arrays (4, 2), counter-clockwise or not (oriented later)
# ---------------------------------------------------------------------------

def _cos_sin(th):
    # exact zeros at multiples of 90 degrees, so rotated copies match exactly
    th = np.atleast_1d(np.asarray(th, float))
    c, s = np.cos(th), np.sin(th)
    c[np.abs(c) < 1e-15] = 0.0
    s[np.abs(s) < 1e-15] = 0.0
    return c, s


def _pt(r, th):
    c, s = _cos_sin(th)
    return np.array([r * c[0], r * s[0]])


def _ring(r0, r1, theta0, angle, n):
    th = theta0 + angle * np.arange(n + 1) / n
    c, s = _cos_sin(th)
    p0 = np.column_stack([r0 * c, r0 * s])
    p1 = np.column_stack([r1 * c, r1 * s])
    return [np.array([p0[i], p0[i + 1], p1[i + 1], p1[i]]) for i in range(n)]


def _tripling_ring(r0, r1, theta0, angle, n):
    """One ring from n cells (radius r0) to 3n cells (radius r1)."""
    # angles computed as in _ring, so the vertices shared with the rings
    # inside (n cells) and outside (3n cells) are bit for bit the same
    quads = []
    rm = r0 + (r1 - r0) / 3.0
    c, s = _cos_sin(theta0 + angle * np.arange(n + 1) / n)
    inner = np.column_stack([r0 * c, r0 * s])
    c, s = _cos_sin(theta0 + angle * np.arange(3 * n + 1) / (3 * n))
    outer = np.column_stack([r1 * c, r1 * s])
    mid = np.column_stack([rm * c, rm * s])
    for i in range(n):
        A, B = inner[i], inner[i + 1]
        D, M1, M2, C = outer[3 * i], outer[3 * i + 1], outer[3 * i + 2], outer[3 * i + 3]
        P1, P2 = mid[3 * i + 1], mid[3 * i + 2]
        quads += [np.array([A, B, P2, P1]), np.array([A, P1, M1, D]),
                  np.array([P1, P2, M2, M1]), np.array([P2, B, C, M2])]
    return quads


def polar(theta0, angle, r0, r_out, ntheta, grading, nr, max_cell):
    """Polar quads from r0 to r_out; returns (quads, list of (radius, ntheta)
    where each refinement ring ends)."""
    if grading == "uniform":
        if not nr:
            sys.exit("--grading uniform needs --nr")
        r = np.linspace(r0, r_out, nr + 1)
        return [q for k in range(nr) for q in _ring(r[k], r[k + 1], theta0, angle, ntheta)], []
    quads, rings, r, n = [], [], r0, ntheta
    while r < r_out * (1 - 1e-12):
        d = angle / n
        step = r * (math.exp(d) - 1.0)          # about square
        nxt = r + step
        if max_cell and r * d >= max_cell * (1 - 1e-9) and r_out - nxt > 3 * step:
            quads += _tripling_ring(r, nxt, theta0, angle, n)
            n *= 3
            rings.append((nxt, n))
        else:
            if nxt > r_out or r_out - nxt < 0.5 * step:
                nxt = r_out
            quads += _ring(r, nxt, theta0, angle, n)
        r = nxt
    return quads, rings


def core(theta0, angle, a, r_polar, ntheta, ntrans):
    """Quadrilateral core at the axis and the transition ring to r_polar."""
    if ntheta % 2:
        sys.exit("--core needs an even --ntheta")
    n = ntheta // 2
    P1, P2 = _pt(a, theta0), _pt(a, theta0 + angle)
    M = _pt(a / math.cos(angle / 2), theta0 + angle / 2)
    s = np.linspace(0.0, 1.0, n + 1)
    S, T = np.meshgrid(s, s, indexing="ij")
    blk = ((S * (1 - T))[..., None] * P1 + (S * T)[..., None] * M
           + ((1 - S) * T)[..., None] * P2)
    quads = [np.array([blk[i, j], blk[i + 1, j], blk[i + 1, j + 1], blk[i, j + 1]])
             for i in range(n) for j in range(n)]
    edge = np.concatenate([blk[n, :, :], blk[::-1, n, :][1:]], axis=0)
    th = theta0 + angle * np.arange(ntheta + 1) / ntheta
    c, s_ = _cos_sin(th)
    arc = np.column_stack([r_polar * c, r_polar * s_])
    if not ntrans:
        ntrans = max(1, int(round((r_polar - a) / (r_polar * angle / ntheta))))
    w = np.linspace(0.0, 1.0, ntrans + 1)
    tr = (1 - w)[None, :, None] * edge[:, None, :] + w[None, :, None] * arc[:, None, :]
    quads += [np.array([tr[i, j], tr[i + 1, j], tr[i + 1, j + 1], tr[i, j + 1]])
              for i in range(ntheta) for j in range(ntrans)]
    return quads


def quality(quads):
    """Smallest angle (deg), largest aspect ratio, smallest and largest edge."""
    q = np.array(quads)
    e = np.roll(q, -1, axis=1) - q
    le = np.linalg.norm(e, axis=2)
    u, v = e, -np.roll(e, 1, axis=1)
    cos = (u * v).sum(-1) / (le * np.roll(le, 1, axis=1))
    ang = np.degrees(np.arccos(np.clip(cos, -1, 1)))
    asp = (le.max(1) / le.min(1)).max()
    return ang.min(), asp, le.min(), le.max()


def rotate90(q, k):
    x, y = q[..., 0], q[..., 1]
    for _ in range(k % 4):
        x, y = -y, x
    return np.stack([x, y], axis=-1)


# ---------------------------------------------------------------------------
# 3D mesh and files
# ---------------------------------------------------------------------------

def build(quads, z0, z1, nz, theta0, angle, r_out, r_in, full):
    index, xy2, cells = {}, [], []
    for q in quads:
        ids = []
        for p in q:
            key = (round(p[0], 15), round(p[1], 15))
            if key not in index:
                index[key] = len(xy2)
                xy2.append(p)
            ids.append(index[key])
        cells.append(ids)
    xy2 = np.array(xy2)
    n2 = len(xy2)
    z = np.linspace(z0, z1, nz + 1)
    xyz = np.column_stack([np.tile(xy2, (nz + 1, 1)), np.repeat(z, n2)])

    # 2D boundary edges (used by one cell) and their patch
    count = {}
    for ids in cells:
        for m in range(4):
            e = tuple(sorted((ids[m], ids[(m + 1) % 4])))
            count[e] = count.get(e, 0) + 1
    tol = 1e-9 * r_out
    n1 = _pt(1.0, theta0 + math.pi / 2)            # normal of plane 1
    n2v = _pt(1.0, theta0 + angle + math.pi / 2)   # normal of plane 2
    edges = {"plane1": [], "plane2": [], "outer": [], "inner": []}
    for e, c in count.items():
        if c != 1:
            continue
        a, b = xy2[e[0]], xy2[e[1]]
        ra, rb = np.hypot(*a), np.hypot(*b)
        if abs(ra - r_out) < tol and abs(rb - r_out) < tol:
            edges["outer"].append(e)
        elif r_in is not None and abs(ra - r_in) < tol and abs(rb - r_in) < tol:
            edges["inner"].append(e)
        elif not full and abs(a @ n1) < tol and abs(b @ n1) < tol:
            edges["plane1"].append(e)
        elif not full and abs(a @ n2v) < tol and abs(b @ n2v) < tol:
            edges["plane2"].append(e)
        else:
            sys.exit("boundary edge not on a patch: %s %s" % (a, b))

    v3 = lambda v, k: k * n2 + v
    hexes = []
    for ids in cells:
        p = xy2[ids]
        area2 = sum(p[m, 0] * p[(m + 1) % 4, 1] - p[(m + 1) % 4, 0] * p[m, 1] for m in range(4))
        if area2 < 0:                 # bottom face counter-clockwise from +z
            ids = ids[::-1]
        for k in range(nz):
            hexes.append([v3(v, k) for v in ids] + [v3(v, k + 1) for v in ids])
    hexes = np.array(hexes)
    centroid = xyz[hexes].mean(axis=1)
    v2h = {}
    for ih, hv in enumerate(hexes):
        for v in hv:
            v2h.setdefault(v, []).append(ih)

    def outward(q):
        ih = set(v2h[q[0]]).intersection(*[v2h[v] for v in q[1:]]).pop()
        a, b, c, d = xyz[q]
        return q if np.dot(np.cross(c - a, d - b), xyz[q].mean(0) - centroid[ih]) > 0 else q[::-1]

    patches = []
    for name in ("plane1", "plane2", "outer", "inner"):
        if not edges[name]:
            continue
        # along the patch: by radius on planes, by angle on arcs
        key = (lambda e: np.hypot(*xy2[e[0]]) + np.hypot(*xy2[e[1]])) if name.startswith("plane") \
            else (lambda e: math.atan2(*xy2[e[0]][::-1]) + math.atan2(*xy2[e[1]][::-1]))
        quads3 = []
        for va, vb in sorted(edges[name], key=key):
            for k in range(nz):
                quads3.append(outward([v3(va, k), v3(vb, k), v3(vb, k + 1), v3(va, k + 1)]))
        patches.append((name, quads3))
    nc = len(cells)
    patches += [("zlow", [outward(list(hexes[i * nz][:4])) for i in range(nc)]),
                ("zhigh", [outward(list(hexes[i * nz + nz - 1][4:])) for i in range(nc)])]
    return xyz, hexes, patches


def write_hyb(path, title, xyz, hexes, patches):
    with open(path, "w") as f:
        w = f.write

        def ints(vals):
            vals = list(vals)
            for s in range(0, len(vals), 10):
                w("".join("%16d" % v for v in vals[s:s + 10]) + "\n")

        def reals(vals):
            vals = list(vals)
            for s in range(0, len(vals), 5):
                w("".join("%16.9E" % v for v in vals[s:s + 5]) + "\n")

        w("%-80s\n" % title[:80])
        w("%16d\n" % len(xyz))
        for c in range(3):
            reals(xyz[:, c])
        ints([0] * len(xyz))
        w("%16d\n" % 0)                       # tetrahedra
        w("%16d\n" % len(hexes))
        for c in range(8):
            ints(hexes[:, c] + 1)
        ints([0] * len(hexes))
        w("%16d\n" % 0)                       # prisms
        w("%16d\n" % 0)                       # pyramids
        w("%16d\n" % len(patches))
        ints([400] * len(patches))
        ints([0] * len(patches))
        ints(np.cumsum([len(q) for _, q in patches]))
        for name, _ in patches:
            w("%-80s\n" % name)
        w("%16d\n" % 0)                       # boundary triangles
        allq = np.array([q for _, qs in patches for q in qs])
        w("%16d\n" % len(allq))
        for c in range(4):
            ints(allq[:, c] + 1)


def write_bc(path, names, angle, nz, z_bc):
    num = {n: i + 1 for i, n in enumerate(names)}
    out = []
    if "plane1" in num:
        for this, other, sign in (("plane1", "plane2", -1), ("plane2", "plane1", 1)):
            out.append("# BC_PERIODIC\nPATCH     %d %d\nRELPATCH  %d %d\nANGLE     %.10g\n"
                       "AXIS      3\nNAME      %s\n#\n" % (num[this], num[this], num[other],
                                                         num[other], sign * angle, this))
    out.append("# BC_OUTFLOW\nPATCH      %d %d\nKIND       1\nNAME       Outer\nNSCBCK     0.0\n"
               "ORDER      1\nREFLECT    0\nTHRUSTFLAG 0\nTYPE       0\n#\n" % (num["outer"], num["outer"]))
    if "inner" in num:
        out.append("# BC_SLIPW\nPATCH      %d %d\nNAME       Inner\nPARTWALL   1\n#\n"
                   % (num["inner"], num["inner"]))
    lo, hi = num["zlow"], num["zhigh"]
    if nz == 1:
        out.append("# BC_VIRTUAL\nPATCH      %d %d\nNAME       Z\n#\n" % (lo, hi))
    elif z_bc == "periodic":
        for this, other in ((lo, hi), (hi, lo)):
            out.append("# BC_PERIODIC\nPATCH     %d %d\nRELPATCH  %d %d\nANGLE     0.0\n"
                       "AXIS      3\nNAME      Z-Periodic\n#\n" % (this, this, other, other))
    else:
        out.append("# BC_SLIPW\nPATCH      %d %d\nNAME       Z\nPARTWALL   1\n#\n" % (lo, hi))
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n# END\n")


def size_report(quads, dp):
    """Smallest cell edge from each radius outward."""
    q = np.array(quads)
    r = np.hypot(q[..., 0], q[..., 1]).min(1)
    le = np.linalg.norm(np.roll(q, -1, axis=1) - q, axis=2).min(1)
    order = np.argsort(r)
    tail_min = np.minimum.accumulate(le[order][::-1])[::-1]
    lines = []
    for rr in np.unique(np.round(np.geomspace(max(r[r > 0].min(), 1e-6), r.max(), 8), 6)):
        k = min(np.searchsorted(r[order], rr), len(order) - 1)
        lines.append("  from r = %8.4f m outward: smallest cell edge %.3g m (%.1f dp)"
                     % (rr, tail_min[k], tail_min[k] / dp))
    return "\n".join(lines)


def preview(path, quads, rings, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.add_collection(PolyCollection(quads, facecolors="none", edgecolors="#55a868", linewidths=0.3))
    ax.autoscale()
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def make(a):
    t0, ang = math.radians(a.theta0), math.radians(a.angle)
    if a.r_in is not None:
        quads, rings = polar(t0, ang, a.r_in, a.r_out, a.ntheta, a.grading, a.nr, a.max_cell)
    else:
        quads = core(t0, ang, a.core[0], a.core[1], a.ntheta, a.ntrans)
        if quality(quads)[0] < MIN_ANGLE:
            sys.exit("core: smallest cell angle %.1f deg < %.0f; use --r-in for this sector"
                     % (quality(quads)[0], MIN_ANGLE))
        pq, rings = polar(t0, ang, a.core[1], a.r_out, a.ntheta, a.grading, a.nr, a.max_cell)
        quads += pq
    if a.full:
        quads = [rotate90(q, k) for k in range(4) for q in quads]
    return quads, rings


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__[__doc__.index("Examples:"):])
    p.add_argument("case", help="writes <case>.hyb.asc and <case>.bc")
    p.add_argument("--angle", type=float, required=True, help="sector angle, degrees (<= 180)")
    p.add_argument("--theta0", type=float, default=0.0, help="angle of plane 1, degrees")
    p.add_argument("--r-out", type=float, required=True)
    p.add_argument("--ntheta", type=int, required=True, help="cells across the sector at the start")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--r-in", type=float, help="annulus from this radius (wall)")
    g.add_argument("--core", type=float, nargs=2, metavar=("A", "R_POLAR"),
                   help="core reaching A on each plane, transition to R_POLAR")
    p.add_argument("--ntrans", type=int, default=0, help="core transition rings (0: about square)")
    p.add_argument("--grading", choices=("square", "uniform"), default="square")
    p.add_argument("--nr", type=int, default=0, help="polar rings (uniform grading)")
    p.add_argument("--max-cell", type=float, default=0.0,
                   help="refinement rings keep polar cells below this size (square grading)")
    p.add_argument("--z", type=float, nargs=2, default=(0.0, 1.0e-4), metavar=("Z0", "Z1"))
    p.add_argument("--nz", type=int, default=1, help="cell layers in z (1: DIMENS 2)")
    p.add_argument("--z-bc", choices=("periodic", "slip"), default="periodic",
                   help="z patches in the .bc template when nz > 1")
    p.add_argument("--full", action="store_true", help="whole cylinder (90 degree sector, theta0 0)")
    p.add_argument("--dp", type=float, default=0.0, help="particle diameter, for the size report")
    p.add_argument("--preview", action="store_true", help="also write <case>_mesh.png")
    a = p.parse_args(argv)
    if a.angle > 180.0:
        sys.exit("sectors wider than 180 degrees are not supported")
    if a.full and (a.angle != 90.0 or a.theta0 != 0.0):
        sys.exit("--full needs --angle 90 --theta0 0")

    quads, rings = make(a)
    xyz, hexes, patches = build(quads, a.z[0], a.z[1], a.nz, math.radians(a.theta0),
                                math.radians(a.angle), a.r_out, a.r_in, a.full)
    names = [n for n, _ in patches]
    title = "sector %.6g deg%s, sector_mesh.py" % (a.angle, " (full)" if a.full else "")
    write_hyb(a.case + ".hyb.asc", title, xyz, hexes, patches)
    write_bc(a.case + ".bc", names, a.angle, a.nz, a.z_bc)
    amin, asp, emin, emax = quality(quads)
    print("%s.hyb.asc: %d vertices, %d hexes (%d layers)" % (a.case, len(xyz), len(hexes), a.nz))
    print("  patches: " + ", ".join("%d %s (%d)" % (i + 1, n, len(q)) for i, (n, q) in enumerate(patches)))
    print("  smallest cell angle %.1f deg, largest aspect ratio %.1f, in-plane edges %.3g to %.3g m"
          % (amin, asp, emin, emax))
    for r, n in rings:
        print("  refinement ring ending at r = %.4f m: %d cells around after it" % (r, n))
    if a.dp:
        print(size_report(quads, a.dp))
    print("%s.bc: template written (check it for your case)" % a.case)
    if a.preview:
        preview(a.case + "_mesh.png", quads, rings, title)
        print("%s_mesh.png: preview" % a.case)


if __name__ == "__main__":
    main()
