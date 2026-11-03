#!/usr/bin/env python3
"""Build a one-cell-thick DIMENS 2 wedge case from the 3D one in ../FluidParticles.

Usage: make_wedge_2d.py <dz [m]> <seed> <output directory> [maxtime [s]]

Mesh: the first layer of wedge.hyb.asc at thickness dz, same in-plane cells.
Bed: same annulus and volume fraction as the 3D bed, so N scales with dz. Random
placement without in-plane overlaps, on the mid-plane, at rest.
"""
import math
import os
import random
import re
import sys

SOURCE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "FluidParticles")
N_3D, LZ_3D = 40781, 1.0e-3
DP = 40.0e-6
R_IN, R_OUT, HALF_ANGLE = 0.250, 0.255, math.radians(0.5)
FIRST_LAYER_TOP = 2.0e-4


def rows(n, per_line):
    return (n + per_line - 1) // per_line


def read_block(f, n, per_line, convert):
    values = []
    for _ in range(rows(n, per_line)):
        line = f.readline().rstrip("\n")
        values += [convert(line[i:i + 16]) for i in range(0, len(line), 16)
                   if line[i:i + 16].strip()]
    if len(values) != n:
        sys.exit("grid read: expected %d values, got %d" % (n, len(values)))
    return values


def write_block(f, values, per_line, form):
    for i in range(0, len(values), per_line):
        f.write("".join(form % v for v in values[i:i + per_line]) + "\n")


def make_mesh(dz, path):
    with open(os.path.join(SOURCE, "wedge.hyb.asc")) as f:
        title = f.readline()
        nvert = int(f.readline())
        x = read_block(f, nvert, 5, float)
        y = read_block(f, nvert, 5, float)
        z = read_block(f, nvert, 5, float)
        read_block(f, nvert, 10, int)
        if int(f.readline()) != 0:
            sys.exit("expected no tetrahedra")
        nhex = int(f.readline())
        hexes = [read_block(f, nhex, 10, int) for _ in range(8)]
        read_block(f, nhex, 10, int)
        for kind in ("prisms", "pyramids"):
            if int(f.readline()) != 0:
                sys.exit("expected no " + kind)
        npatch = int(f.readline())
        kinds = read_block(f, npatch, 10, int)
        read_block(f, npatch, 10, int)
        last_quad = read_block(f, npatch, 10, int)
        names = [f.readline() for _ in range(npatch)]
        if int(f.readline()) != 0:
            sys.exit("expected no boundary triangles")
        nquad = int(f.readline())
        quads = [read_block(f, nquad, 10, int) for _ in range(4)]

    tol = 1.0e-9
    at_bottom = lambda v: abs(z[v - 1]) < tol
    at_layer_top = lambda v: abs(z[v - 1] - FIRST_LAYER_TOP) < tol
    keep_hex = [i for i in range(nhex) if all(z[h[i] - 1] <= FIRST_LAYER_TOP + tol for h in hexes)]
    top_of_layer = {(x[v - 1], y[v - 1]): v for v in range(1, nvert + 1) if at_layer_top(v)}

    patch_quads = []
    first = 0
    for p in range(npatch):
        faces = []
        for q in range(first, last_quad[p]):
            corners = [c[q] for c in quads]
            if all(abs(z[v - 1] - LZ_3D) < tol for v in corners):     # the 3D top: first-layer top
                corners = [top_of_layer[(x[v - 1], y[v - 1])] for v in corners]
            elif not all(z[v - 1] <= FIRST_LAYER_TOP + tol for v in corners):
                continue                                               # a side face of a higher layer
            faces.append(corners)
        patch_quads.append(faces)
        first = last_quad[p]

    used = sorted({h[i] for h in hexes for i in keep_hex})
    new_id = {v: k + 1 for k, v in enumerate(used)}
    if not all(at_bottom(v) or at_layer_top(v) for v in used):
        sys.exit("first layer has vertices off its two planes")
    with open(path, "w") as f:
        f.write(title)
        f.write("%16d\n" % len(used))
        write_block(f, [x[v - 1] for v in used], 5, "%16.9E")
        write_block(f, [y[v - 1] for v in used], 5, "%16.9E")
        write_block(f, [0.0 if at_bottom(v) else dz for v in used], 5, "%16.9E")
        write_block(f, [0] * len(used), 10, "%16d")
        f.write("%16d\n" % 0)
        f.write("%16d\n" % len(keep_hex))
        for h in hexes:
            write_block(f, [new_id[h[i]] for i in keep_hex], 10, "%16d")
        write_block(f, [0] * len(keep_hex), 10, "%16d")
        f.write("%16d\n%16d\n" % (0, 0))
        f.write("%16d\n" % npatch)
        write_block(f, kinds, 10, "%16d")
        write_block(f, [0] * npatch, 10, "%16d")
        ends, total = [], 0
        for faces in patch_quads:
            total += len(faces)
            ends.append(total)
        write_block(f, ends, 10, "%16d")
        f.writelines(names)
        f.write("%16d\n" % 0)
        f.write("%16d\n" % total)
        for corner in range(4):
            write_block(f, [new_id[q[corner]] for faces in patch_quads for q in faces], 10, "%16d")
    return len(used), len(keep_hex), [len(faces) for faces in patch_quads]


def make_bed(dz, seed, path):
    n = round(N_3D * dz / LZ_3D)
    rng = random.Random(seed)
    cell = DP
    grid = {}
    points = []
    r0, r1 = R_IN + 0.5 * DP, R_OUT - 0.5 * DP
    while len(points) < n:
        r = math.sqrt(r0 ** 2 + rng.random() * (r1 ** 2 - r0 ** 2))
        theta = (2.0 * rng.random() - 1.0) * HALF_ANGLE
        if r * math.sin(HALF_ANGLE - abs(theta)) < 0.5 * DP:
            continue                                   # would cut a side wall
        px, py = r * math.cos(theta), r * math.sin(theta)
        i, j = int(math.floor(px / cell)), int(math.floor(py / cell))
        if any((px - a) ** 2 + (py - b) ** 2 < DP ** 2
               for ii in range(i - 1, i + 2) for jj in range(j - 1, j + 2)
               for a, b in grid.get((ii, jj), [])):
            continue
        grid.setdefault((i, j), []).append((px, py))
        points.append((px, py))
    with open(path, "w") as f:
        f.write("! make_wedge_2d.py: dz = %g m, seed = %d, dp = %g m, r %g to %g, +-0.5 deg, "
                "no in-plane overlaps\n" % (dz, seed, DP, R_IN, R_OUT))
        f.write("%d\n" % n)
        for px, py in points:
            f.write("%.16E %.16E %.16E 0 0 0 %.16E\n" % (px, py, 0.5 * dz, DP))
    return n


def make_inp(path, maxtime):
    with open(os.path.join(SOURCE, "wedge.inp")) as f:
        text = f.read()
    text, count = re.subn(r"(?m)^DIMENS .*$", "DIMENS   2      ! one cell thick, BC_VIRTUAL z patches", text)
    if count != 1:
        sys.exit("no DIMENS line in wedge.inp")
    if maxtime:
        for key in ("MAXTIME", "RUNTIME"):
            text, _ = re.subn(r"(?m)^%s .*$" % key, "%s    %s" % (key, maxtime), text)
        for key in ("WRITIME", "PRNTIME"):
            text, _ = re.subn(r"(?m)^%s .*$" % key, "%s    %s" % (key, maxtime), text)
    with open(path, "w") as f:
        f.write(text)


BC = """# BC_SLIPW
PATCH      1 2
NAME       InnerOuterWalls
PARTWALL   1
#

# BC_SLIPW
PATCH      3 4
NAME       Sides
PARTWALL   1
#

# BC_VIRTUAL
PATCH      5 6
NAME       ZVirtual
#

# END
"""


def main():
    dz, seed, out = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    maxtime = sys.argv[4] if len(sys.argv) > 4 else None
    os.makedirs(out, exist_ok=True)
    nv, nh, faces = make_mesh(dz, os.path.join(out, "wedge.hyb.asc"))
    n = make_bed(dz, seed, os.path.join(out, "points.dat"))
    make_inp(os.path.join(out, "wedge.inp"), maxtime)
    with open(os.path.join(out, "wedge.bc"), "w") as f:
        f.write(BC)
    area = 0.5 * 2.0 * HALF_ANGLE * (R_OUT ** 2 - R_IN ** 2)
    vp = math.pi * DP ** 3 / 6.0
    print("%s: %d vertices, %d hexes, patch faces %s" % (out, nv, nh, faces))
    print("  N = %d, phi = %.4f, area fraction = %.3f" % (n, n * vp / (area * dz), n * math.pi * DP ** 2 / 4.0 / area))


if __name__ == "__main__":
    main()
