#!/usr/bin/env python3
"""Refine the one-layer wedge grid k times in r and theta; same geometry, patches and corner order.
Usage: refine.py <source hyb.asc> <k> <output>"""
import math, sys
src, k, dst = sys.argv[1], int(sys.argv[2]), sys.argv[3]
NR, NT = 1499, 24
R0, R1, TH = 0.1, 0.75, math.radians(0.5)
def rows(n, per): return (n + per - 1) // per
def block(f, n, per, conv):
    v = []
    for _ in range(rows(n, per)):
        l = f.readline().rstrip("\n")
        v += [conv(l[i:i+16]) for i in range(0, len(l), 16) if l[i:i+16].strip()]
    assert len(v) == n
    return v
def put(f, v, per, form):
    for i in range(0, len(v), per):
        f.write("".join(form % a for a in v[i:i+per]) + "\n")
with open(src) as f:
    title = f.readline(); nv = int(f.readline())
    x = block(f, nv, 5, float); y = block(f, nv, 5, float); z = block(f, nv, 5, float); block(f, nv, 10, int)
    assert int(f.readline()) == 0
    nh = int(f.readline()); hexes = [block(f, nh, 10, int) for _ in range(8)]; block(f, nh, 10, int)
    assert int(f.readline()) == 0 and int(f.readline()) == 0
    np_ = int(f.readline()); kinds = block(f, np_, 10, int); block(f, np_, 10, int); ends = block(f, np_, 10, int)
    names = [f.readline() for _ in range(np_)]
    assert int(f.readline()) == 0
    nq = int(f.readline()); quads = [block(f, nq, 10, int) for _ in range(4)]
dz = max(z)
assert nh == NR * NT
def ijk(v):
    r, t = math.hypot(x[v-1], y[v-1]), math.atan2(y[v-1], x[v-1])
    i, j = (r - R0) / (R1 - R0) * NR, (t + TH) / (2 * TH) * NT
    assert abs(i - round(i)) < 1e-3 and abs(j - round(j)) < 1e-3, (v, i, j)
    return round(i), round(j), round(z[v-1] / dz)
def offsets(vs):
    p = [ijk(v) for v in vs]; m = [min(c[a] for c in p) for a in range(3)]
    return tuple(tuple(c[a] - m[a] for a in range(3)) for c in p), m
hex_t = {offsets([h[c] for h in hexes])[0] for c in range(nh)}
assert len(hex_t) == 1, hex_t
hex_t = hex_t.pop()
patches, first = [], 0
for p in range(np_):
    ts, axis = set(), set()
    for q in range(first, ends[p]):
        t, m = offsets([c[q] for c in quads]); ts.add(t)
        for a in range(3):
            if all(o[a] == 0 for o in t): axis.add((a, m[a]))
    assert len(ts) == 1 and len(axis) == 1, (p, ts, axis)
    patches.append((ts.pop(), axis.pop())); first = ends[p]
nr, nt = NR * k, NT * k
idx = lambda i, j, l: l * (nr + 1) * (nt + 1) + j * (nr + 1) + i + 1
X, Y, Z = [], [], []
for l in (0, 1):
    for j in range(nt + 1):
        t = -TH + 2 * TH * j / nt
        for i in range(nr + 1):
            r = R0 + (R1 - R0) * i / nr
            X.append(r * math.cos(t)); Y.append(r * math.sin(t)); Z.append(l * dz)
cells = [(i, j) for j in range(nt) for i in range(nr)]
out_quads, out_ends, tot = [[] for _ in range(4)], [], 0
for t, (a, val) in patches:
    if a == 0: base = [(val * k, j, 0) for j in range(nt)]
    elif a == 1: base = [(i, val * k, 0) for i in range(nr)]
    else: base = [(i, j, val) for j in range(nt) for i in range(nr)]
    for b in base:
        for c, o in enumerate(t):
            out_quads[c].append(idx(b[0] + o[0], b[1] + o[1], b[2] + o[2]))
    tot += len(base); out_ends.append(tot)
with open(dst, "w") as f:
    f.write(title); f.write("%16d\n" % len(X))
    for v in (X, Y, Z): put(f, v, 5, "%16.9E")
    put(f, [0] * len(X), 10, "%16d")
    f.write("%16d\n%16d\n" % (0, len(cells)))
    for o in hex_t:
        put(f, [idx(i + o[0], j + o[1], o[2]) for i, j in cells], 10, "%16d")
    put(f, [0] * len(cells), 10, "%16d")
    f.write("%16d\n%16d\n%16d\n" % (0, 0, np_))
    put(f, kinds, 10, "%16d"); put(f, [0] * np_, 10, "%16d"); put(f, out_ends, 10, "%16d")
    f.writelines(names); f.write("%16d\n%16d\n" % (0, tot))
    for q in out_quads: put(f, q, 10, "%16d")
print(dst, "points", len(X), "cells", len(cells), "patch quads", out_ends, "corners ok")
