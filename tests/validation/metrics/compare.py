"""Measures and pass/fail checks for rung 1 (Sod, Lax), from any solver's output in plain columns.

Result file format (dimensionless units, as in ../problems/):

    # t 0.2                  end time reached
    # steps 1400             number of time steps taken
    # columns x rho u p v w  names of the columns below; v and w only in 2D and 3D
    0.00125 1.0 0.0 1.0 0.0 0.0
    ...

One line per cell, in any order, on a uniform grid over 0 <= x <= 1 (in 2D and 3D every cell, with
equal volumes). Lines starting with # other than the three above are ignored.

    python3 compare.py sod --order 1 sod_n200.txt sod_n400.txt sod_n800.txt
    python3 compare.py sod --same ranks sod_n200_r1.txt sod_n200_r4.txt
    python3 compare.py sod --same dims  sod_n800_2d.txt sod_n800_3d.txt

The limits are those of criteria.md.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "exact"))
from riemann import PROBLEMS, problem  # noqa: E402
import godunov  # noqa: E402

LIMITS = dict(conservation=1e-10, end_cells=1e-12, ranks=1e-10, transverse=1e-10,
              rate={1: 0.5, 2: 0.67}, shock_cells=1.5, godunov=1.5, dims=0.01)


def read(path):
    """Header values and columns of a result file."""
    head, rows = {}, []
    for line in open(path):
        if line.startswith("#"):
            parts = line[1:].split()
            if parts and parts[0] in ("t", "steps"):
                head[parts[0]] = float(parts[1])
            elif parts and parts[0] == "columns":
                head["columns"] = parts[1:]
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    for key in ("t", "steps", "columns"):
        if key not in head:
            raise ValueError("%s: no '# %s' line" % (path, key))
    data = np.array(rows)
    cols = {name: data[:, i] for i, name in enumerate(head["columns"])}
    for name in ("x", "rho", "u", "p"):
        if name not in cols:
            raise ValueError("%s: no column %s" % (path, name))
    return head, cols


def columns(x, q):
    """Unique cell-centre positions and the mean of q over each cross-section."""
    key = np.round(x, 9)
    xs, inv = np.unique(key, return_inverse=True)
    return xs, np.bincount(inv, q) / np.bincount(inv)


def crossing(x, f, level, lo, hi):
    """Rightmost place in [lo, hi] where f crosses level, interpolated linearly."""
    m = np.where((x >= lo) & (x <= hi))[0]
    above = f[m] > level
    k = np.where(above[:-1] != above[1:])[0]
    if len(k) == 0:
        return float("nan")
    i, j = m[k[-1]], m[k[-1] + 1]
    return float(x[i] + (level - f[i]) / (f[j] - f[i]) * (x[j] - x[i]))


def flux(state, g):
    r, u, p = state
    E = p / (g - 1) + 0.5 * r * u * u
    return np.array([r * u, r * u * u + p, u * (E + p)])


def totals(state, g):
    r, u, p = state
    return np.array([r, r * u, p / (g - 1) + 0.5 * r * u * u])


def measure(name, path):
    """All single-run measures of criteria.md for one result file."""
    head, c = read(path)
    d = PROBLEMS[name]
    ex, x0, t_end = problem(name)
    g, t = d["gamma"], head["t"]
    w = ex.waves()
    out = dict(path=path, t=t, steps=int(head["steps"]), cells=len(c["x"]))
    if abs(t / t_end - 1) > 1e-10:
        out["warning"] = "end time %.10g, problem end time %.10g" % (t, t_end)
    xs, rho = columns(c["x"], c["rho"])
    n = len(xs)
    out["n"] = n
    exact = ex.sample(c["x"], t, x0)
    out["L1"] = [float(np.abs(c[q] - e).mean()) for q, e in zip(("rho", "u", "p"), exact)]
    # shock: rightmost crossing of the mid level between the contact and the shock and the end
    xc, xsh = x0 + w["u_star"] * t, x0 + w["shock_R"] * t
    level = 0.5 * (w["rho_star_R"] + d["right"][0])
    out["shock_cells"] = (crossing(xs, rho, level, 0.5 * (xc + xsh), 1.0) - xsh) * n
    # conservation: cell totals against the initial totals plus t times the exact boundary fluxes
    U = np.array([c["rho"], c["rho"] * c["u"], c["p"] / (g - 1) + 0.5 * c["rho"] * c["u"] ** 2])
    if "v" in c:
        U[2] += 0.5 * c["rho"] * (c["v"] ** 2 + (c["w"] ** 2 if "w" in c else 0))
    got = U.mean(axis=1)
    want = x0 * totals(d["left"], g) + (1 - x0) * totals(d["right"], g) + t * (flux(d["left"], g) - flux(d["right"], g))
    out["conservation"] = [float(e) for e in (got - want) / np.abs(want)]
    # end cells: do the first and last columns still hold the initial states?
    a_max = max(np.sqrt(g * s[2] / s[0]) for s in (d["left"], d["right"]))
    end = 0.0
    for xe, s in ((xs[0], d["left"]), (xs[-1], d["right"])):
        m = np.round(c["x"], 9) == xe
        end = max(end, np.abs(c["rho"][m] / s[0] - 1).max(), np.abs(c["p"][m] / s[2] - 1).max(),
                  np.abs(c["u"][m] - s[1]).max() / a_max)
    out["end_cells"] = float(end)
    if "v" in c:
        vt = np.sqrt(np.mean(c["v"] ** 2 + (c["w"] ** 2 if "w" in c else 0)))
        out["transverse"] = float(vt / a_max)
    return out


def report(name, order, paths):
    """Measures, rates and verdicts for runs of one scheme on several grids."""
    runs = sorted((measure(name, p) for p in paths), key=lambda r: r["n"])
    rows, fails = [], 0

    def verdict(what, value, ok, shown=None):
        nonlocal fails
        fails += ok is False
        rows.append((what, shown if shown is not None else "%.3g" % value,
                     {True: "pass", False: "FAIL", None: "reported"}[ok]))

    print("%s, order %d" % (name, order))
    print("%6s %7s %11s %11s %11s %8s %26s %9s %9s" % ("n", "steps", "L1 rho", "L1 u", "L1 p", "shock", "conservation", "end cells", "transv."))
    for r in runs:
        print("%6d %7d %11.4e %11.4e %11.4e %8.3f %8.1e %8.1e %8.1e %9.1e %9s" % (
            r["n"], r["steps"], *r["L1"], r["shock_cells"], *[abs(e) for e in r["conservation"]], r["end_cells"],
            "%.1e" % r["transverse"] if "transverse" in r else "-"))
        if "warning" in r:
            print("       warning: " + r["warning"])
    for r in runs:
        tag = "n %d " % r["n"]
        cons = max(abs(e) for e in r["conservation"])
        clean = r["end_cells"] <= LIMITS["end_cells"]
        verdict(tag + "conservation" + ("" if clean else " (end cells disturbed)"), cons,
                cons < LIMITS["conservation"] if clean else None)
        verdict(tag + "shock position, cells", abs(r["shock_cells"]), abs(r["shock_cells"]) <= LIMITS["shock_cells"])
        if "transverse" in r:
            verdict(tag + "transverse velocity / a", r["transverse"], r["transverse"] < LIMITS["transverse"])
        if order == 1:
            ref = godunov.run(name, r["n"], steps=r["steps"])
            verdict(tag + "L1(rho) / Godunov (%d steps)" % r["steps"], r["L1"][0] / ref["L1"][0],
                    r["L1"][0] / ref["L1"][0] <= LIMITS["godunov"])
    for a, b in zip(runs, runs[1:]):
        rates = [np.log2(p / q) for p, q in zip(a["L1"], b["L1"])]
        finest = b is runs[-1]
        verdict("rate rho %d-%d" % (a["n"], b["n"]), rates[0],
                rates[0] >= LIMITS["rate"][order] if finest else None)
        verdict("rates u, p %d-%d" % (a["n"], b["n"]), 0, None, "%.2f %.2f" % tuple(rates[1:]))
    print()
    for what, value, ok in rows:
        print("  %-44s %12s  %s" % (what, value, ok))
    print("%d failed" % fails)
    return fails


def same(name, kind, path_a, path_b):
    """Rank independence (ranks) or 2D/3D consistency (dims) between two runs."""
    if kind == "ranks":
        (_, a), (_, b) = read(path_a), read(path_b)
        ka, kb = [np.lexsort(np.round(np.array([c[k] for k in ("z", "y", "x") if k in c]), 9)) for c in (a, b)]
        diff = max(float(np.abs(b[q][kb] / a[q][ka] - 1).max()) for q in ("rho", "p"))
        a_max = max(np.sqrt(PROBLEMS[name]["gamma"] * s[2] / s[0]) for s in (PROBLEMS[name]["left"], PROBLEMS[name]["right"]))
        diff = max(diff, float(np.abs(b["u"][kb] - a["u"][ka]).max() / a_max))
        ok = diff < LIMITS["ranks"]
        print("rank independence: largest relative difference %.2e  %s" % (diff, "pass" if ok else "FAIL"))
    else:
        la, lb = measure(name, path_a)["L1"][0], measure(name, path_b)["L1"][0]
        diff = abs(la - lb) / lb
        ok = diff < LIMITS["dims"]
        print("dimensions: |L1 a - L1 b| / L1 b = %.2e  %s" % (diff, "pass" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) >= 4 and args[0] in PROBLEMS and args[1] == "--order" and args[2] in ("1", "2"):
        sys.exit(1 if report(args[0], int(args[2]), args[3:]) else 0)
    if len(args) == 5 and args[0] in PROBLEMS and args[1] == "--same" and args[2] in ("ranks", "dims"):
        sys.exit(same(args[0], args[2], args[3], args[4]))
    sys.exit(__doc__)
