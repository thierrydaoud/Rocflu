"""Godunov's first-order scheme with the exact Riemann solver: the reference for first-order schemes.

1D Euler equations, ideal gas, uniform cells on 0 <= x <= 1, zero-gradient ends, forward Euler in
time. The face states use a vectorised exact solver, written separately from riemann.py so that
each checks the other (see --test). Needs only numpy.

    python3 godunov.py sod 800                 # CFL 0.5 in the 1D sense: dt = 0.5 dx / max(|u| + a)
    python3 godunov.py sod 800 --steps 2400    # 2400 equal steps, to match another solver's step count
    python3 godunov.py --test

Prints x rho u p at the cell centres in the format metrics/compare.py reads, with the L1 errors in a comment.
"""
import sys
import numpy as np
from riemann import PROBLEMS, Riemann, cell_centres, problem


def _f(p, r, pk, a, g):
    """Toro's f_K(p) and its derivative, for arrays."""
    A, B = 2 / ((g + 1) * r), (g - 1) / (g + 1) * pk
    q = np.sqrt(A / (p + B))
    ratio = np.maximum(p / pk, 1e-300)
    shock = p > pk
    f = np.where(shock, (p - pk) * q, 2 * a / (g - 1) * (ratio ** ((g - 1) / (2 * g)) - 1))
    df = np.where(shock, q * (1 - 0.5 * (p - pk) / (p + B)), ratio ** (-(g + 1) / (2 * g)) / (r * a))
    return f, df


def face_state(rl, ul, pl, rr, ur, pr, g):
    """Exact Riemann solution at x/t = 0 for arrays of left and right states: rho, u, p."""
    al, ar = np.sqrt(g * pl / rl), np.sqrt(g * pr / rr)
    e = (g - 1) / (2 * g)  # two-rarefaction guess: exact when both waves are rarefactions
    p = ((al + ar - 0.5 * (g - 1) * (ur - ul)) / (al / pl ** e + ar / pr ** e)) ** (1 / e)
    for _ in range(100):
        fl, dl = _f(p, rl, pl, al, g)
        fr, dr = _f(p, rr, pr, ar, g)
        new = np.maximum(p - (fl + fr + ur - ul) / (dl + dr), 0.1 * p)
        done = np.all(np.abs(new - p) <= 1e-12 * 0.5 * (new + p))
        p = new
        if done:
            break
    else:
        raise RuntimeError("face pressure did not converge")
    us = 0.5 * (ul + ur) + 0.5 * (_f(p, rr, pr, ar, g)[0] - _f(p, rl, pl, al, g)[0])
    # Take the side of the contact that holds x/t = 0, mirrored so it is always a left side.
    left = us >= 0
    s = np.where(left, 1.0, -1.0)
    r0, u0, p0, a0 = np.where(left, rl, rr), s * np.where(left, ul, ur), np.where(left, pl, pr), np.where(left, al, ar)
    u_star = s * us
    k = (g - 1) / (g + 1)
    shock = p > p0
    S = u0 - a0 * np.sqrt((g + 1) / (2 * g) * p / p0 + (g - 1) / (2 * g))
    head, tail = u0 - a0, u_star - a0 * (p / p0) ** ((g - 1) / (2 * g))
    c = np.maximum(2 / (g + 1) + (g - 1) / ((g + 1) * a0) * u0, 1e-300)  # inside the fan at x/t = 0
    fan = (r0 * c ** (2 / (g - 1)), 2 / (g + 1) * (a0 + (g - 1) / 2 * u0), p0 * c ** (2 * g / (g - 1)))
    ahead = np.where(shock, S >= 0, head >= 0)
    behind = np.where(shock, S < 0, tail <= 0)
    star = (np.where(shock, r0 * (p / p0 + k) / (k * p / p0 + 1), r0 * (p / p0) ** (1 / g)), u_star, p)
    out = []
    for q0, qs, qf in zip((r0, u0, p0), star, fan):
        out.append(np.where(ahead, q0, np.where(behind, qs, qf)))
    return out[0], s * out[1], out[2]


def solve(left, right, gamma, x0, t_end, n, cfl=0.5, steps=None):
    """Run to t_end. With steps, take that many equal steps; else dt = cfl dx / max(|u| + a).
    Returns x, rho, u, p, number of steps, largest 1D CFL number used, conservation error."""
    g, dx = gamma, 1.0 / n
    x = cell_centres(n)
    W = [np.where(x < x0, a, b).astype(float) for a, b in zip(left, right)]
    U = np.array([W[0], W[0] * W[1], W[2] / (g - 1) + 0.5 * W[0] * W[1] ** 2])
    total0, inflow = U.sum(axis=1) * dx, np.zeros(3)
    t, taken, cfl_max = 0.0, 0, 0.0
    while t < t_end * (1 - 1e-14):
        r, u = U[0], U[1] / U[0]
        p = (g - 1) * (U[2] - 0.5 * r * u ** 2)
        speed = np.max(np.abs(u) + np.sqrt(g * p / r))
        dt = t_end / steps if steps else min(cfl * dx / speed, t_end - t)
        cfl_max = max(cfl_max, dt * speed / dx)
        if cfl_max > 1:
            raise ValueError("time step too large: 1D CFL number %.3f" % cfl_max)
        ext = [np.concatenate(([a[0]], a, [a[-1]])) for a in (r, u, p)]  # zero-gradient ends
        rf, uf, pf = face_state(ext[0][:-1], ext[1][:-1], ext[2][:-1], ext[0][1:], ext[1][1:], ext[2][1:], g)
        E = pf / (g - 1) + 0.5 * rf * uf ** 2
        F = np.array([rf * uf, rf * uf ** 2 + pf, uf * (E + pf)])
        U = U - dt / dx * (F[:, 1:] - F[:, :-1])
        inflow += dt * (F[:, 0] - F[:, -1])
        t += dt
        taken += 1
    r, u = U[0], U[1] / U[0]
    p = (g - 1) * (U[2] - 0.5 * r * u ** 2)
    error = (U.sum(axis=1) * dx - total0 - inflow) / np.abs(U.sum(axis=1) * dx)
    return x, r, u, p, taken, cfl_max, error


def run(name, n, cfl=0.5, steps=None):
    """Solve a named problem; also returns the L1 errors of rho, u, p against the exact solution."""
    d = PROBLEMS[name]
    x, r, u, p, taken, cfl_max, error = solve(d["left"], d["right"], d["gamma"], d["x0"], d["t"], n, cfl, steps)
    exact = problem(name)[0].sample(x, d["t"], d["x0"])
    L1 = [float(np.abs(q - e).sum() / n) for q, e in zip((r, u, p), exact)]
    return dict(x=x, rho=r, u=u, p=p, steps=taken, cfl_max=cfl_max, conservation=error, L1=L1)


def _test():
    fails = 0

    def check(what, ok, detail):
        nonlocal fails
        fails += not ok
        print("%-4s %-44s %s" % ("ok" if ok else "FAIL", what, detail))

    # 1. face_state against riemann.Riemann at x/t = 0, random states of all four wave patterns
    rng = np.random.default_rng(1)
    m = 2000
    rl, rr = rng.uniform(0.1, 10, m), rng.uniform(0.1, 10, m)
    pl, pr = 10 ** rng.uniform(-2, 3, m), 10 ** rng.uniform(-2, 3, m)
    ul, ur = rng.uniform(-3, 3, m), rng.uniform(-3, 3, m)
    keep = 2 / 0.4 * (np.sqrt(1.4 * pl / rl) + np.sqrt(1.4 * pr / rr)) > (ur - ul) * 1.05  # no vacuum
    got = np.array(face_state(rl[keep], ul[keep], pl[keep], rr[keep], ur[keep], pr[keep], 1.4))
    want = np.array([np.ravel(Riemann((a, b, c), (d, e, f)).sample(0.0, 1.0))
                     for a, b, c, d, e, f in zip(rl[keep], ul[keep], pl[keep], rr[keep], ur[keep], pr[keep])]).T
    scale = np.array([np.maximum(rl, rr), np.sqrt(1.4 * np.maximum(pl / rl, pr / rr)), np.maximum(pl, pr)])[:, keep]
    worst = float(np.max(np.abs(got - want) / scale))
    check("face states vs riemann.py (%d random pairs)" % keep.sum(), worst < 1e-10, "max scaled diff %.1e" % worst)
    # 2. conservation and convergence on Sod and Lax
    for name in PROBLEMS:
        res = [run(name, n) for n in (200, 400, 800, 1600)]
        cons = max(float(np.abs(r["conservation"]).max()) for r in res)
        check("%s conservation, n 200-1600" % name, cons < 1e-12, "max %.1e" % cons)
        rates = [np.log2(a["L1"][0] / b["L1"][0]) for a, b in zip(res, res[1:])]
        check("%s L1(rho) rates" % name, min(rates[1:]) >= 0.5, " ".join("%.2f" % q for q in rates))
    print("%d failed" % fails)
    return fails


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["--test"]:
        sys.exit(1 if _test() else 0)
    if len(args) not in (2, 4) or args[0] not in PROBLEMS or (len(args) == 4 and args[2] not in ("--steps", "--cfl")):
        sys.exit(__doc__)
    opt = {} if len(args) == 2 else ({"steps": int(args[3])} if args[2] == "--steps" else {"cfl": float(args[3])})
    res = run(args[0], int(args[1]), **opt)
    print("# Godunov %s n %s: %d steps, max 1D CFL %.4f; L1 rho %.6e u %.6e p %.6e; conservation %.1e %.1e %.1e"
          % (args[0], args[1], res["steps"], res["cfl_max"], *res["L1"], *res["conservation"]))
    print("# t %.17g\n# steps %d\n# columns x rho u p" % (PROBLEMS[args[0]]["t"], res["steps"]))
    for row in zip(res["x"], res["rho"], res["u"], res["p"]):
        print("%.10f %.12e %.12e %.12e" % row)
