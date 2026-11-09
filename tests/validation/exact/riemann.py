"""Exact solution of the Riemann problem for the 1D Euler equations of an ideal gas.

Follows E. F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, 3rd ed.
(Springer 2009), ch. 4: all four wave patterns (shock or rarefaction on either side).
Needs only numpy.

As a script, prints the exact solution at cell centres in plain columns:

    python3 riemann.py sod 200 > sod_exact_n200.txt      # x rho u p
    python3 riemann.py --test                            # self-test
"""
import sys
import numpy as np

# Problems on 0 <= x <= 1 with the jump at x0; see ../problems/<name>.md.
PROBLEMS = {
    "sod": dict(left=(1.0, 0.0, 1.0), right=(0.125, 0.0, 0.1), gamma=1.4, x0=0.5, t=0.2),
    "lax": dict(left=(0.445, 0.698, 3.528), right=(0.5, 0.0, 0.571), gamma=1.4, x0=0.5, t=0.14),
}


class Riemann:
    """Exact solution for left and right states (rho, u, p) and ratio of specific heats gamma."""

    def __init__(self, left, right, gamma=1.4, tol=1e-12):
        self.left, self.right, self.g = tuple(map(float, left)), tuple(map(float, right)), float(gamma)
        rl, ul, pl = self.left
        rr, ur, pr = self.right
        g = self.g
        self.al, self.ar = np.sqrt(g * pl / rl), np.sqrt(g * pr / rr)
        if 2 / (g - 1) * (self.al + self.ar) <= ur - ul:
            raise ValueError("the initial states create a vacuum")
        self.p = self._star_pressure(tol)
        self.u = 0.5 * (ul + ur) + 0.5 * (self._f(self.p, *self._side("R"))[0] - self._f(self.p, *self._side("L"))[0])

    def _side(self, s):
        r, _, p = self.left if s == "L" else self.right
        return r, p, (self.al if s == "L" else self.ar)

    def _f(self, p, r, pk, a):
        """Toro's f_K(p) and its derivative: velocity change across the wave on side K."""
        g = self.g
        if p > pk:  # shock
            A, B = 2 / ((g + 1) * r), (g - 1) / (g + 1) * pk
            q = np.sqrt(A / (p + B))
            return (p - pk) * q, q * (1 - 0.5 * (p - pk) / (p + B))
        e = (g - 1) / (2 * g)  # rarefaction
        return 2 * a / (g - 1) * ((p / pk) ** e - 1), (p / pk) ** (-(g + 1) / (2 * g)) / (r * a)

    def _star_pressure(self, tol):
        """Newton iteration on f_L(p) + f_R(p) + u_R - u_L = 0, from the two-rarefaction guess."""
        rl, ul, pl = self.left
        rr, ur, pr = self.right
        g, e = self.g, (self.g - 1) / (2 * self.g)
        p = ((self.al + self.ar - 0.5 * (g - 1) * (ur - ul)) / (self.al / pl ** e + self.ar / pr ** e)) ** (1 / e)
        for _ in range(100):
            fl, dl = self._f(p, *self._side("L"))
            fr, dr = self._f(p, *self._side("R"))
            new = max(p - (fl + fr + ur - ul) / (dl + dr), 0.1 * p)
            if abs(new - p) <= tol * 0.5 * (new + p):
                return new
            p = new
        raise RuntimeError("star pressure did not converge")

    def waves(self):
        """Star states and wave speeds, as a dict."""
        g = self.g
        out = dict(p_star=self.p, u_star=self.u)
        for s, sign in (("L", -1), ("R", 1)):
            r, pk, a = self._side(s)
            uk = (self.left if s == "L" else self.right)[1]
            if self.p > pk:
                k = (g - 1) / (g + 1)
                out["rho_star_" + s] = r * (self.p / pk + k) / (k * self.p / pk + 1)
                out["shock_" + s] = uk + sign * a * np.sqrt((g + 1) / (2 * g) * self.p / pk + (g - 1) / (2 * g))
            else:
                out["rho_star_" + s] = r * (self.p / pk) ** (1 / g)
                a_star = a * (self.p / pk) ** ((g - 1) / (2 * g))
                out["head_" + s] = uk + sign * a
                out["tail_" + s] = self.u + sign * a_star
        return out

    def sample(self, x, t, x0=0.0):
        """rho, u, p at positions x (scalar or array) at time t > 0."""
        s = (np.atleast_1d(np.asarray(x, dtype=float)) - x0) / t
        out = np.empty((3, s.size))
        for i, si in enumerate(s):
            out[:, i] = self._at(si)
        return out[0], out[1], out[2]

    def _at(self, s):
        g = self.g
        left = s <= self.u
        r, uk, pk = self.left if left else self.right
        a = self.al if left else self.ar
        sign = -1 if left else 1
        if self.p > pk:  # shock on this side
            k = (g - 1) / (g + 1)
            S = uk + sign * a * np.sqrt((g + 1) / (2 * g) * self.p / pk + (g - 1) / (2 * g))
            if sign * (s - S) >= 0:
                return r, uk, pk
            return r * (self.p / pk + k) / (k * self.p / pk + 1), self.u, self.p
        head = uk + sign * a  # rarefaction on this side
        tail = self.u + sign * a * (self.p / pk) ** ((g - 1) / (2 * g))
        if sign * (s - head) >= 0:
            return r, uk, pk
        if sign * (s - tail) <= 0:
            return r * (self.p / pk) ** (1 / g), self.u, self.p
        c = 2 / (g + 1) - sign * (g - 1) / ((g + 1) * a) * (uk - s)
        return r * c ** (2 / (g - 1)), 2 / (g + 1) * (-sign * a + (g - 1) / 2 * uk + s), pk * c ** (2 * g / (g - 1))


def problem(name):
    """The Riemann solution for a named problem, plus its jump position and end time."""
    d = PROBLEMS[name]
    return Riemann(d["left"], d["right"], d["gamma"]), d["x0"], d["t"]


def cell_centres(n):
    return (np.arange(n) + 0.5) / n


def _test():
    """Checks against the star values in ../problems and in Toro (2009), Table 4.3."""
    fails = 0

    def check(what, got, want, rel=0.0, absolute=0.0):
        nonlocal fails
        ok = abs(got - want) <= max(rel * abs(want), absolute)
        fails += not ok
        print("%-4s %-28s %.6g (expected %.6g)" % ("ok" if ok else "FAIL", what, got, want))

    star = {"sod": (0.30313, 0.92745, 0.42632, 0.26557, 1.75216),
            "lax": (2.46610, 1.52872, 0.34457, 1.30408, 2.47932)}
    for name, (p, u, rl, rr, S) in star.items():
        w = problem(name)[0].waves()
        for key, want in (("p_star", p), ("u_star", u), ("rho_star_L", rl), ("rho_star_R", rr), ("shock_R", S)):
            check("%s %s" % (name, key), w[key], want, absolute=5e-6)  # values given to 5 decimals
    # Toro, Table 4.3: tests 1-5 (p*, u*); test 2 is the 123 problem (two rarefactions)
    toro = [((1.0, 0.0, 1.0), (0.125, 0.0, 0.1), 0.30313, 0.92745),
            ((1.0, -2.0, 0.4), (1.0, 2.0, 0.4), 0.00189, 0.0),
            ((1.0, 0.0, 1000.0), (1.0, 0.0, 0.01), 460.894, 19.5975),
            ((1.0, 0.0, 0.01), (1.0, 0.0, 100.0), 46.0950, -6.19633),
            ((5.99924, 19.5975, 460.894), (5.99242, -6.19633, 46.0950), 1691.64, 8.68975)]
    for i, (L, R, p, u) in enumerate(toro, 1):
        r = Riemann(L, R)
        check("Toro test %d p_star" % i, r.p, p, rel=1e-5, absolute=5e-6)  # 6 significant figures
        check("Toro test %d u_star" % i, r.u, u, rel=1e-5, absolute=5e-6)
    # Rankine-Hugoniot across the Sod shock and continuity at every wave edge
    r, x0, t = problem("sod")
    w, g = r.waves(), r.g
    rho1, u1, p1 = r.right
    rho2, u2, p2 = w["rho_star_R"], r.u, r.p
    S = w["shock_R"]
    E = lambda rho, u, p: p / (g - 1) + 0.5 * rho * u * u
    for what, a, b in (("mass", rho2 * (u2 - S), rho1 * (u1 - S)),
                       ("momentum", rho2 * u2 * (u2 - S) + p2, rho1 * u1 * (u1 - S) + p1),
                       ("energy", E(rho2, u2, p2) * (u2 - S) + p2 * u2, E(rho1, u1, p1) * (u1 - S) + p1 * u1)):
        check("Sod shock jump in " + what, a, b, rel=1e-12)
    for edge in ("head_L", "tail_L"):
        x = x0 + w[edge] * t
        a, b = np.array(r.sample([x - 1e-9, x + 1e-9], t, x0))[:, :2].T
        check("Sod continuity at " + edge, np.abs(a - b).max(), 0.0, absolute=1e-6)
    print("%d failed" % fails)
    return fails


if __name__ == "__main__":
    if sys.argv[1:] == ["--test"]:
        sys.exit(1 if _test() else 0)
    if len(sys.argv) != 3 or sys.argv[1] not in PROBLEMS:
        sys.exit(__doc__)
    r, x0, t = problem(sys.argv[1])
    x = cell_centres(int(sys.argv[2]))
    rho, u, p = r.sample(x, t, x0)
    print("# %s exact at t = %g: x rho u p" % (sys.argv[1], t))
    for row in zip(x, rho, u, p):
        print("%.10f %.12e %.12e %.12e" % row)
