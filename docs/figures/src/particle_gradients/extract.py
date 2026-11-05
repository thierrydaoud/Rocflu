#!/usr/bin/env python3
"""Write data.json for plot.py from the diagnostic runs (HiPerGator scratch, 2026-11-05).

The runs used a diagnostic build (main e7f09523 plus a dump of each gradient variant on the
same flow state, a particle force log and a switch for LS-only particle gradients) and MFC
master 7b16d32f with a dump of its reconstructed face states. Paths below are where they ran.
"""
import glob, json, os, sys
import numpy as np

RUNS = "/blue/bala1s/t.daoud/scratch_claude/grad_compare/runs"
sys.path.insert(0, "/blue/bala1s/t.daoud/scratch_claude/numerics/scripts")
from riemann import Riemann  # exact Riemann solver used for the numerics validation

HERE = os.path.dirname(os.path.abspath(__file__))
LEFT, RIGHT, X0 = (4.0, 0.0, 6.0e5), (0.9721771595969, 0.0, 82700.0), -0.02
R = Riemann(LEFT, RIGHT)


def load_rocflu(run, k):
    """Cells of one dump: coordinates, volume, rho_g, p, phi_g, and gradients
    [cell, field (rho_g, p), variant (LS, LS+WENO, LS+lim, LS+WENO+lim), component]."""
    rows, t = [], None
    for f in sorted(glob.glob("%s/%s/gradcmp_r*_k%06d.bin" % (RUNS, run, k))):
        b = open(f, "rb").read()
        n = int(np.frombuffer(b[:4], "i4")[0]); t = float(np.frombuffer(b[4:12], "f8")[0])
        rows.append(np.frombuffer(b[12:], "f8").reshape(n, 37))
    a = np.vstack(rows)
    return t, a, a[:, 7:31].reshape(-1, 2, 4, 3)


def planes(x, dx, *fields):
    ix = np.round((x - x.min()) / dx).astype(int); c = np.bincount(ix)
    return [np.bincount(ix, f) / c for f in (x,) + fields]


def exact_avg(x, t, dx):
    s = np.linspace(-0.5, 0.5, 41)
    r = np.mean([R.sample(x + q * dx, t, X0)[0] for q in s], axis=0)
    p = np.mean([R.sample(x + q * dx, t, X0)[2] for q in s], axis=0)
    return r, p, np.gradient(r, dx), np.gradient(p, dx)


out = {"states": {"left": LEFT, "right": RIGHT, "x0": X0}}

# Riemann tube, Rocflu, 1200 x 4 x 4 cells of 0.25 mm, HLLC, RECONST 2, t = 60 us
t, a, g = load_rocflu("tube2", 3)
x, rho, p, gr_ls, gr_w, gp_ls, gp_w = planes(a[:, 0], 2.5e-4, a[:, 4], a[:, 5], g[:, 0, 0, 0],
                                             g[:, 0, 3, 0], g[:, 1, 0, 0], g[:, 1, 3, 0])
er, ep, egr, egp = exact_avg(x, t, 2.5e-4)
m = (x > -0.055) & (x < 0.02)
out["tube"] = {"t": t, "shock": X0 + R.shock * t, "contact": X0 + R.us * t,
               "head": X0 + R.head * t, "tail": X0 + R.tail * t,
               "x": x[m].tolist(), "p": p[m].tolist(), "rho": rho[m].tolist(),
               "dpdx": {"LS": gp_ls[m].tolist(), "LS+WENO": gp_w[m].tolist(), "exact": egp[m].tolist()},
               "drhodx": {"LS": gr_ls[m].tolist(), "LS+WENO": gr_w[m].tolist(), "exact": egr[m].tolist()},
               "limiter_cells_changed": int((g[:, :, 3] != g[:, :, 1]).any(axis=(1, 2)).sum())}

# MFC, same problem, 1200 x 26 cells, WENO5 mapped, HLLC, RK3, CFL 0.4
def load_mfc(run):
    f = "%s/%s/gradcmp_mfc.dat" % (RUNS, run)
    tt = float(open(f).readline().split()[2]); d = np.loadtxt(f)
    xs = np.unique(d[:, 0]); j = np.searchsorted(xs, d[:, 0]); c = np.bincount(j)
    av = lambda v: np.bincount(j, v) / c
    pm, dxm = av(d[:, 3]), av(d[:, 8])
    fd = np.r_[np.nan, (pm[2:] - pm[:-2]) / (2 * dxm[1:-1]), np.nan]
    return tt, xs, av((d[:, 5] - d[:, 4]) / d[:, 8]), fd
win = lambda xx, tt: np.abs(xx - (X0 + R.shock * tt)) < 2.5e-3
codes = {"t": t}
mm = win(x, t)
codes["Rocflu LS"] = [x[mm].tolist(), gp_ls[mm].tolist()]
codes["Rocflu LS+WENO"] = [x[mm].tolist(), gp_w[mm].tolist()]
codes["exact"] = [x[mm].tolist(), egp[mm].tolist()]
for run, label in (("mfc_tube", "MFC WENO faces, mp_weno"), ("mfc_tube_nomp", "MFC WENO faces, no mp_weno")):
    tt, xs, gw, fd = load_mfc(run); w = win(xs, tt)
    codes[label] = [xs[w].tolist(), gw[w].tolist()]
    if run == "mfc_tube":
        codes["MFC central difference"] = [xs[w].tolist(), np.nan_to_num(fd[w]).tolist()]
out["codes"] = codes

# Planar curtain copy (1.8M cells, bed 0-2 mm), plane averages near the bed
curt = {}
for k in (2, 3, 6):
    t, a, g = load_rocflu("curtain", k)
    u = np.unique(a[:, 0]); d = np.diff(u); dxm = np.median(d[d > 1e-6])
    ix = np.round((a[:, 0] - a[:, 0].min()) / dxm).astype(int); c = np.bincount(ix)
    av = lambda v: np.bincount(ix, v) / c
    xx = av(a[:, 0]); m = (xx > -3e-3) & (xx < 5e-3)
    phip = 1 - a[:, 6]; w = phip * a[:, 3]; bed = phip > 1e-3
    ratio = lambda fi: float((w[bed] * np.linalg.norm(g[bed, fi, 1], axis=1)).sum()
                             / (w[bed] * np.linalg.norm(g[bed, fi, 0], axis=1)).sum())
    curt[str(k)] = {"t": t, "x": xx[m].tolist(), "phip": av(phip)[m].tolist(), "p": av(a[:, 5])[m].tolist(),
                    "dpdx_LS": av(g[:, 1, 0, 0])[m].tolist(), "dpdx_WENO": av(g[:, 1, 3, 0])[m].tolist(),
                    "pw_ratio_p": ratio(1), "pw_ratio_rho": ratio(0),
                    "limiter_cells_changed": int((g[:, :, 3] != g[:, :, 1]).any(axis=(1, 2)).sum())}
out["curtain"] = curt

# Fixed particles in the tube, one-way, PGFORCE 1 + AMFORCE 1 (force log at stage 1)
def forces(tag):
    a = np.vstack([np.loadtxt(f, ndmin=2) for f in glob.glob("%s/force_tube_%s/fort.91*" % (RUNS, tag))
                   if open(f).read().strip()])
    res = {}
    for xp in np.unique(np.round(a[:, 1], 7)):
        b = a[np.abs(a[:, 1] - xp) < 1e-8]; b = b[np.argsort(b[:, 0])]
        _, iu = np.unique(b[:, 0], return_index=True); res[round(float(xp), 6)] = b[iu]
    return res
B, L = forces("base"), forces("ls")
vp = 4.0 / 3.0 * np.pi * (0.5 * 1.15e-4) ** 3; h = 1.25e-4
fz = {}
for xp in sorted(B):
    b, l = B[xp], L[xp]
    tf = np.linspace(1e-9, b[:, 0].max(), 200001)
    pp = lambda q: R.sample(np.full_like(tf, q), tf, X0)[2]
    ex = -vp * (pp(xp + h) - pp(xp - h)) / (2 * h)
    sl = slice(None, None, max(1, len(b) // 400))
    fz[str(xp)] = {"t": b[sl, 0].tolist(), "pg_old": b[sl, 6].tolist(), "pg_ls": l[sl, 6].tolist(),
                   "am_old": b[sl, 5].tolist(), "am_ls": l[sl, 5].tolist(),
                   "t_exact": tf[::500].tolist(), "pg_exact": ex[::500].tolist(),
                   "I_pg": [float(np.trapz(ex, tf)), float(np.trapz(b[:, 6], b[:, 0])), float(np.trapz(l[:, 6], l[:, 0]))],
                   "I_am": [float(np.trapz(b[:, 5], b[:, 0])), float(np.trapz(l[:, 5], l[:, 0]))]}
out["forces"] = fz

def rounded(v):
    if isinstance(v, float):
        return float("%.6g" % v)
    if isinstance(v, dict):
        return {k: rounded(q) for k, q in v.items()}
    if isinstance(v, (list, tuple)):
        return [rounded(q) for q in v]
    return v


json.dump(rounded(out), open(os.path.join(HERE, "data.json"), "w"), separators=(",", ":"))
print("data.json", os.path.getsize(os.path.join(HERE, "data.json")) // 1024, "kB")
