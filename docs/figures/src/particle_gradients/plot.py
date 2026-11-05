#!/usr/bin/env python3
"""Draw the particle-gradient figures from data.json (written by extract.py), into docs/figures."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..")
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
BLUE, ORANGE, GREEN, RED, GREY = "#2a78d6", "#eb6834", "#2e9e5b", "#c2364b", "#9a9893"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE, "legend.frameon": False})
D = json.load(open(os.path.join(HERE, "data.json")))
mm = lambda v: 1e3 * np.asarray(v)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=130, bbox_inches="tight")
    plt.close(fig)


# 1. Riemann tube: least squares against least squares + WENO, and exact cell averages
T = D["tube"]
x = mm(T["x"]); s = 1e3 * T["shock"]
fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
for ax, key, lab, lo, hi, scale in ((axes[0], "dpdx", r"$\partial p/\partial x$ (MPa/m)", -55, 20, 1e-6),
                                    (axes[1], "dpdx", r"$\partial p/\partial x$ (MPa/m)", s - 2.2, s + 2.2, 1e-6),
                                    (axes[2], "drhodx", r"$\partial \rho_g/\partial x$ (kg/m$^4$)", s - 2.2, s + 2.2, 1.0)):
    m = (x >= lo) & (x <= hi); zoom = hi - lo < 10
    ax.plot(x[m], scale * np.asarray(T[key]["exact"])[m], color=GREY, ls="--", lw=1.1, label="exact, cell average")
    ax.plot(x[m], scale * np.asarray(T[key]["LS+WENO"])[m], color=ORANGE, lw=1.6,
            marker="s" if zoom else None, ms=4, label="least squares + WENO (before)")
    ax.plot(x[m], scale * np.asarray(T[key]["LS"])[m], color=BLUE, lw=1.6,
            marker="o" if zoom else None, ms=4, label="least squares (now)")
    ax.set_xlabel("x (mm)"); ax.set_ylabel(lab); ax.grid(color=GRID, lw=0.6)
axes[0].set_title("Rarefaction, contact and shock", fontsize=10)
axes[1].set_title("At the shock", fontsize=10); axes[2].set_title("At the shock", fontsize=10)
axes[0].legend(loc="lower left", fontsize=8)
save(fig, "fig_grad_tube.png")

# 2. The same shock in Rocflu and MFC
C = D["codes"]
fig, ax = plt.subplots(figsize=(7.2, 4.2))
style = {"exact": (GREY, "--", None), "Rocflu LS": (BLUE, "-", "o"), "Rocflu LS+WENO": (ORANGE, "-", "s"),
         "MFC WENO faces, mp_weno": (RED, "-", "^"), "MFC WENO faces, no mp_weno": (RED, ":", "v"),
         "MFC central difference": (GREEN, "-", "D")}
names = {"exact": "exact, cell average", "Rocflu LS": "Rocflu, least squares (now)",
         "Rocflu LS+WENO": "Rocflu, least squares + WENO (before)",
         "MFC WENO faces, mp_weno": "MFC particles: WENO face states, mp_weno",
         "MFC WENO faces, no mp_weno": "MFC particles: WENO face states, no mp_weno",
         "MFC central difference": "MFC gas source: central difference"}
for k, (c, ls, mk) in style.items():
    xx, yy = C[k]
    ax.plot(mm(xx), 1e-6 * np.asarray(yy), color=c, ls=ls, marker=mk, ms=4, lw=1.4, label=names[k])
ax.set_xlabel("x (mm)"); ax.set_ylabel(r"$\partial p/\partial x$ (MPa/m)"); ax.grid(color=GRID, lw=0.6)
ax.legend(fontsize=8, loc="lower left")
save(fig, "fig_grad_codes.png")

# 3. Curtain: plane averages through the bed
K = D["curtain"]; keys = sorted(K, key=int)
fig, axes = plt.subplots(3, len(keys), figsize=(4.2 * len(keys), 7.4), sharex=True,
                         gridspec_kw={"height_ratios": [1, 1.3, 1.8]})
for j, k in enumerate(keys):
    c = K[k]; xx = mm(c["x"])
    a0 = axes[0, j]; a0.fill_between(xx, 0, c["phip"], color=GRID); a0.plot(xx, c["phip"], color=MUTED, lw=1.0)
    a0.set_ylim(0, 0.24); a0.set_title("t = %.0f µs" % (1e6 * c["t"]), fontsize=10)
    a1 = axes[1, j]; a1.plot(xx, 1e-3 * np.asarray(c["p"]), color=INK, lw=1.4); a1.grid(color=GRID, lw=0.6)
    a2 = axes[2, j]
    a2.plot(xx, 1e-6 * np.asarray(c["dpdx_WENO"]), color=ORANGE, lw=1.6, label="least squares + WENO (before)")
    a2.plot(xx, 1e-6 * np.asarray(c["dpdx_LS"]), color=BLUE, lw=1.6, label="least squares (now)")
    a2.set_xlabel("x (mm)"); a2.grid(color=GRID, lw=0.6)
axes[0, 0].set_ylabel(r"$\phi_p$"); axes[1, 0].set_ylabel("p (kPa)")
axes[2, 0].set_ylabel(r"$\partial p/\partial x$ (MPa/m)"); axes[2, 0].legend(fontsize=8)
fig.align_ylabels(axes[:, 0])
save(fig, "fig_grad_curtain.png")

# 4. Forces on fixed particles
F = D["forces"]
picks = [("2.6", "Crossed by the shock, x = 2.6 mm", (38, 43)), ("-35.0", "In the rarefaction, x = -35 mm", (20, 60))]
fig, axes = plt.subplots(2, 2, figsize=(11, 6.2))
for j, (k, title, (t0, t1)) in enumerate(picks):
    key = [q for q in F if abs(float(q) * 1e3 - float(k)) < 1e-3][0]; f = F[key]
    t, te = mm(f["t"]) * 1e3, mm(f["t_exact"]) * 1e3
    a = axes[0, j]
    a.plot(te, 1e3 * np.asarray(f["pg_exact"]), color=GREY, ls="--", lw=1.1, label="exact, averaged over one cell (0.25 mm)")
    a.plot(t, 1e3 * np.asarray(f["pg_old"]), color=ORANGE, lw=1.6, label="least squares + WENO (before)")
    a.plot(t, 1e3 * np.asarray(f["pg_ls"]), color=BLUE, lw=1.6, label="least squares (now)")
    a.set_xlim(t0, t1); a.set_title(title, fontsize=10); a.set_ylabel("pressure-gradient force (mN)")
    a.grid(color=GRID, lw=0.6)
    b = axes[1, j]
    b.plot(t, 1e3 * np.asarray(f["am_old"]), color=ORANGE, lw=1.6)
    b.plot(t, 1e3 * np.asarray(f["am_ls"]), color=BLUE, lw=1.6)
    b.set_xlim(t0, t1); b.set_xlabel("t (µs)"); b.set_ylabel("added-mass force (mN)"); b.grid(color=GRID, lw=0.6)
axes[0, 0].legend(fontsize=8)
save(fig, "fig_grad_forces.png")
print("written")
