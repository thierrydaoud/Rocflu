#!/usr/bin/python3
"""Radial profiles of the 3D and 2D wedge runs. Usage: plot_wedge.py <dir with profiles*.json>"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = sys.argv[1]
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
SERIES = [
    ("3D", "3D, 1 mm, 40,781 particles", "#2a78d6", "-", 1.8),
    ("2D100s1", "2D, 100 µm, bed 1", "#eb6834", "-", 1.3),
    ("2D100s2", "2D, 100 µm, bed 2", "#1baf7a", "--", 1.3),
    ("2D125s1", "2D, 125 µm, bed 1", "#eda100", ":", 1.6),
    ("2Dcfl37", "2D, 100 µm, bed 1, CFL 0.37", "#e87ba4", "-.", 1.5),
]
TIMES = ("5", "10", "20")
FIELDS = [("rho", "Gas density (kg/m³)", 1.0, False),
          ("p", "Pressure (MPa)", 1.0e-6, True),
          ("ur", "Radial velocity (km/s)", 1.0e-3, False)]
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

data = json.load(open(os.path.join(HERE, "profiles.json")))
extra = os.path.join(HERE, "profiles_cfl.json")
if os.path.exists(extra):
    data.update(json.load(open(extra)))
runs = [s for s in SERIES if s[0] in data]


def front(run, t):
    """Largest radius where the pressure is above twice ambient."""
    d = data[run][t]
    return max(r for r, p in zip(d["r"], d["p"]) if p > 2.0 * 101325.0)


def lines(ax, t, field, scale):
    for key, _, color, style, width in runs:
        d = data[key][t]
        ax.plot(d["r"], [v * scale for v in d[field]], color=color, ls=style, lw=width,
                solid_capstyle="round")


def finish(fig, ax_for_legend, name, top=0.90):
    handles = [plt.Line2D([], [], color=c, ls=s, lw=w) for _, _, c, s, w in runs]
    fig.legend(handles, [lab for _, lab, _, _, _ in runs], loc="upper center",
               ncol=len(runs), frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, top))
    fig.savefig(os.path.join(HERE, name), dpi=170)
    plt.close(fig)


# Figure 1: whole domain
fig, axes = plt.subplots(3, 3, figsize=(13, 9), sharex="col")
for j, t in enumerate(TIMES):
    axes[0, j].set_title("t = %s µs" % t, color=INK, fontsize=11)
    for i, (field, label, scale, log) in enumerate(FIELDS):
        ax = axes[i, j]
        lines(ax, t, field, scale)
        ax.grid(color=GRID, lw=0.7)
        if log:
            ax.set_yscale("log")
        if j == 0:
            ax.set_ylabel(label)
        if i == 2:
            ax.set_xlabel("Radius r (m)")
        ax.set_xlim(0.10, 0.56)
finish(fig, axes, "wedge_profiles_all.png", top=0.95)

# Figure 2: leading shock, same width in every panel
HALF = 0.010
fig, axes = plt.subplots(2, 3, figsize=(13, 6.6), sharex="col")
for j, t in enumerate(TIMES):
    centre = front("3D", t)
    axes[0, j].set_title("t = %s µs, leading shock" % t, color=INK, fontsize=11)
    for i, (field, label, scale, log) in enumerate((FIELDS[1], FIELDS[2])):
        ax = axes[i, j]
        lines(ax, t, field, scale)
        ax.set_xlim(centre - 1.5 * HALF, centre + 0.5 * HALF)
        ax.grid(color=GRID, lw=0.7)
        window = [v * scale for k in [r[0] for r in runs]
                  for r_, v in zip(data[k][t]["r"], data[k][t][field])
                  if centre - 1.5 * HALF <= r_ <= centre + 0.5 * HALF]
        ax.set_ylim(min(0.0, min(window)), 1.08 * max(window))
        if j == 0:
            ax.set_ylabel(label)
        if i == 1:
            ax.set_xlabel("Radius r (m)")
finish(fig, axes, "wedge_leading_shock.png", top=0.93)

# Figure 3: where the particles are
fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
for j, t in enumerate(TIMES):
    ax = axes[j]
    all_r = data["3D"][t]["pr"]
    lo, hi = min(all_r) - 0.002, max(all_r) + 0.002
    nbin = 40
    width = (hi - lo) / nbin
    for key, _, color, style, lw in runs:
        counts = [0] * nbin
        rs = data[key][t]["pr"]
        for r in rs:
            k = int((r - lo) / width)
            if 0 <= k < nbin:
                counts[k] += 1
        ax.plot([lo + (k + 0.5) * width for k in range(nbin)],
                [100.0 * c / len(rs) for c in counts], color=color, ls=style, lw=lw)
    ax.set_title("t = %s µs" % t, color=INK, fontsize=11)
    ax.set_xlabel("Radius r (m)")
    ax.grid(color=GRID, lw=0.7)
    if j == 0:
        ax.set_ylabel("Share of particles per bin (%)")
finish(fig, axes, "wedge_particles.png", top=0.88)
print("written to", HERE)
