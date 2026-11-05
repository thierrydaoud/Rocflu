#!/usr/bin/env python3
"""Draw fig_shock_planar.png and fig_shock_wedge.png from results.json, into docs/figures."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..")
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})
results = json.load(open(os.path.join(HERE, "results.json")))

# Planar tube at 20 us: (label, run prefix)
ROWS = [("Aligned, Roe, order 1", "flux1_o1_aligned"),
        ("Aligned, HLLC, order 1", "flux3_o1_aligned"),
        ("Aligned, AUSM+UP, order 1", "o1_aligned"),
        ("Aligned, AUSM+, order 1", "flux4_o1_aligned"),
        ("Aligned, AUSM+UP, order 2", "aligned"),
        ("Aligned, AUSM+, order 2", "flux4_o2_aligned"),
        ("Sheared, AUSM+UP, order 1", "o1_sheared"),
        ("Sheared, AUSM+UP, order 2", "sheared"),
        ("Sheared, AUSM+, order 2", "flux4_o2_sheared")]
planar = {d["run"]: d for d in results["planar"] if d["T"] == 20.0}
panels = [("Front ahead of exact, twice-ambient (mm)", lambda d: 1e3 * (d["x2"] - d["xe"])),
          ("Transverse velocity behind the shock, rms (m/s)", lambda d: d["vrms"])]
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
for ax, (label, value) in zip(axes, panels):
    for i, (_, run) in enumerate(ROWS):
        a, b = value(planar[run + "_d2"]), value(planar[run + "_d3"])
        ax.plot([a, b], [i, i], color=GRID, lw=3, zorder=1, solid_capstyle="round")
        ax.scatter(a, i, s=52, color=BLUE, marker="o", zorder=3, edgecolors=SURFACE, linewidths=1)
        ax.scatter(b, i, s=56, color=ORANGE, marker="D", zorder=2, edgecolors=SURFACE, linewidths=1)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel(label)
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
axes[0].set_yticks(range(len(ROWS)), [r[0] for r in ROWS])
axes[0].invert_yaxis()
handles = [plt.Line2D([], [], color=BLUE, marker="o", ls="", ms=7),
           plt.Line2D([], [], color=ORANGE, marker="D", ls="", ms=7)]
fig.legend(handles, ["2D, one layer", "3D, five layers"], loc="upper center", ncol=2, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(OUT, "fig_shock_planar.png"), dpi=170)
plt.close(fig)

# Wedge at 20 us against refinement
wedge = {d["name"]: d for d in results["wedge"] if d["T"] == 20.0}
SCHEMES = [("AUSM+", BLUE, "o"), ("AUSM+UP", ORANGE, "s")]
LEVELS = (1, 2, 4)
panels = [("Front at half the peak pressure (m)", lambda d: d["fh"]),
          ("Peak pressure (MPa)", lambda d: d["pk"] * 1e-6),
          ("Tangential velocity behind the shock, rms (m/s)", lambda d: d["ut"])]
fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9))
for ax, (label, value) in zip(axes, panels):
    for scheme, color, marker in SCHEMES:
        ax.plot(LEVELS, [value(wedge["%s 2D %dx" % (scheme, k)]) for k in LEVELS], color=color,
                marker=marker, ms=7, lw=2, solid_capstyle="round", mec=SURFACE, mew=1)
        ax.scatter(1, value(wedge["%s 3D 1x" % scheme]), s=70, marker=marker, facecolors="none",
                   edgecolors=color, linewidths=1.8, zorder=3)
    ax.set_xscale("log", base=2)
    ax.set_xticks(LEVELS, ["1×", "2×", "4×"])
    ax.minorticks_off()
    ax.set_xlim(0.85, 4.6)
    ax.set_xlabel("In-plane refinement of the 2D grid")
    ax.set_title(label, fontsize=10, color=INK, loc="left")
    ax.grid(color=GRID, lw=0.7)
handles = []
names = []
for scheme, color, marker in SCHEMES:
    handles.append(plt.Line2D([], [], color=color, marker=marker, lw=2, ms=7))
    names.append(scheme + ", 2D")
    handles.append(plt.Line2D([], [], color=color, marker=marker, ls="", ms=8, mfc="none", mew=1.8))
    names.append(scheme + ", 3D at 1×")
fig.legend(handles, names, loc="upper center", ncol=4, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(os.path.join(OUT, "fig_shock_wedge.png"), dpi=170)
