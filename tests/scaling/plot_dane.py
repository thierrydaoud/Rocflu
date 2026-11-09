#!/usr/bin/env python3
"""Plot Dane fluid-only scaling from collect.py's CSV.

Usage: plot_dane.py <results.csv> <outdir>
strong_scaling.png: 1-16 full nodes, baseline 1 full node.
node_scaling.png:   1-112 ranks on one exclusive node, baseline 1 rank.
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MEASURED = "#2a6fdb"
IDEAL = "#8a8f98"
INK = "#1f2328"
MUTED = "#57606a"

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e4e7eb", "grid.linewidth": 0.6, "lines.linewidth": 2,
})


def load(path):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        for k in ("ranks", "nodes", "steps", "finished", "cells"):
            r[k] = int(r[k])
        r["t_step_s"] = float(r["t_step_s"])
    return rows


def label_points(ax, xs, ys, texts, dy=6):
    for x, y, t in zip(xs, ys, texts):
        ax.annotate(t, (x, y), textcoords="offset points", xytext=(0, dy), ha="center",
                    fontsize=8.5, color=INK)


def strong(rows, cells, out):
    runs = sorted([r for r in rows if r["study"] == "strong_fluid" and r["partition"] == "pbatch"],
                  key=lambda r: r["nodes"])
    nodes = [r["nodes"] for r in runs]
    t = [r["t_step_s"] for r in runs]
    eff = [t[0] / (ti * n / nodes[0]) for ti, n in zip(t, nodes)]
    per_rank = [cells / r["ranks"] for r in runs]

    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    a.loglog(nodes, [t[0] * nodes[0] / n for n in nodes], "--", color=IDEAL, lw=1.5, label="ideal")
    a.loglog(nodes, t, "o-", color=MEASURED, ms=7, label="measured")
    label_points(a, nodes, t, ["%.3g s" % v for v in t], dy=8)
    a.set_xlabel("nodes (112 ranks each)")
    a.set_ylabel("wall time per step [s]")
    a.set_title("Time per step", loc="left", color=INK)
    a.legend(frameon=False)
    for ax in (a, b):
        ax.set_xscale("log", base=2)
        ax.set_xticks(nodes)
        ax.set_xticklabels([str(n) for n in nodes])
        ax.minorticks_off()
        top = ax.secondary_xaxis("top")
        top.set_xticks(nodes)
        top.set_xticklabels(["%.0fk" % (c / 1e3) for c in per_rank], fontsize=8.5, color=MUTED)
        top.set_xlabel("cells per rank", color=MUTED)
        top.tick_params(colors=MUTED)
    a.set_ylim(0.08, 2.4)
    a.set_yticks([0.1, 0.25, 0.5, 1, 1.5])
    a.set_yticklabels(["0.1", "0.25", "0.5", "1", "1.5"])

    b.plot(nodes, [100] * len(nodes), "--", color=IDEAL, lw=1.5)
    b.plot(nodes, [100 * e for e in eff], "o-", color=MEASURED, ms=7)
    label_points(b, nodes, [100 * e for e in eff], ["%.0f%%" % (100 * e) for e in eff], dy=8)
    b.set_ylim(0, 110)
    b.set_xlabel("nodes (112 ranks each)")
    b.set_ylabel("parallel efficiency vs 1 node [%]")
    b.set_title("Efficiency, baseline one full node", loc="left", color=INK)
    fig.suptitle("Rocflu strong scaling on Dane: fluid only, %.1fM cells, 199 steps" % (cells / 1e6),
                 x=0.01, ha="left", color=INK, fontsize=11)
    fig.savefig(os.path.join(out, "strong_scaling.png"), dpi=200)


def on_node(rows, cells, out):
    runs = sorted([r for r in rows if r["study"] == "node_fluid" and r["nodes"] == 1],
                  key=lambda r: r["ranks"])
    spread = [r for r in rows if r["study"] == "node_fluid" and r["nodes"] == 2]
    base = runs[0]
    ranks = [r["ranks"] for r in runs]
    speed = [base["t_step_s"] / r["t_step_s"] for r in runs]
    eff = [s / n for s, n in zip(speed, ranks)]

    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    a.plot([1, 112], [1, 112], "--", color=IDEAL, lw=1.5, label="ideal")
    a.plot(ranks, speed, "o-", color=MEASURED, ms=7, label="one node")
    label_points(a, ranks[2:], speed[2:], ["%.0f" % s for s in speed[2:]], dy=-14)
    a.annotate("%.0f" % speed[1], (ranks[1], speed[1]), textcoords="offset points", xytext=(-6, 4),
               ha="right", fontsize=8.5, color=INK)
    for r in spread:
        s = base["t_step_s"] / r["t_step_s"]
        a.plot([r["ranks"]], [s], "D", ms=8, mfc="white", mec=MEASURED, mew=2,
               label="112 ranks on 2 nodes (56 each)")
        a.annotate("%.0f" % s, (r["ranks"], s), textcoords="offset points", xytext=(-12, 4),
                   ha="right", fontsize=8.5, color=INK)
        b.plot([r["ranks"]], [100 * s / r["ranks"]], "D", ms=8, mfc="white", mec=MEASURED, mew=2)
        b.annotate("%.0f%%" % (100 * s / r["ranks"]), (r["ranks"], 100 * s / r["ranks"]),
                   textcoords="offset points", xytext=(-10, 4), ha="right", fontsize=8.5, color=INK)
    a.set_xlabel("MPI ranks")
    a.set_ylabel("speed-up vs 1 rank")
    a.set_title("Speed-up", loc="left", color=INK)
    a.set_xlim(0, 118)
    a.set_ylim(0, 118)
    a.set_xticks([1, 8, 16, 28, 56, 112])
    a.legend(frameon=False, loc="upper left")

    b.plot([1, 112], [100, 100], "--", color=IDEAL, lw=1.5)
    b.plot(ranks, [100 * e for e in eff], "o-", color=MEASURED, ms=7)
    label_points(b, ranks[1:], [100 * e for e in eff[1:]], ["%.0f%%" % (100 * e) for e in eff[1:]], dy=8)
    b.set_xlim(0, 118)
    b.set_ylim(0, 110)
    b.set_xticks([1, 8, 16, 28, 56, 112])
    b.set_xlabel("MPI ranks")
    b.set_ylabel("efficiency vs 1 rank [%]")
    b.set_title("Efficiency", loc="left", color=INK)
    note = "" if base["finished"] else "  (1-rank time from step %d of 199)" % base["steps"]
    fig.suptitle("Rocflu on one Dane node (2 x 56-core Sapphire Rapids): fluid only, %.1fM cells%s"
                 % (cells / 1e6, note), x=0.01, ha="left", color=INK, fontsize=11)
    fig.savefig(os.path.join(out, "node_scaling.png"), dpi=200)


def main():
    rows = load(sys.argv[1])
    out = sys.argv[2]
    cells = rows[0]["cells"]
    strong(rows, cells, out)
    on_node(rows, cells, out)


if __name__ == "__main__":
    main()
