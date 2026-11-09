#!/usr/bin/env python3
"""Collect scaling results into a CSV.

Usage: collect.py <study dir> [<study dir> ...] > results.csv
A study dir holds <grid>/r<ranks>[_<tag>]/ run directories; each
has the solver log rflump.<jobid>.log and the prep files <step>.time.
Nodes and partition come from sacct. A run that has not finished is listed
with finished = 0 and its time per step so far.
"""

import csv
import glob
import os
import re
import subprocess
import sys

PROGRESS = re.compile(r"^\s*\[\s*(\d+)%\].*@ step\s+(\d+)(?:.*avg ([0-9.E+-]+) s\).*total ([0-9:]+)\))?")
MEMORY = re.compile(r"Peak memory:\s+([0-9.]+) (MB|GB)/rank \(max\), ([0-9.]+) (MB|GB) total")
TOTAL = re.compile(r"Total-time:\s+([0-9.]+)s")
PREP_STEPS = ("mesh", "rflupart", "rfluinit")


def seconds(hms):
    """Seconds from h:mm:ss, m:ss or m:ss.ss."""
    parts = [float(p) for p in hms.split(":")]
    total = 0.0
    for p in parts:
        total = 60 * total + p
    return total


def gigabytes(value, unit):
    return float(value) / (1024.0 if unit == "MB" else 1.0)


def cell_count(rundir):
    """Cells of the box from mesh.in (vertex counts on lines 7-9)."""
    lines = open(os.path.join(rundir, "mesh.in")).read().split()
    return (int(lines[6]) - 1) * (int(lines[7]) - 1) * (int(lines[8]) - 1)


def job_info(jobid):
    """Nodes and partition of a Slurm job (empty if sacct does not know it)."""
    try:
        out = subprocess.check_output(["sacct", "-X", "-n", "-P", "-j", jobid, "-o", "NNodes,Partition"],
                                      universal_newlines=True).strip().splitlines()
        nodes, partition = out[0].split("|")
        return int(nodes), partition
    except (OSError, subprocess.CalledProcessError, IndexError, ValueError):
        return "", ""


def prep_row(rundir):
    row = {}
    for step in PREP_STEPS:
        path = os.path.join(rundir, step + ".time")
        wall = mem = ""
        if os.path.exists(path):
            for line in open(path):
                if "Elapsed (wall" in line:
                    wall = seconds(line.split()[-1])
                elif "Maximum resident" in line:
                    mem = int(line.split()[-1]) / 1048576.0
        row[step + "_s"] = wall
        row[step + "_GB"] = mem
    return row


def run_rows(study):
    for griddir in sorted(glob.glob(os.path.join(study, "*", ""))):
        for rundir in sorted(glob.glob(os.path.join(griddir, "r[0-9]*"))):
            name = os.path.basename(rundir)
            ranks = int(name[1:5])
            for log in sorted(glob.glob(os.path.join(rundir, "rflump.*.log"))):
                jobid = log.split(".")[-2]
                nodes, partition = job_info(jobid)
                row = {"study": os.path.basename(os.path.normpath(study)),
                       "grid": os.path.basename(os.path.normpath(griddir)), "run": name,
                       "cells": cell_count(rundir), "ranks": ranks, "nodes": nodes, "partition": partition,
                       "jobid": jobid, "steps": "", "finished": 0, "t_step_s": "", "loop_s": "",
                       "total_s": "", "mem_rank_GB": "", "mem_total_GB": ""}
                for line in open(log, errors="replace"):
                    m = PROGRESS.match(line)
                    if m:
                        row["steps"] = int(m.group(2))
                        if m.group(3):
                            row["t_step_s"] = float(m.group(3))
                            row["loop_s"] = seconds(m.group(4))
                        row["finished"] = int(m.group(1) == "100")
                    m = TOTAL.search(line)
                    if m:
                        row["total_s"] = float(m.group(1))
                    m = MEMORY.search(line)
                    if m:
                        row["mem_rank_GB"] = gigabytes(m.group(1), m.group(2))
                        row["mem_total_GB"] = gigabytes(m.group(3), m.group(4))
                row.update(prep_row(rundir))
                yield row


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    fields = ["study", "grid", "run", "cells", "ranks", "nodes", "partition", "jobid", "steps", "finished",
              "t_step_s", "loop_s", "total_s", "mem_rank_GB", "mem_total_GB"]
    fields += [s + suffix for s in PREP_STEPS for suffix in ("_s", "_GB")]
    out = csv.DictWriter(sys.stdout, fieldnames=fields)
    out.writeheader()
    for study in sys.argv[1:]:
        for row in run_rows(study):
            out.writerow({k: ("%.4g" % v if isinstance(v, float) else v) for k, v in row.items()})


if __name__ == "__main__":
    main()
