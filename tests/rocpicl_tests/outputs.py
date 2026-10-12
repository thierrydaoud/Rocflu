"""Readers for Rocflu/ppiclF output files.

- ppiclF particle files (par*.vtu): VTK XML with raw appended Float32 data.
- Rocflu ASCII flow files (SOL_<t>/<case>.mixt.cva_<reg>_<t>): double precision.
"""

import glob
import os
import re
import struct


def read_particle_vtu(path):
    """Return (time, fields) where fields maps name -> list of floats per particle.

    Position is split into x, y, z. Names are stripped of padding.
    """
    with open(path, "rb") as f:
        data = f.read()
    i = data.find(b"<AppendedData")
    header = data[:i].decode()
    # Raw payload starts after the '_' marker
    payload = data[data.index(b"_", i) + 1:]

    time = float(re.search(r'Name="TIME"[^>]*>\s*([-+.\dEe]+)', header).group(1))
    npts = int(re.search(r'NumberOfPoints="(\d+)"', header).group(1))

    fields = {}
    for m in re.finditer(r"<DataArray ([^>]*)/>", header):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        if attrs.get("format") != "append":
            continue
        name = attrs["Name"].strip()
        ncomp = int(attrs.get("NumberOfComponents", "1"))
        off = int(attrs["offset"])
        fmt = {"Float32": "f", "Int32": "i"}[attrs["type"]]
        nbytes = struct.unpack_from("<i", payload, off)[0]
        n = nbytes // 4
        vals = struct.unpack_from("<%d%s" % (n, fmt), payload, off + 4)
        if name == "Position":
            for c, axis in enumerate("xyz"):
                fields[axis] = list(vals[c::3])
        elif ncomp == 1 and n == npts:
            fields[name] = list(vals)
    return time, fields


def latest_particle_vtu(rundir):
    files = sorted(glob.glob(os.path.join(rundir, "par[0-9]*.vtu")))
    return files[-1] if files else None


_NUM = re.compile(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?")


def particle_ids(fields, npart, nranks):
    """Global particle ids = 1-based line order in points.dat.

    tag02 is the rank that read the particle, tag01 its index on that rank,
    and each rank reads ceil(npart/nranks) consecutive lines.
    """
    per_rank = -(-npart // nranks)
    return [int(r) * per_rank + int(i) for r, i in zip(fields["tag02"], fields["tag01"])]


def read_flow_ascii(path):
    """Return (time, list of variable arrays) from a Rocflu ASCII flow file."""
    with open(path) as f:
        lines = f.read().split("\n")
    time = None
    arrays = []
    cur = None
    k = 0
    while k < len(lines):
        line = lines[k].strip()
        if line.startswith("# Physical time"):
            time = float(lines[k + 1])
            k += 2
            continue
        if re.match(r"^#\s+\d+$", line):
            cur = []
            arrays.append(cur)
        elif line.startswith("#"):
            cur = None
        elif cur is not None and line:
            # Fixed-width E23.16: negative values abut the previous one
            cur.extend(float(v) for v in _NUM.findall(line))
        k += 1
    return time, arrays


def latest_flow(rundir, casename, kind="cva"):
    """Return {region: (time, arrays)} for the latest SOL_* directory."""
    sols = glob.glob(os.path.join(rundir, "SOL_*"))
    if not sols:
        return None
    latest = max(sols, key=lambda p: float(os.path.basename(p)[4:]))
    out = {}
    for path in sorted(glob.glob(os.path.join(latest, "%s.mixt.%s_*" % (casename, kind)))):
        region = int(os.path.basename(path).split("_")[1])
        out[region] = read_flow_ascii(path)
    return out
