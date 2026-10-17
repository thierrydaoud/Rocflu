"""Readers for Rocflu/ppiclF output files.

- ppiclF particle files (par*.vtu): VTK XML with raw appended Float32 data.
- Rocflu ASCII flow files (SOL_<t>/<case>.mixt.cva_<reg>_<t>): double precision.
- Fluid and particle dumps (restart_data/fluid_<step>_<region>.dat,
  restart_data/particles_<step>.dat): raw stream, layout in
  modflu/RFLU_ModRestartData.F90.
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


def read_vtu_cell_field(path, name):
    """One DataArray (e.g. a CellData field) from a Rocflu ParaView region file (.vtu).

    Arrays are base64 'binary' with a UInt32 byte-count header, uncompressed.
    Multi-component arrays come back flat.
    """
    import base64
    with open(path, "rb") as f:
        text = f.read().decode("latin-1")
    m = re.search(r'<DataArray type="(\w+)"[^>]* Name="%s"[^>]*format="binary">\s*([^<]*)<'
                  % re.escape(name), text)
    if m is None:
        raise KeyError("%s not in %s" % (name, path))
    fmt = {"Float32": "f", "Float64": "d", "Int32": "i"}[m.group(1)]
    raw = base64.b64decode(m.group(2).strip())
    nbytes = struct.unpack_from("<I", raw, 0)[0]
    n = nbytes // struct.calcsize(fmt)
    return list(struct.unpack_from("<%d%s" % (n, fmt), raw, 4))


def read_vtu_cell_centroids(path):
    """Cell centroids (mean of vertices) from a Rocflu ParaView region file."""
    pts = read_vtu_cell_field(path, "Points")
    conn = read_vtu_cell_field(path, "connectivity")
    offsets = read_vtu_cell_field(path, "offsets")
    out, start = [], 0
    for end in offsets:
        ids = conn[start:end]
        out.append([sum(pts[3 * i + l] for i in ids) / len(ids) for l in range(3)])
        start = end
    return out


def latest_particle_vtu(rundir):
    files = sorted(glob.glob(os.path.join(rundir, "par[0-9]*.vtu")))
    return files[-1] if files else None


_NUM = re.compile(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?")


def particle_ids(fields, npart, nranks):
    """Global particle ids = 1-based line order in points.dat.

    id_rank is the rank that read the particle, id its index on that rank,
    and each rank reads ceil(npart/nranks) consecutive lines.
    """
    per_rank = -(-npart // nranks)
    return [int(r) * per_rank + int(i) for r, i in zip(fields["id_rank"], fields["id"])]


FLUID_DUMP_HEADER = struct.Struct("<8siqd6i")
FORM_PURE, FORM_VOLUME_WEIGHTED = 1, 2


def read_fluid_dump(path):
    """Return a dict: step, time, region, form, cv (list of nCv arrays over the
    real cells) and, for PICL runs, phi and grad_rhog (3 arrays)."""
    with open(path, "rb") as f:
        data = f.read()
    magic, version, step, time, region, ncells, ncv, form, with_picl, with_pt = \
        FLUID_DUMP_HEADER.unpack_from(data, 0)
    if magic != b"RFLUDUMP" or version != 3:
        raise ValueError("%s: not a version 3 fluid dump" % path)
    pos = FLUID_DUMP_HEADER.size

    def take():
        nonlocal pos
        a = struct.unpack_from("<%dd" % ncells, data, pos)
        pos += 8 * ncells
        return a

    out = {"step": step, "time": time, "region": region, "form": form,
           "cv": [take() for _ in range(ncv)]}
    if with_picl:
        out["phi"] = take()
        out["grad_rhog"] = [take() for _ in range(3)]
    if with_pt:
        out["div_phi_rsg"] = [take() for _ in range(3)]
        out["f_cell"] = [take() for _ in range(3)]
    # NSCBC face states: {global patch: list of nCv arrays over the faces}
    npatches, = struct.unpack_from("<i", data, pos)
    pos += 4
    out["nscbc"] = {}
    for _ in range(npatches):
        patch, nfaces = struct.unpack_from("<2i", data, pos)
        pos += 8
        out["nscbc"][patch] = []
        for _ in range(ncv):
            out["nscbc"][patch].append(struct.unpack_from("<%dd" % nfaces, data, pos))
            pos += 8 * nfaces
    if pos != len(data):
        raise ValueError("%s: %d bytes left after the data" % (path, len(data) - pos))
    return out


def fluid_dumps(rundir):
    """All fluid dumps of a run: {(step, region): path}."""
    out = {}
    for path in glob.glob(os.path.join(rundir, "restart_data", "fluid_*_*.dat")):
        _, step, region = os.path.basename(path)[:-4].split("_")
        out[(int(step), int(region))] = path
    return out


PARTICLE_DUMP_HEADER = struct.Struct("<8siqdiq7i")
OUT_Y, OUT_RPROP, OUT_TAG = 1, 2, 3


def read_particle_dump(path):
    """Return a dict: step, time, counts (per rank), sizes, history times, bins,
    the output-field table and, per particle in rank order, y, rprop, iprop and
    hist (lists of slot values)."""
    with open(path, "rb") as f:
        data = f.read()
    (magic, version, step, time, nranks, ntotal,
     lrs, lrp, lip, lrp3, lrppro, nhist, ntimebh) = PARTICLE_DUMP_HEADER.unpack_from(data, 0)
    if magic != b"PPICLDMP" or version != 3:
        raise ValueError("%s: not a version 3 particle dump" % path)
    pos = PARTICLE_DUMP_HEADER.size

    def take(fmt, n):
        nonlocal pos
        out = struct.unpack_from("<%d%s" % (n, fmt), data, pos)
        pos += struct.calcsize("<%d%s" % (n, fmt))
        return out

    time_bh = take("d", nhist)
    n_bins = take("i", 3)
    bin_box = take("d", 6)
    bin_dx = take("d", 3)
    nout, = take("i", 1)
    fields = []
    for _ in range(nout):
        source, slot = take("i", 2)
        name = data[pos:pos + 24].decode().strip()
        pos += 24
        fields.append((source, slot, name))
    nsaved, = take("i", 1)
    saved = take("i", nsaved)              # the rprop slots in the file
    counts = take("q", nranks)
    cell_counts = take("q", nranks)
    if sum(counts) != ntotal:
        raise ValueError("%s: rank counts add up to %d, not %d" % (path, sum(counts), ntotal))
    y, rprop, iprop, hist, feedback, projected = [], [], [], [], [], []
    for n, ncells in zip(counts, cell_counts):
        for arr, fmt, width in ((y, "d", lrs), (rprop, "d", nsaved), (iprop, "i", lip),
                                (hist, "d", lrp3), (feedback, "d", lrppro)):
            flat = take(fmt, width * n)
            arr.extend(flat[k * width:(k + 1) * width] for k in range(n))
        # fields projected to the rank's fluid cells: one array per field
        projected.append([take("d", ncells) for _ in range(lrppro)])
    if pos != len(data):
        raise ValueError("%s: %d bytes left after the data" % (path, len(data) - pos))
    # rprop as full rows (lrp slots, 0 where not saved)
    full = []
    for row in rprop:
        r = [0.0] * lrp
        for j, slot in enumerate(saved):
            r[slot - 1] = row[j]
        full.append(tuple(r))
    rprop = full
    return {"step": step, "time": time, "counts": counts, "lrs": lrs, "lrp": lrp, "lip": lip,
            "saved_slots": saved,
            "lrp3": lrp3, "time_bh": time_bh, "ntimebh": ntimebh, "n_bins": n_bins,
            "bin_box": bin_box, "bin_dx": bin_dx, "fields": fields,
            "y": y, "rprop": rprop, "iprop": iprop, "hist": hist, "feedback": feedback,
            "projected": projected}


def particle_dumps(rundir):
    """All particle dumps of a run: {step: path}."""
    return {int(os.path.basename(p)[10:18]): p
            for p in glob.glob(os.path.join(rundir, "restart_data", "particles_*.dat"))}


def read_post_vtu(path):
    """Return (time, step, points, fields) of a post/ .vtu (rflupost: raw appended
    data, UInt64 block sizes). fields maps name -> list of values, one tuple per
    cell or point (one value for one component)."""
    with open(path, "rb") as f:
        data = f.read()
    i = data.index(b"<AppendedData")
    header = data[:i].decode()
    base = data.index(b"_", i) + 1
    time = float(re.search(r'Name="TIME"[^>]*>\s*([-+.\dEe]+)', header).group(1))
    step = int(re.search(r'Name="STEP"[^>]*>\s*(\d+)', header).group(1))
    arrays = {}
    for m in re.finditer(r'<DataArray type="(\w+)" Name="([^"]*)" NumberOfComponents="(\d+)" '
                         r'format="appended" offset="(\d+)"/>', header):
        dtype, name, ncomp, offset = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        nbytes, = struct.unpack_from("<Q", data, base + offset)
        fmt = {"Float64": "d", "Int64": "q", "UInt8": "B"}[dtype]
        size = struct.calcsize(fmt)
        flat = struct.unpack_from("<%d%s" % (nbytes // size, fmt), data, base + offset + 8)
        arrays[name] = flat if ncomp == 1 else [flat[k:k + ncomp] for k in range(0, len(flat), ncomp)]
    points = arrays.pop("Points")
    for name in ("connectivity", "offsets", "types"):
        arrays.pop(name)
    return time, step, points, arrays


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
