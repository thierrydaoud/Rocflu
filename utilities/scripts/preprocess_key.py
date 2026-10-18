#!/usr/bin/env python3
"""Keys of a case's pre-processing, to tell whether it must run again.

    preprocess_key.py <dir>/<case>.inp <regions>

Prints two lines:

    partition <key>    what rflumap and rflupart depend on: the number of
                       regions, the grid, the .bc file and the input file
    init <key>         what rfluinit depends on: the partition and the
                       input file

The input file counts without comments, blank lines and the settings that
cannot change these outputs (IGNORED and PARTITION_IGNORED below). The
executables do not count: ./rocflu.sh run --repartition redoes everything.
Without a grid the script knows (CENTAUR ASCII, GRIDSRC 0) it prints no
keys, so the pre-processing always runs. ./rocflu.sh run keeps the keys of
the last pre-processing that succeeded in <dir>/.rocflu_pre.stamp.
"""

import hashlib
import os
import sys

# Sections, or (section, key), that change neither the partition nor the
# initial fields: run length, output and probes
IGNORED = {
    "POST", "PROBE",
    ("NUMERICS", "CFL"),
    ("TIMESTEP", "MAXTIME"), ("TIMESTEP", "RUNTIME"), ("TIMESTEP", "WRITIME"),
    ("TIMESTEP", "PRNTIME"), ("TIMESTEP", "FOMWRITIME"), ("TIMESTEP", "TIMESTEP"),
    ("TIMESTEP", "DTMINLIMIT"), ("TIMESTEP", "WALLTIME"), ("TIMESTEP", "SAFEWRTIME"),
    ("TIMESTEP", "START_STEP"),
}

# ... and those only the initial fields depend on (rflupart does not read them)
PARTITION_IGNORED = IGNORED | {
    "INITFLOW", "PICL", "REFERENCE", "MIXTURE", ("TIMESTEP", "STARTTIME"),
}


def input_lines(path):
    """[(section, key, text)] of the input file: comments, blank lines and
    the section end lines '#' left out, spaces made single."""
    lines = []
    section = None
    with open(path) as f:
        for line in f:
            text = " ".join(line.split("!")[0].split())
            if not text:
                continue
            if text.startswith("#"):
                name = text[1:].strip()
                section = name or None
                if name:
                    lines.append((name, None, text))
                continue
            lines.append((section, text.split()[0], text))
    return lines


def input_key(lines, ignored):
    kept = [text for section, key, text in lines
            if section not in ignored and (section, key) not in ignored]
    return "\n".join(kept)


def file_digest(path):
    digest = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def keys(inp, regions):
    """(partition key, init key), or None when the grid is not known."""
    base = os.path.splitext(inp)[0]
    lines = input_lines(inp)
    source = [text.split()[1:2] for section, key, text in lines
              if section == "FORMATS" and key == "GRIDSRC"]
    grid = base + ".hyb.asc"
    if (source and source[0] != ["0"]) or not os.path.exists(grid):
        return None
    with open(base + ".bc", "rb") as f:
        bc = f.read()
    partition = hashlib.sha1("\n".join([
        "regions %d" % regions, "grid " + file_digest(grid),
        "bc " + hashlib.sha1(bc).hexdigest(),
        input_key(lines, PARTITION_IGNORED)]).encode()).hexdigest()
    init = hashlib.sha1("\n".join([
        "partition " + partition, input_key(lines, IGNORED)]).encode()).hexdigest()
    return partition, init


def main(argv):
    if len(argv) != 3 or argv[1] in ("-h", "--help"):
        print(__doc__.strip())
        return 0 if len(argv) == 2 else 2
    found = keys(argv[1], int(argv[2]))
    if found:
        print("partition %s\ninit %s" % found)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
