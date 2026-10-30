"""Test case definition, set-up, and execution.

A case is a box grid from tests/mesh_generator, the base input deck in
tests/templates with per-case overrides, a boundary-condition layout, and a
list of particles. Setting up runs rflumap/rflupart/rfluinit serially;
the solver runs through srun.
"""

import filecmp
import hashlib
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

TESTS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATES = os.path.join(TESTS_DIR, "templates")
CASENAME = "rectshktb"  # Rocflu requires files named after the case


@dataclass
class Grid:
    """Box from (x0,y0,z0) with lengths (lx,ly,lz) and nx*ny*nz cells."""
    lx: float = 0.04
    ly: float = 0.004
    lz: float = 0.004
    x0: float = -0.02
    y0: float = 0.0
    z0: float = 0.0
    nx: int = 40
    ny: int = 4
    nz: int = 4

    @property
    def cell_volume(self):
        return self.lx * self.ly * self.lz / (self.nx * self.ny * self.nz)


@dataclass
class Case:
    name: str
    description: str
    nranks: int = 1
    grid: Grid = field(default_factory=Grid)
    # "shocktube": inflow x-, outflow x+, periodic y,z (templates/rectshktb.bc)
    # "periodic":  periodic in x, y and z
    # "closed":    slip walls on all six faces
    # "partwall":  slip walls on all six faces, particle walls too (PARTWALL 1)
    # "symmetry":  symmetry patches on all six faces
    bc: str = "shocktube"
    # {("SECTION", "KEY"): value} applied to templates/rectshktb.inp
    inp: Dict[Tuple[str, str], object] = field(default_factory=dict)
    # (x, y, z, dp) or (x, y, z, dp, u, v, w); material is always AL.
    # Particles without a velocity start at rest.
    particles: List[Tuple[float, ...]] = field(default_factory=list)
    # Physics checks run on every result (see checks.py)
    checks: List[str] = field(default_factory=list)
    # Name of another case whose final state must match this one
    same_as: Optional[str] = None
    # Particle fields compared with same_as (name prefixes, e.g. ("vp_",) for
    # the particle velocity only); None compares every field
    same_as_fields: Optional[Tuple[str, ...]] = None
    # Known code problem that makes the physics/cross checks fail. The case
    # reports XFAIL while they fail (golden still guards against changes) and
    # XPASS once they pass, so the marker can be removed.
    known_issue: Optional[str] = None
    # After the run: restart it from its first dump after step 0 and require
    # the final dumps to be byte-identical (restart_test below)
    restart_test: bool = False
    timeout_min: int = 15

    @property
    def uuid(self):
        """Short id that changes whenever the case definition changes."""
        spec = json.dumps(self.spec(), sort_keys=True, default=str)
        return hashlib.sha1(spec.encode()).hexdigest()[:8].upper()

    def spec(self):
        return {
            "nranks": self.nranks,
            "grid": self.grid.__dict__,
            "bc": self.bc,
            "inp": {"%s/%s" % k: v for k, v in self.inp.items()},
            "particles": self.particles,
        }


# ---------------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------------

def edit_inp(text, overrides):
    """Replace KEY values inside '# SECTION' blocks of a Rocflu input deck."""
    lines = text.split("\n")
    for (section, key), value in overrides.items():
        in_section = False
        found = False
        for i, line in enumerate(lines):
            s = line.strip()
            if s.startswith("#"):
                in_section = s[1:].strip().split() == [section]
                continue
            if in_section and s.split()[:1] == [key]:
                lines[i] = "%-16s %s" % (key, value)
                found = True
                break
        if not found:
            raise KeyError("No %s in section # %s of the input deck" % (key, section))
    return "\n".join(lines)


def _periodic_block(patch, rel, axis, name):
    return ("# BC_PERIODIC\nPATCH    %d %d\nRELPATCH %d %d\nAXIS     %d\n"
            "ANGLE    0.0\nNAME     %s\n#\n\n" % (patch, patch, rel, rel, axis, name))


def bc_text(kind):
    """Boundary-condition file. Patches 1-6 are x-, x+, y-, y+, z-, z+."""
    if kind == "shocktube":
        with open(os.path.join(TEMPLATES, CASENAME + ".bc")) as f:
            return f.read()
    if kind == "periodic":
        text = ""
        for axis, (lo, hi) in enumerate([(1, 2), (3, 4), (5, 6)], start=1):
            name = "XYZ"[axis - 1] + "-Periodic"
            text += _periodic_block(lo, hi, axis, name) + _periodic_block(hi, lo, axis, name)
        return text + "# END\n"
    if kind == "closed":
        return "# BC_SLIPW\nPATCH      1 6\nNAME       SlipWalls\n#\n\n# END\n"
    if kind == "partwall":
        return ("# BC_SLIPW\nPATCH      1 6\nNAME       SlipWalls\nPARTWALL   1\n#\n\n"
                "# END\n")
    if kind == "symmetry":
        return "# BC_SYMMETRY\nPATCH      1 6\nNAME       Symmetry\n#\n\n# END\n"
    raise ValueError("Unknown bc layout: " + kind)


def points_text(particles):
    """points.dat: x y z vx vy vz dp per line. Particles are (x, y, z, dp) or
    (x, y, z, dp, u, v, w) in the case definitions; without u v w they start
    at rest."""
    lines = ["! generated by tests/rocpicl_tests", str(len(particles))]
    for p in particles:
        velocity = tuple(p[4:7]) if len(p) > 4 else (0.0, 0.0, 0.0)
        lines.append("  ".join("%.10e" % v for v in tuple(p[:3]) + velocity + (p[3],)))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Set-up and run
# ---------------------------------------------------------------------------

class RunError(RuntimeError):
    pass


def _run(cmd, cwd, log, stdin=None, timeout=None):
    with open(os.path.join(cwd, log), "w") as f:
        p = subprocess.run(cmd, cwd=cwd, input=stdin, stdout=f, stderr=subprocess.STDOUT,
                           text=True, timeout=timeout)
    return p.returncode


def setup(case, rundir, exe_dir):
    """Create rundir with grid, inputs, and pre-processed region files."""
    if os.path.exists(rundir):
        shutil.rmtree(rundir)
    os.makedirs(rundir)

    g = case.grid
    answers = "\n".join(str(v) for v in (
        g.lx, g.ly, g.lz, g.x0, g.y0, g.z0, g.nx + 1, g.ny + 1, g.nz + 1, 2)) + "\n"
    if _run([os.path.join(TESTS_DIR, "mesh_generator")], rundir, "mesh.log", stdin=answers):
        raise RunError("mesh_generator failed (see mesh.log)")
    os.rename(os.path.join(rundir, "hass.hyb.asc"), os.path.join(rundir, CASENAME + ".hyb.asc"))

    with open(os.path.join(TEMPLATES, CASENAME + ".inp")) as f:
        inp = edit_inp(f.read(), case.inp)
    with open(os.path.join(rundir, CASENAME + ".inp"), "w") as f:
        f.write(inp)
    with open(os.path.join(rundir, CASENAME + ".bc"), "w") as f:
        f.write(bc_text(case.bc))
    with open(os.path.join(rundir, "points.dat"), "w") as f:
        f.write(points_text(case.particles))

    n = str(case.nranks)
    steps = [
        ("rflumap", ["-c", CASENAME, "-m", "1", "-p", n, "-r", n, "-v", "1"]),
        ("rflupart", ["-c", CASENAME, "-v", "1"]),
        ("rfluinit", ["-c", CASENAME, "-v", "1"]),
    ]
    for exe, args in steps:
        rc = _run([os.path.join(exe_dir, exe)] + args, rundir, exe + ".log")
        # Rocflu utilities can stop on ErrorStop with exit code 0
        with open(os.path.join(rundir, exe + ".log"), errors="replace") as f:
            failed = any(line.lstrip().startswith("ERROR") for line in f)
        if rc or failed:
            raise RunError("%s failed (see %s.log)" % (exe, exe))


def srun_prefix(nranks, timeout_min, slurm):
    """srun inside an existing allocation, or as its own job otherwise.
    Small cases (nranks <= slurm["local_ranks"]) run with mpirun on this node."""
    if nranks <= slurm.get("local_ranks", 0):
        return ["mpirun", "-np", str(nranks), "--bind-to", "none"]
    cmd = ["srun", "-n", str(nranks), "--mpi=pmix"]
    if not os.environ.get("SLURM_JOB_ID"):
        cmd += ["-A", slurm["account"], "-p", slurm["partition"],
                "-t", "00:%02d:00" % timeout_min, "--mem-per-cpu=2gb"]
    return cmd


def run_solver(case, rundir, exe_dir, slurm):
    cmd = srun_prefix(case.nranks, case.timeout_min, slurm)
    cmd += [os.path.join(exe_dir, "rflump"), "-c", CASENAME, "-v", "1"]
    rc = _run(cmd, rundir, "rflump.log", timeout=60 * (case.timeout_min + 5))
    with open(os.path.join(rundir, "rflump.log"), errors="replace") as f:
        log = f.read()
    if rc != 0 or "| Finished" not in log:
        # srun on a shared partition can fail for reasons unrelated to the
        # code; report its messages so such failures are recognisable
        slurm_lines = [l.strip() for l in log.splitlines()
                       if "srun:" in l or "slurm" in l.lower()][:3]
        raise RunError("rflump failed (rc=%d; see rflump.log)%s"
                       % (rc, ("; slurm: " + " | ".join(slurm_lines)) if slurm_lines else ""))


def restart_test(case, rundir, exe_dir, slurm):
    """Restart the finished run in a copy of its directory (rundir + "_restart")
    from its first dump after step 0, to the same end time, and require every
    final dump (fluid of each region, particles) to be byte-identical to the
    uninterrupted run's. Returns a list of failure messages."""
    data = os.path.join(rundir, "restart_data")
    steps = sorted({int(n[6:14]) for n in os.listdir(data)
                    if n.startswith("fluid_") and n.endswith(".dat")})
    if len(steps) < 3:
        return ["restart test: needs dumps at step 0, a middle step and the end "
                "(found steps %s)" % steps]
    start, last = steps[1], steps[-1]

    again = rundir + "_restart"
    shutil.rmtree(again, ignore_errors=True)
    shutil.copytree(rundir, again)
    inp = os.path.join(again, CASENAME + ".inp")
    with open(inp) as f:
        text = f.read()
    text, n = re.subn(r"(?m)^START_STEP\b.*$", "START_STEP  %d" % start, text)
    if n != 1:
        return ["restart test: no START_STEP line in the input"]
    with open(inp, "w") as f:
        f.write(text)

    cmd = srun_prefix(case.nranks, case.timeout_min, slurm)
    cmd += [os.path.join(exe_dir, "rflump"), "-c", CASENAME, "-v", "1"]
    rc = _run(cmd, again, "rflump.log", timeout=60 * (case.timeout_min + 5))
    with open(os.path.join(again, "rflump.log"), errors="replace") as f:
        log = f.read()
    if rc != 0 or "| Finished" not in log:
        return ["restart test: rflump failed restarting from step %d (rc=%d; see %s)"
                % (start, rc, os.path.join(again, "rflump.log"))]

    final = sorted(n for n in os.listdir(data) if n.endswith(".dat") and
                   (n.startswith("fluid_%08d_" % last) or n == "particles_%08d.dat" % last))
    errors = []
    for name in final:
        restarted = os.path.join(again, "restart_data", name)
        if not os.path.exists(restarted):
            errors.append("restart test: restarted from step %d, no %s" % (start, name))
        elif not filecmp.cmp(os.path.join(data, name), restarted, shallow=False):
            errors.append("restart test: restarted from step %d, %s differs from the "
                          "uninterrupted run" % (start, name))
    return errors


def run_post(rundir, exe_dir):
    """Convert every dump step of the run into post/ with rflupost (one process).
    Returns a list of failure messages."""
    rc = _run([os.path.join(exe_dir, "rflupost"), "-c", CASENAME, "-s", "all", "-v", "0"],
              rundir, "rflupost.log")
    with open(os.path.join(rundir, "rflupost.log"), errors="replace") as f:
        log = f.read()
    if rc != 0 or "collection.pvd and fluid.pvd list" not in log:
        return ["rflupost failed (rc=%d; see %s)" % (rc, os.path.join(rundir, "rflupost.log"))]
    return []
