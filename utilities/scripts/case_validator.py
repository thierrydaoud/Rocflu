#!/usr/bin/env python3
"""Checks a Rocflu + ppiclF case before it runs.

    case_validator.py <dir>/<case>.inp

Reads the input file the way Rocflu does and reports every problem it finds
at once: unknown, removed or misplaced keys in the # PICL section, values out
of range, combinations of particle models that do not work together, and
files in the case directory that would stop the run. Exit status 0 when the
case is fine, 1 when it is not.

./rocflu.sh run calls it before building, ./rocflu.sh validate on its own,
and tests/run_tests.py on every case. The allowed values and defaults below
mirror rocpicl/PICL_ReadPiclSection.F90 and libflu/RFLU_InitGlobal.F90: keep
them in step when either changes.
"""

import difflib
import os
import sys

# Keys of the # PICL section, in the order PICL_ReadPiclSection lists them.
# Rocflu matches a line to the first key not yet read that the line starts
# with, so the order decides how QSFLUCT/QSFLUCTFILT/QSFLUCTFILTADAPT lines
# are read.
PICL_KEYS = [
    "QSFORCE", "AMFORCE", "PGFORCE", "COLLISIONS", "HEATTRANSFER", "FEEDBACK",
    "VELOCITYMODEL", "USED", "QSFLUCT", "TEMP", "DU HEATTRANSFER",
    "NEIGHBORWIDTH", "KSP", "EREST", "QSFLUCTFILT", "QSFLUCTFILTADAPT",
    "VISCOUSUNSTEADY", "PERIODICX", "PERIODICY", "PERIODICZ", "SUBBIN",
    "PSEUDOTURB", "VUSTEPS", "OUTPARTICLE", "OUTGAS", "OUTGRADGAS", "OUTMODEL",
    "DENSITY", "SPECIFICHEAT",
]

# Values used when a key is not in the file (RFLU_InitGlobal). KSP and EREST
# have none: they must be given when collisions are on.
PICL_DEFAULTS = {
    "QSFORCE": 0, "AMFORCE": 0, "PGFORCE": 0, "COLLISIONS": 0,
    "HEATTRANSFER": 0, "FEEDBACK": 0, "VELOCITYMODEL": 1, "USED": 0,
    "QSFLUCT": 0, "TEMP": 0.0, "DU HEATTRANSFER": 0, "NEIGHBORWIDTH": 0.0,
    "QSFLUCTFILT": 0, "QSFLUCTFILTADAPT": 0, "VISCOUSUNSTEADY": 0,
    "PERIODICX": 0, "PERIODICY": 0, "PERIODICZ": 0, "SUBBIN": 0,
    "PSEUDOTURB": 0, "VUSTEPS": 50, "OUTPARTICLE": 0, "OUTGAS": 0,
    "OUTGRADGAS": 0, "OUTMODEL": 0, "DENSITY": 0.0, "SPECIFICHEAT": 0.0,
}

# Whole-number flags and their allowed range
PICL_FLAGS = {
    "USED": (0, 1),
    "VELOCITYMODEL": (0, 1),        # 0 particles keep their state, 1 they move
    "QSFORCE": (0, 3),              # none, Gidaspow, Parmar, Osnes
    "QSFLUCT": (0, 2),              # none, Lattanzi, Osnes
    "AMFORCE": (0, 2),              # none, Parmar, Briney
    "PGFORCE": (0, 1),
    "COLLISIONS": (0, 4),
    "VISCOUSUNSTEADY": (0, 4),      # none, trapezoidal (uniform, variable dt),
                                    # modified trapezoidal, Hinsberg
    "HEATTRANSFER": (0, 5),         # none, Stokes, Ranz-Marshall, Gunn, Fox, Sun
    "DU HEATTRANSFER": (0, 1),
    "FEEDBACK": (0, 1),
    "PSEUDOTURB": (0, 1),
    "SUBBIN": (0, 1),
    "PERIODICX": (0, 1), "PERIODICY": (0, 1), "PERIODICZ": (0, 1),
    "OUTPARTICLE": (0, 1), "OUTGAS": (0, 1), "OUTGRADGAS": (0, 1),
    "OUTMODEL": (0, 1),
}

# Keys only read when a model is on: left in the file with the model off,
# they are ignored, so they are only checked when it is on.
DEPENDENT_FLAGS = {
    "QSFLUCTFILT": (0, 1), "QSFLUCTFILTADAPT": (0, 1),
}

# Keys that were removed, with what to do instead
REMOVED_KEYS = {
    "STATIONARY": "replaced by VELOCITYMODEL (STATIONARY 1 is VELOCITYMODEL 0, "
                  "STATIONARY 0 is VELOCITYMODEL 1)",
    "FILTERWIDTH": "no longer read: the filter width is set from the grid "
                   "(2 x the largest cell length around the particles); remove the line",
    "REACTIVE": "the particle burn-rate model was removed; remove the line",
    "PICLDEBUG": "removed; remove the line",
    "ANGULARPERIODIC": "angular periodicity is not available; remove the line",
    "ANGLE": "angular periodicity is not available; remove the line",
    "XANGLE": "angular periodicity is not available; remove the line",
    "RIN": "angular periodicity is not available; remove the line",
    "ROUT": "angular periodicity is not available; remove the line",
}

RK3 = 2   # RKSCHEME of the 3-stage Runge-Kutta scheme (RK_SCHEME_3_WRAY)
RK4 = 1   # Rocflu's default (RK_SCHEME_4_CLASSICAL)


def read_number(text):
    """First value of a key line, as Fortran's list-directed READ takes it."""
    token = text.replace(",", " ").split("!")[0].split()
    if not token:
        return None
    try:
        return float(token[0].replace("d", "e").replace("D", "e"))
    except ValueError:
        return None


def read_sections(lines):
    """Sections of the input file: {name: [(line number, text), ...]}.

    A section starts at a line '# NAME' and ends at the next line that starts
    with '#', as in ReadInputFile and ReadSection. The second copy of a section
    is kept under 'NAME#2' so that it can be reported.
    """
    sections = {}
    i = 0
    while i < len(lines):
        text = lines[i].rstrip("\n").rstrip()
        # Rocflu compares the whole trimmed line: '# PICL', not '#PICL'
        if len(text.split()) == 2 and text == "# " + text.split()[1]:
            name = text.split()[1]
            body = []
            i += 1
            while i < len(lines) and not lines[i].startswith("#"):
                body.append((i + 1, lines[i].rstrip("\n")))
                i += 1
            key = name
            while key in sections:
                key += "#2"
            sections[key] = body
        i += 1
    return sections


def simple_value(body, key):
    """Value of KEY in a section read with ReadSection: (line number, value)."""
    for number, text in body:
        if text.startswith(key) and (len(text) == len(key) or text[len(key)] in " \t"):
            return number, read_number(text[len(key):])
    return None, None


class Validator:
    def __init__(self, inp):
        self.inp = inp
        self.case_dir = os.path.dirname(os.path.abspath(inp))
        self.errors = []

    def prohibit(self, condition, message, line=None):
        """Records an error when condition holds (as MFC's case validator)."""
        if condition:
            where = "%s:%d" % (os.path.basename(self.inp), line) if line else \
                os.path.basename(self.inp)
            self.errors.append("%s: %s" % (where, message))

    # -- # PICL section ------------------------------------------------------

    def read_picl(self, body):
        """Values given in the # PICL section ({key: (line, value)}), with
        every line Rocflu would ignore or misread reported."""
        known = PICL_KEYS + list(REMOVED_KEYS)
        given = {}
        for number, text in body:
            stripped = text.strip()
            if not stripped or stripped.startswith("!"):
                continue
            # The key the line means: the longest known key it starts with
            meant = None
            for key in known:
                if stripped.startswith(key) and \
                        (len(stripped) == len(key) or stripped[len(key)] in " \t"):
                    if meant is None or len(key) > len(meant):
                        meant = key
            if meant is None:
                word = stripped.split()[0]
                close = difflib.get_close_matches(word.upper(), PICL_KEYS, n=1)
                hint = " (did you mean %s?)" % close[0] if close else ""
                self.prohibit(True, "unknown key %s in # PICL%s; Rocflu ignores it"
                              % (word, hint), number)
                continue
            if meant in REMOVED_KEYS:
                self.prohibit(True, "%s: %s" % (meant, REMOVED_KEYS[meant]), number)
                continue
            if text[:1] in (" ", "\t"):
                self.prohibit(True, "%s starts with a space; Rocflu only reads keys "
                              "at the start of the line, so it ignores it" % meant, number)
                continue
            # The key Rocflu takes: the first key not yet read the line starts with
            taken = next((k for k in PICL_KEYS if k not in given and text.startswith(k)),
                         None)
            if taken is None:
                self.prohibit(True, "%s is given twice; Rocflu uses the first one "
                              "(line %d)" % (meant, given[meant][0]), number)
                continue
            if taken != meant:
                self.prohibit(True, "%s comes before %s; Rocflu reads this line as %s. "
                              "Put %s first" % (meant, taken, taken, taken), number)
                given[taken] = (number, None)
                continue
            value = read_number(text[len(meant):])
            self.prohibit(value is None, "%s has no number after it" % meant, number)
            given[meant] = (number, value)
        return given

    def check_picl(self, sections):
        bodies = [k for k in sections if k.split("#")[0] == "PICL"]
        self.prohibit(len(bodies) > 1, "the # PICL section appears more than once")
        if not bodies:
            return None
        given = self.read_picl(sections["PICL"])
        value = {k: (given[k][1] if k in given and given[k][1] is not None else d)
                 for k, d in PICL_DEFAULTS.items()}
        line = {k: given[k][0] for k in given}

        def flag_ok(key, lo, hi):
            v = value[key]
            if v != int(v):
                self.prohibit(True, "%s %g: must be a whole number" % (key, v), line.get(key))
                return False
            self.prohibit(not lo <= v <= hi, "%s %d: must be %s" % (
                key, v, "%d or %d" % (lo, hi) if hi == lo + 1 else
                "from %d to %d" % (lo, hi)), line.get(key))
            return lo <= v <= hi

        for key, (lo, hi) in PICL_FLAGS.items():
            flag_ok(key, lo, hi)
        if value["USED"] != 1:
            return value

        v = value
        # Values that are only read with their model on
        if v["QSFLUCT"] >= 1 or v["PSEUDOTURB"] >= 1:
            for key, (lo, hi) in DEPENDENT_FLAGS.items():
                flag_ok(key, lo, hi)
        if v["COLLISIONS"] >= 1:
            for key in ("KSP", "EREST"):
                self.prohibit(key not in given, "%s must be given when COLLISIONS >= 1"
                              % key, line.get("COLLISIONS"))
            ksp = given.get("KSP", (None, None))[1]
            erest = given.get("EREST", (None, None))[1]
            self.prohibit(ksp is not None and ksp <= 0,
                          "KSP %g: must be > 0 (time steps a contact lasts)" % (ksp or 0),
                          line.get("KSP"))
            self.prohibit(erest is not None and not 0 < erest <= 1,
                          "EREST %g: must be in (0, 1]" % (erest or 0), line.get("EREST"))
        if v["VISCOUSUNSTEADY"] >= 1:
            self.prohibit(v["VUSTEPS"] < 1 or v["VUSTEPS"] != int(v["VUSTEPS"]),
                          "VUSTEPS %g: must be a whole number >= 1 when "
                          "VISCOUSUNSTEADY >= 1" % v["VUSTEPS"], line.get("VUSTEPS"))
        self.prohibit(v["TEMP"] <= 0, "TEMP %g: the initial particle temperature must be "
                      "given and > 0" % v["TEMP"], line.get("TEMP"))
        self.prohibit(v["DENSITY"] <= 0, "DENSITY %g: the particle density [kg/m^3] must "
                      "be given and > 0 (it replaced DENS in # MATERIAL)" % v["DENSITY"],
                      line.get("DENSITY"))
        self.prohibit(v["SPECIFICHEAT"] <= 0, "SPECIFICHEAT %g: the particle specific heat "
                      "[J/(kg K)] must be given and > 0 (it replaced SPHT in # MATERIAL)"
                      % v["SPECIFICHEAT"], line.get("SPECIFICHEAT"))
        self.prohibit(v["NEIGHBORWIDTH"] < 0, "NEIGHBORWIDTH %g: must be >= 0"
                      % v["NEIGHBORWIDTH"], line.get("NEIGHBORWIDTH"))

        # Models that need another model
        self.prohibit(v["QSFLUCT"] >= 1 and v["QSFORCE"] == 0,
                      "QSFLUCT %d needs a quasi-steady force (QSFORCE >= 1): the "
                      "fluctuation is added to it" % v["QSFLUCT"], line.get("QSFLUCT"))
        self.prohibit(v["PSEUDOTURB"] == 1 and v["QSFORCE"] == 0,
                      "PSEUDOTURB 1 needs a quasi-steady force (QSFORCE >= 1): the "
                      "pseudo-turbulence model uses the drag", line.get("PSEUDOTURB"))
        self.prohibit(v["PSEUDOTURB"] == 1 and v["QSFLUCT"] != 2,
                      "PSEUDOTURB 1 needs QSFLUCT 2 (Osnes fluctuations)",
                      line.get("PSEUDOTURB"))
        self.prohibit(v["PSEUDOTURB"] == 1 and v["HEATTRANSFER"] == 0,
                      "PSEUDOTURB 1 needs HEATTRANSFER >= 1 for now: its heat flux "
                      "divides by the Nusselt number (until the PT refactor)",
                      line.get("PSEUDOTURB"))
        self.prohibit(v["DU HEATTRANSFER"] == 1 and v["VISCOUSUNSTEADY"] == 0,
                      "DU HEATTRANSFER 1 needs VISCOUSUNSTEADY >= 1 (they share "
                      "the history)", line.get("DU HEATTRANSFER"))
        self.prohibit(v["VELOCITYMODEL"] == 0 and v["COLLISIONS"] >= 1,
                      "COLLISIONS %d with VELOCITYMODEL 0: particles that keep their "
                      "position cannot collide" % v["COLLISIONS"], line.get("COLLISIONS"))
        self.prohibit(v["AMFORCE"] == 2,
                      "AMFORCE 2 (Briney added mass) is not validated yet",
                      line.get("AMFORCE"))
        return v

    # -- other sections and files -----------------------------------------------

    def check_with_particles(self, sections):
        body = sections.get("TIMESTEP", [])
        number, rk = simple_value(body, "RKSCHEME")
        rk = RK4 if rk is None else int(abs(rk) + 0.5)
        self.prohibit(rk != RK3, "particles need the 3-stage Runge-Kutta scheme: "
                      "RKSCHEME 2 in # TIMESTEP (%s)" % ("now %d" % rk if number else
                                                         "not given, Rocflu uses 1, RK4"),
                      number)
        self.prohibit("VISCMODEL" not in sections,
                      "particles need a # VISCMODEL section with VISCOSITY")
        if "VISCMODEL" in sections:
            number, visc = simple_value(sections["VISCMODEL"], "VISCOSITY")
            self.prohibit(visc is None, "particles need VISCOSITY in # VISCMODEL")
        number, axi = simple_value(sections.get("NUMERICS", []), "AXIFLAG")
        self.prohibit(axi == 1, "axisymmetric runs (AXIFLAG 1) are not available "
                      "with particles", number)

        def present(name):
            return os.path.exists(os.path.join(self.case_dir, name))
        self.prohibit(not present("points.dat"),
                      "no points.dat (the particles) in %s" % self.case_dir)
        self.prohibit(present("filein.vtk"),
                      "filein.vtk is no longer read: set PARTWALL 1 on the wall "
                      "patches in the .bc file and remove filein.vtk")
        self.prohibit(present("par00002.vtu"),
                      "par00002.vtu from an earlier run is in the case directory: it "
                      "starts a particle restart, which is not available yet. Remove "
                      "the par*.vtu files (and the .rin file) to start a new run")

    def run(self):
        try:
            with open(self.inp) as f:
                lines = f.readlines()
        except OSError as e:
            self.prohibit(True, "cannot read: %s" % e)
            return self.errors
        sections = read_sections(lines)
        picl = self.check_picl(sections)
        if picl is not None and picl["USED"] == 1:
            self.check_with_particles(sections)
        return self.errors


def validate(inp):
    """List of error messages for the case <dir>/<case>.inp (empty: fine)."""
    return Validator(inp).run()


def main(argv):
    if len(argv) != 2 or argv[1] in ("-h", "--help"):
        print(__doc__.strip())
        return 0 if len(argv) == 2 else 2
    errors = validate(argv[1])
    for message in errors:
        print("case_validator: ERROR: " + message, file=sys.stderr)
    if errors:
        print("case_validator: %d error%s in %s" % (len(errors), "s" if len(errors) > 1
                                                    else "", argv[1]), file=sys.stderr)
        return 1
    print("case_validator: %s: OK" % argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
