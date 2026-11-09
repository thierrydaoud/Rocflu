"""Write the last Rocflu dump of a run in the plain columns metrics/compare.py reads.

    python3 to_columns.py <run directory> > result.txt

Reads post/p*/<step>.vtu (written by rflupost) and converts from SI to the dimensionless units of
../problems/: pressure scale 1e5 Pa, density scale 1 kg/m3, length scale 1 m, so velocity scale
sqrt(1e5) m/s and time scale 1/sqrt(1e5) s. Refuses runs whose rflump.log does not say "Finished",
because the solver can stop on an error with exit status 0.
"""
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from rocpicl_tests import outputs  # noqa: E402
from units import P_SCALE, RHO_SCALE, U_SCALE  # noqa: E402


def main(rd):
    log = os.path.join(rd, "rflump.log")
    if not os.path.exists(log) or "Finished" not in open(log, errors="replace").read():
        sys.exit("%s: rflump.log does not say Finished" % rd)
    last = max(int(os.path.basename(f)[:-4]) for f in glob.glob(os.path.join(rd, "post", "p0", "*.vtu")))
    rows, t, step = [], None, None
    for f in sorted(glob.glob(os.path.join(rd, "post", "p*", "%08d.vtu" % last))):
        t, step, _, arr = outputs.read_post_vtu(f)
        centres = outputs.read_post_vtu(f, centroids=True)
        for c, r, p, v in zip(centres, arr["Gas Density"], arr["Pressure"], arr["Velocity"]):
            rows.append((c[0], c[1], c[2], r / RHO_SCALE, v[0] / U_SCALE, p / P_SCALE, v[1] / U_SCALE, v[2] / U_SCALE))
    print("# Rocflu run %s, dump %d" % (os.path.abspath(rd), last))
    print("# t %.17g\n# steps %d\n# columns x y z rho u p v w" % (t * U_SCALE, step))
    for row in rows:
        print("%.10f %.10f %.10f %.15e %.15e %.15e %.15e %.15e" % row)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
