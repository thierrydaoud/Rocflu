"""Markdown tables of rung 1 from a case root written by make_cases.py and run by run_cases.sh.

    python3 summarize.py <root> > rung1_tables.md

Every number comes from metrics/compare.py; a set counts as failed if compare.py fails any check.
"""
import contextlib
import io
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "metrics"))
import compare  # noqa: E402
import godunov  # noqa: E402  (compare.py puts ../exact on the path)

SCHEMES = {1: "Roe", 3: "HLLC", 4: "AUSM+", 5: "AUSM+UP"}


def main(root):
    metas = [json.load(open(os.path.join(root, c, "meta.json"))) | {"case": c}
             for c in open(os.path.join(root, "cases.txt")).read().split()]
    sets = sorted({(m["problem"], m["discr"], m["order"], m["dim"]) for m in metas})
    print("| Problem | Scheme | Order | Dim | L1(ρ), n = 800 | Rate 400–800 | Shock, cells (max) "
          "| Conservation (max, undisturbed ends) | Transverse / a | ÷ Godunov, n = 800 | Result |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for prob, discr, order, dim in sets:
        cases = sorted((m for m in metas if (m["problem"], m["discr"], m["order"], m["dim"]) == (prob, discr, order, dim)),
                       key=lambda m: m["n"])
        paths = [os.path.join(root, m["case"], "result.txt") for m in cases]
        missing = [m["n"] for m, p in zip(cases, paths) if not os.path.exists(p)]
        row = "| %s | %s | %d | %dD |" % (prob, SCHEMES[discr], order, dim)
        if missing:
            print(row + " | | | | | | did not finish: n = %s |" % ", ".join(map(str, missing)))
            continue
        with contextlib.redirect_stdout(io.StringIO()):
            fails = compare.report(prob, order, paths)
        runs = [compare.measure(prob, p) for p in paths]
        fine = runs[-1]
        rate = np.log2(runs[-2]["L1"][0] / fine["L1"][0])
        clean = [max(abs(e) for e in r["conservation"]) for r in runs if r["end_cells"] <= compare.LIMITS["end_cells"]]
        cons = "%.1e" % max(clean) if clean else "–"
        trans = "%.1e" % max(r["transverse"] for r in runs)
        ratio = "%.3f" % (fine["L1"][0] / godunov.run(prob, fine["n"], steps=fine["steps"])["L1"][0]) if order == 1 else "–"
        shock = max(abs(r["shock_cells"]) for r in runs)
        print(row + " %.4e | %.2f | %.2f | %s | %s | %s | %s |" % (
            fine["L1"][0], rate, shock, cons, trans, ratio, "pass" if fails == 0 else "**%d failed**" % fails))
    print()
    print("| Problem | Scheme | Order | 2D vs 3D, n = 800 |")
    print("|---|---|---|---|")
    for prob, discr, order in sorted({s[:3] for s in sets}):
        p2, p3 = [os.path.join(root, "%s_d%d_o%d_dim%d_n800" % (prob, discr, order, d), "result.txt") for d in (2, 3)]
        if os.path.exists(p2) and os.path.exists(p3):
            l2, l3 = compare.measure(prob, p2)["L1"][0], compare.measure(prob, p3)["L1"][0]
            print("| %s | %s | %d | %.1e |" % (prob, SCHEMES[discr], order, abs(l2 - l3) / l3))
        else:
            print("| %s | %s | %d | missing |" % (prob, SCHEMES[discr], order))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
