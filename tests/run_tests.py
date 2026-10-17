#!/usr/bin/env python3
"""Rocflu + ppiclF regression tests.

    ./run_tests.py list                      # show cases
    ./run_tests.py run                       # run all, compare to goldens
    ./run_tests.py run -o shocktube/two-way  # cases whose name contains a pattern
    ./run_tests.py run --generate            # (re)write goldens from this build
    ./run_tests.py run --exe ../../other_tree/build/bin

Runs from a login node submit one short srun job per case (-A/-p below);
inside an allocation (salloc) they run as job steps. Set-up is serial.
"""

import argparse
import datetime
import json
import os
import platform
import shutil
import subprocess
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "utilities", "scripts"))

from rocpicl_tests import case as rcase          # noqa: E402
from rocpicl_tests import checks, state           # noqa: E402
from rocpicl_tests.cases import CASES             # noqa: E402
import case_validator                             # noqa: E402

TESTS_DIR = rcase.TESTS_DIR

# Colours only on a terminal, and not when NO_COLOR is set (as rocflu.sh)
COLOURS = (sys.stdout.isatty() and not os.environ.get("NO_COLOR")
           and os.environ.get("TERM", "") != "dumb")
STATUS_COLOUR = {"PASS": "32", "FAIL": "1;31", "XFAIL": "33", "XPASS": "36"}


def cases_text(n):
    return "%d case%s" % (n, "" if n == 1 else "s")


def paint(text, code):
    """text in an ANSI colour (e.g. "32" green, "1;31" bold red) on a terminal."""
    return "\033[%sm%s\033[0m" % (code, text) if COLOURS else text
GOLDENS = os.path.join(TESTS_DIR, "goldens")
WORK = os.path.join(TESTS_DIR, "_work", "runs")


def golden_dir(case):
    return os.path.join(GOLDENS, case.name)


def default_exe_dir():
    # The last configuration built by ./rocflu.sh build
    return os.path.join(os.path.dirname(TESTS_DIR), "build", "bin")


def git_info(exe_dir):
    def git(*args):
        try:
            return subprocess.run(["git", "-C", exe_dir] + list(args), capture_output=True,
                                  text=True).stdout.strip()
        except OSError:
            return ""
    # Sources changed since the last commit
    changed = git("diff", "--name-only", "HEAD").split("\n")
    sources = [f for f in changed if f.endswith((".F90", ".f", ".F", ".h", ".C"))]
    return {"commit": git("rev-parse", "HEAD"), "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "modified_sources": sources}


def run_case(case, args):
    """Set up, run and check one case. Returns a result dict."""
    rundir = os.path.join(WORK, case.name)
    # errors: run failures and golden mismatches (always fatal)
    # physics: physics/cross-check failures (expected while known_issue is set)
    res = {"case": case, "errors": [], "physics": [], "final": None}
    slurm = {"account": args.account, "partition": args.partition,
             "local_ranks": args.local_ranks}
    try:
        # One retry from a clean set-up: a failed srun can leave partial
        # dumps that a rerun in the same directory would restart from
        for attempt in (1, 2):
            rcase.setup(case, rundir, args.exe)
            # The case's input must pass the checks a user's case gets
            invalid = case_validator.validate(os.path.join(rundir, rcase.CASENAME + ".inp"))
            if invalid:
                res["errors"] += ["case validator: " + e for e in invalid]
                return res
            try:
                rcase.run_solver(case, rundir, args.exe, slurm)
                break
            except rcase.RunError as e:
                if attempt == 2:
                    raise
                res["note"] = "retried after: %s" % e
        init = state.snapshot(rundir, case, 0)
        final = state.snapshot(rundir, case, -1)
        res["final"] = final
        res["errors"] += checks.check_restart_data(rundir)
        if case.particles:
            res["errors"] += checks.check_particle_dump(rundir)
        # The post program must give what the old-style output gives
        post_errors = rcase.run_post(rundir, args.exe)
        res["errors"] += post_errors or checks.check_post(rundir)
        if case.restart_test:
            res["errors"] += rcase.restart_test(case, rundir, args.exe, slurm)

        if not state.all_finite(final):
            res["errors"].append("non-finite values in final state")
        for name in case.checks:
            res["physics"] += checks.PHYSICS[name](init, final, case)

        gpath = os.path.join(golden_dir(case), "golden.json")
        if args.generate:
            if not res["errors"] and (not res["physics"] or case.known_issue or args.force):
                write_golden(case, final, args)
                res["note"] = "golden written"
            else:
                res["note"] = "golden NOT written (checks failed; --force to override)"
        elif not os.path.exists(gpath):
            res["errors"].append("no golden (run with --generate)")
        else:
            with open(gpath) as f:
                golden = json.load(f)
            if golden["uuid"] != case.uuid:
                res["errors"].append("case definition changed since golden was written "
                                     "(regenerate with --generate)")
            else:
                res["errors"] += checks.compare_states(final, golden["state"])
    except Exception as e:  # report and keep going with other cases
        res["errors"].append("%s: %s" % (type(e).__name__, e))
        if args.verbose:
            traceback.print_exc()
    return res


def status(res):
    if res["errors"]:
        return "FAIL"
    if res["case"].known_issue:
        return "XFAIL" if res["physics"] else "XPASS"
    return "FAIL" if res["physics"] else "PASS"


def write_golden(case, final, args):
    d = golden_dir(case)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "golden.json"), "w") as f:
        state_out = {k: v for k, v in final.items() if not k.startswith("_")}
        json.dump({"uuid": case.uuid, "spec": case.spec(), "state": state_out}, f, indent=0)
    meta = {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "host": platform.node(),
        "exe_dir": os.path.abspath(args.exe),
        "git": git_info(args.exe),
        "modules": os.environ.get("LOADEDMODULES", ""),
    }
    with open(os.path.join(d, "metadata.json"), "w") as f:
        json.dump(meta, f, indent=1)


def cross_checks(results):
    """same_as: particles must match the referenced case (other rank count)."""
    by_name = {r["case"].name: r for r in results}
    for r in results:
        ref_name = r["case"].same_as
        if not ref_name or ref_name not in by_name or r["final"] is None:
            continue
        ref = by_name[ref_name]
        if ref["final"] is None:
            continue
        errs = checks.compare_particles_only(r["final"], ref["final"], ref_name,
                                             r["case"].same_as_fields)
        r["physics"] += ["differs from %s: %s" % (ref_name, e) for e in errs]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("command", choices=["list", "run"])
    p.add_argument("-o", "--only", nargs="*", default=[], help="name substrings to select")
    p.add_argument("--generate", action="store_true", help="write goldens instead of comparing")
    p.add_argument("--force", action="store_true", help="with --generate: write even if checks fail")
    p.add_argument("--exe", default=default_exe_dir(),
                   help="directory with the rflu* executables (default: build/bin of this tree)")
    p.add_argument("-j", "--jobs", type=int, default=8, help="cases run concurrently")
    p.add_argument("--account", default=os.environ.get("ROCPICL_TEST_ACCOUNT", "bala1s"))
    p.add_argument("--partition", default=os.environ.get("ROCPICL_TEST_PARTITION", "hpg-dev"))
    p.add_argument("--local-ranks", type=int, default=0,
                   help="run cases with at most this many ranks with mpirun on this node")
    p.add_argument("--keep", action="store_true", help="keep run directories of passing cases")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()
    args.exe = os.path.abspath(args.exe)

    selected = [c for c in CASES if not args.only or any(o in c.name for o in args.only)]
    # A same_as reference must run too
    names = {c.name for c in selected}
    selected += [c for c in CASES if c.name not in names
                 and any(s.same_as == c.name for s in selected)]

    if args.command == "list":
        for c in selected:
            print("%-48s np=%d  %s" % (c.name, c.nranks, c.description))
        return 0

    for exe in ("rflumap", "rflupart", "rfluinit", "rflump"):
        if not os.access(os.path.join(args.exe, exe), os.X_OK):
            print(paint("Missing executable %s in %s (build with ./rocflu.sh build)"
                        % (exe, args.exe), "1;31"))
            return 2

    print(paint("Running %s" % cases_text(len(selected)), "1") + " with %s (%s)" % (
        args.exe, "generate" if args.generate else "compare"))
    start = datetime.datetime.now()
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda c: run_case(c, args), selected))
    cross_checks(results)

    counts = {}
    for r in sorted(results, key=lambda r: r["case"].name):
        st = status(r)
        counts[st] = counts.get(st, 0) + 1
        print("%s  %s%s" % (paint("%-5s" % st, STATUS_COLOUR[st]), r["case"].name,
                             paint("   [" + r["note"] + "]", "2") if r.get("note") else ""))
        if st in ("XFAIL", "XPASS"):
            print(paint("         known issue: " + r["case"].known_issue, "33"))
        if st == "XPASS":
            print(paint("         physics checks now pass: remove known_issue from the case", "36"))
        for e in r["errors"] + r["physics"]:
            print(paint("         - " + e, "31"))
        if st == "PASS" and not args.keep:
            shutil.rmtree(os.path.join(WORK, r["case"].name), ignore_errors=True)
            shutil.rmtree(os.path.join(WORK, r["case"].name) + "_restart", ignore_errors=True)
        elif st != "PASS":
            print(paint("         run dir: " + os.path.join(WORK, r["case"].name), "2"))
    seconds = (datetime.datetime.now() - start).total_seconds()
    print("\n" + ", ".join(paint("%d %s" % (counts[k], k), STATUS_COLOUR[k])
                           for k in ("PASS", "XFAIL", "XPASS", "FAIL") if k in counts))
    if counts.get("FAIL"):
        print(paint("✗ %d of %s failed" % (counts["FAIL"], cases_text(len(results))), "1;31")
              + " (%.0f s)" % seconds)
    else:
        print(paint("✔ All %s passed" % cases_text(len(results)) if len(results) > 1
                    else "✔ The case passed", "1;32") + " in %.0f s" % seconds)
    return 1 if counts.get("FAIL") else 0


if __name__ == "__main__":
    sys.exit(main())
