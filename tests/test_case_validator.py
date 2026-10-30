#!/usr/bin/env python3
"""Tests of utilities/scripts/case_validator.py (no cluster needed).

    python3 tests/test_case_validator.py

Each test writes the template input (tests/templates/rectshktb.inp) with a
few changes and a points.dat into a temporary directory and checks the
errors the validator reports.
"""

import os
import shutil
import sys
import tempfile
import unittest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TESTS_DIR)
sys.path.insert(0, os.path.join(TESTS_DIR, "..", "utilities", "scripts"))

import case_validator                          # noqa: E402
from rocpicl_tests.case import bc_text, edit_inp   # noqa: E402

with open(os.path.join(TESTS_DIR, "templates", "rectshktb.inp")) as f:
    TEMPLATE = f.read()
with open(os.path.join(TESTS_DIR, "templates", "rectshktb.bc")) as f:
    TEMPLATE_BC = f.read()

# Two particles (x y z vx vy vz dp), one with an initial velocity
POINTS = "! test\n2\n0.001 0.002 0.003 0 0 0 1e-4\n0.004 0.005 0.006 1.0 0.0 0.0 1e-4\n"


def write_box_grid(path, lo, hi):
    """The vertex part of a CENTAUR ASCII grid: the 8 corners of a box."""
    corners = [(x, y, z) for x in (lo, hi) for y in (lo, hi) for z in (lo, hi)]
    with open(path, "w") as f:
        f.write("box\n%16d\n" % len(corners))
        for axis in range(3):
            values = [c[axis] for c in corners]
            for i in range(0, len(values), 5):
                f.write("".join("%16.9E" % v for v in values[i:i + 5]) + "\n")


def write_patch_grid(path, faces):
    """A CENTAUR ASCII grid with no cells and len(faces) patches; faces[i] is
    the number of triangles of patch i + 1 (0: an empty patch)."""
    def block(values):
        return "".join("".join("%16d" % v for v in values[i:i + 10]) + "\n"
                       for i in range(0, len(values), 10))
    ends = [sum(faces[:i + 1]) for i in range(len(faces))]
    with open(path, "w") as f:
        f.write("patches\n%16d\n" % 1)
        f.write(("%16.9E\n" % 0.0) * 3 + block([0]))
        f.write(("%16d\n" % 0) * 4)
        f.write("%16d\n" % len(faces))
        f.write(block([1] * len(faces)) + block(ends) + block([0] * len(faces)))


class CaseValidatorTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir)

    def errors(self, overrides=None, text=None, files=("points.dat",), points=POINTS,
               grid=None, bc=TEMPLATE_BC, patches=None):
        """Errors for the template with {(SECTION, KEY): value} changes, or
        for the given input text, with the given files in the case directory
        (points.dat holds `points`; a grid file is written when grid is a
        (min, max) box, or with the patches of write_patch_grid; the .bc file
        holds `bc`, none when it is None)."""
        if text is None:
            text = edit_inp(TEMPLATE, overrides or {})
        inp = os.path.join(self.dir, "case.inp")
        with open(inp, "w") as f:
            f.write(text)
        if bc is not None:
            with open(os.path.join(self.dir, "case.bc"), "w") as f:
                f.write(bc)
        if patches:
            write_patch_grid(os.path.join(self.dir, "case.hyb.asc"), patches)
        for name in files:
            with open(os.path.join(self.dir, name), "w") as f:
                f.write(points if name == "points.dat" else "")
        if grid:
            write_box_grid(os.path.join(self.dir, "case.hyb.asc"), *grid)
        return case_validator.validate(inp)

    def assertError(self, errors, *words):
        """One of the errors contains all the words."""
        self.assertTrue(any(all(w in e for w in words) for e in errors),
                        "no error with %s in:\n  %s" % (words, "\n  ".join(errors)))

    # -- accepted cases ----------------------------------------------------------

    def test_template_is_fine(self):
        self.assertEqual(self.errors(), [])

    def test_pseudo_turbulence_with_what_it_needs(self):
        for heattransfer in (0, 2):
            self.assertEqual(self.errors({("PICL", "PSEUDOTURB"): 1, ("PICL", "QSFLUCT"): 1,
                                          ("PICL", "HEATTRANSFER"): heattransfer}), [])
        self.assertEqual(self.errors({("PICL", "PSEUDOTURB"): 2, ("PICL", "QSFLUCT"): 1,
                                      ("PICL", "HEATTRANSFER"): 5}), [])

    def test_unused_parameters_left_in_the_file(self):
        # Collisions off: KSP and EREST are not read, any value is fine
        self.assertEqual(self.errors({("PICL", "KSP"): -5, ("PICL", "EREST"): 3,
                                      ("PICL", "QSFLUCTFILT"): 7}), [])

    def test_particles_off_skips_model_checks(self):
        self.assertEqual(self.errors({("PICL", "USED"): 0, ("PICL", "QSFORCE"): 0,
                                      ("PICL", "QSFLUCT"): 1}, files=()), [])

    # -- keys ----------------------------------------------------------------------

    def test_removed_keys(self):
        for key, words in (("FILTERWIDTH 0.002", ("FILTERWIDTH", "remove")),
                           ("STATIONARY 1", ("VELOCITYMODEL 0",)),
                           ("REACTIVE 0", ("burn-rate",)),
                           ("ANGULARPERIODIC 0", ("angular",)),
                           ("OUTGAS 1", ("GAS_WRT",))):
            text = TEMPLATE.replace("# PICL\n", "# PICL\n%s\n" % key)
            self.assertError(self.errors(text=text), *words)

    def test_unknown_key_with_hint(self):
        text = TEMPLATE.replace("# PICL\n", "# PICL\nQSFLUC 1\n")
        self.assertError(self.errors(text=text), "unknown key QSFLUC", "QSFLUCT")

    def test_key_after_a_space_is_ignored_by_rocflu(self):
        text = TEMPLATE.replace("\nPGFORCE", "\n  PGFORCE")
        self.assertError(self.errors(text=text), "PGFORCE starts with a space")

    def test_key_given_twice(self):
        text = TEMPLATE.replace("# PICL\n", "# PICL\nQSFORCE 1\n")
        self.assertError(self.errors(text=text), "QSFORCE is given twice")

    def test_key_read_as_another(self):
        # QSFLUCTFILT before QSFLUCT is read as QSFLUCT by Rocflu
        text = TEMPLATE.replace("# PICL\n", "# PICL\nQSFLUCTFILT 1\n")
        self.assertError(self.errors(text=text), "QSFLUCTFILT comes before QSFLUCT")

    # -- values --------------------------------------------------------------------

    def test_flag_out_of_range(self):
        self.assertError(self.errors({("PICL", "QSFORCE"): 7}), "QSFORCE 7", "from 0 to 3")
        self.assertError(self.errors({("PICL", "VELOCITYMODEL"): 2}), "VELOCITYMODEL 2")

    def test_flag_must_be_whole(self):
        self.assertError(self.errors({("PICL", "QSFORCE"): 2.5}), "whole number")

    def test_collision_parameters(self):
        self.assertError(self.errors({("PICL", "COLLISIONS"): 1, ("PICL", "EREST"): 0}),
                         "EREST")
        self.assertError(self.errors({("PICL", "COLLISIONS"): 1, ("PICL", "KSP"): 0}),
                         "KSP")
        text = edit_inp(TEMPLATE, {("PICL", "COLLISIONS"): 1}).replace("\nKSP ", "\n!KSP ")
        self.assertError(self.errors(text=text), "KSP must be given")

    def test_temperature_needed(self):
        self.assertError(self.errors({("PICL", "TEMP"): 0}), "TEMP")

    def test_particle_material_needed(self):
        self.assertError(self.errors({("PICL", "DENSITY"): 0}), "DENSITY 0")
        text = TEMPLATE.replace("\nSPECIFICHEAT", "\n!SPECIFICHEAT")
        self.assertError(self.errors(text=text), "SPECIFICHEAT 0", "replaced SPHT")

    # -- combinations ----------------------------------------------------------------

    def test_fluctuations_need_quasi_steady_force(self):
        self.assertError(self.errors({("PICL", "QSFORCE"): 0, ("PICL", "QSFLUCT"): 1}),
                         "QSFLUCT 1 needs a quasi-steady force")

    def test_pseudo_turbulence_needs(self):
        errors = self.errors({("PICL", "PSEUDOTURB"): 1, ("PICL", "QSFLUCT"): 0})
        self.assertError(errors, "PSEUDOTURB 1 needs QSFLUCT 1")
        for heattransfer in (0, 2):
            self.assertError(self.errors({("PICL", "PSEUDOTURB"): 2, ("PICL", "QSFLUCT"): 1,
                                          ("PICL", "HEATTRANSFER"): heattransfer}),
                             "PSEUDOTURB 2 needs HEATTRANSFER 5")
        for qsforce in (0, 1, 2):
            self.assertError(self.errors({("PICL", "PSEUDOTURB"): 1, ("PICL", "QSFLUCT"): 1,
                                          ("PICL", "HEATTRANSFER"): 2,
                                          ("PICL", "QSFORCE"): qsforce}),
                             "PSEUDOTURB 1 needs QSFORCE 3")

    def test_old_fluctuation_numbering(self):
        self.assertError(self.errors({("PICL", "QSFLUCT"): 2}),
                         "QSFLUCT 2 (Osnes) is now QSFLUCT 1")

    def test_unsteady_heat_transfer_needs_history(self):
        text = TEMPLATE.replace("# PICL\n", "# PICL\nDU HEATTRANSFER 1\n")
        self.assertError(self.errors(text=text), "DU HEATTRANSFER 1 needs VISCOUSUNSTEADY")

    def test_viscous_unsteady_steps(self):
        self.assertError(self.errors({("PICL", "VISCOUSUNSTEADY"): 1, ("PICL", "VUSTEPS"): 0}),
                         "VUSTEPS")

    def test_fixed_particles_cannot_collide(self):
        self.assertError(self.errors({("PICL", "VELOCITYMODEL"): 0, ("PICL", "COLLISIONS"): 1}),
                         "VELOCITYMODEL 0")

    def test_binary_added_mass_not_validated(self):
        self.assertError(self.errors({("PICL", "AMFORCE"): 2}), "AMFORCE 2")

    # -- other sections and files -------------------------------------------------------

    def test_runge_kutta_scheme(self):
        self.assertError(self.errors({("TIMESTEP", "RKSCHEME"): 1}), "RKSCHEME 2", "now 1")
        text = TEMPLATE.replace("\nRKSCHEME", "\n!RKSCHEME")
        self.assertError(self.errors(text=text), "RKSCHEME 2", "not given")

    def test_viscosity_needed(self):
        text = TEMPLATE.replace("\nVISCOSITY", "\n!VISCOSITY")
        self.assertError(self.errors(text=text), "VISCOSITY")

    # -- points.dat ------------------------------------------------------------------

    def test_points_old_format_with_material(self):
        old = "2\nAL 0.001 0.002 0.003 1e-4\nAL 0.004 0.005 0.006 1e-4\n"
        self.assertError(self.errors(points=old), "points.dat:2", "material name was removed")
        short = "1\n0.001 0.002 0.003 1e-4\n"
        self.assertError(self.errors(points=short), "points.dat:2", "velocity now comes before dp")

    def test_points_count_and_lines(self):
        line = "0.001 0.002 0.003 0 0 0 1e-4\n"
        self.assertError(self.errors(points="3\n" + line), "3 particles announced", "only 1")
        self.assertError(self.errors(points="1\n" + line + line),
                         "1 line after the 1 particles")
        self.assertError(self.errors(points="1\n0.001 0.002 0.003 1.0 0 1e-4\n"),
                         "points.dat:2", "6 entries")
        self.assertError(self.errors(points="1\n0.001 0.002 0.003 0 0 0 0\n"), "diameter 0")
        self.assertError(self.errors(points="x\n"), "number of particles")

    def test_points_in_the_domain(self):
        self.assertEqual(self.errors(grid=(0.0, 0.01)), [])
        self.assertError(self.errors(grid=(0.0, 0.005)), "points.dat:4",
                         "outside the fluid domain")

    def test_species_and_mixtures_disabled(self):
        self.assertError(self.errors({("SPECIES", "USED"): 1}), "species transport is disabled")
        self.assertError(self.errors({("MIXTURE", "GASMODEL"): 7}), "GASMODEL 7")
        self.assertError(self.errors({("PICL", "USED"): 0, ("MIXTURE", "GASMODEL"): 3},
                                     files=()), "GASMODEL 3")

    def test_files_in_case_directory(self):
        self.assertError(self.errors(files=()), "no points.dat")
        self.assertError(self.errors(files=("points.dat", "filein.vtk")), "filein.vtk")

    def test_files_of_an_earlier_run(self):
        # They do not matter: a new run removes them
        self.assertEqual(self.errors(files=("points.dat", "par00002.vtu")), [])

    def test_post_field_keys(self):
        text = TEMPLATE.replace("# POST", "# POST\nTEMPERATURE_WRT 1\nPT_WRT 2", 1)
        errors = self.errors(text=text)
        self.assertError(errors, "PT_WRT", "0 or 1")
        self.assertFalse(any("TEMPERATURE_WRT" in e for e in errors))

    def test_end_of_run_and_start_step(self):
        text = TEMPLATE.replace("# TIMESTEP", "# TIMESTEP\nPRNTENDTIME 2", 1)
        self.assertError(self.errors(text=text), "PRNTENDTIME", "removed")
        self.assertEqual(self.errors({("TIMESTEP", "START_STEP"): 0}), [])
        self.assertEqual(self.errors({("TIMESTEP", "START_STEP"): -1}), [])
        self.assertError(self.errors({("TIMESTEP", "START_STEP"): -2}), "START_STEP")
        self.assertError(self.errors({("TIMESTEP", "START_STEP"): 1.5}), "START_STEP")

    # -- the .bc file ----------------------------------------------------------------

    def bc_errors(self, old, new, **kwargs):
        """Errors of the template .bc with `old` replaced by `new`."""
        self.assertIn(old, TEMPLATE_BC)
        return [e for e in self.errors(bc=TEMPLATE_BC.replace(old, new, 1), **kwargs)
                if "case.bc" in e]

    def test_bc_template_and_suite_files_are_fine(self):
        for kind in ("shocktube", "periodic", "closed", "partwall"):
            errors = self.errors(bc=bc_text(kind), patches=[1] * 6)
            self.assertEqual([e for e in errors if "case.bc" in e], [], kind)

    def test_bc_file_missing(self):
        self.assertError(self.errors(bc=None), "no case.bc")

    def test_patch_needs_two_values(self):
        self.assertEqual(self.bc_errors("PATCH      2 2", "PATCH      2 2 ! outflow"), [])
        for line in ("PATCH      2", "PATCH      2 ! outflow"):
            self.assertError(self.bc_errors("PATCH      2 2", line), "PATCH needs two",
                             "PATCH 3 3")
        self.assertError(self.bc_errors("PATCH      2 2", "PATCH      2 1"), "first <= last")

    def test_bc_keys(self):
        self.assertError(self.bc_errors("NSCBCK     0.0", "PRESSURE 1e5"),
                         "unknown key PRESSURE", "PRESS")
        self.assertError(self.bc_errors("NSCBCK     0.0", "MVPATCH 0"), "MVPATCH", "remove")
        self.assertError(self.bc_errors("NSCBCK     0.0", "COUPLED 2"), "COUPLED", "GENX")
        self.assertError(self.bc_errors("NSCBCK     0.0", "  NSCBCK 0.0"), "starts with a space")
        self.assertError(self.bc_errors("NSCBCK     0.0", "NSCBCK 0.0\nNSCBCK 0.1"),
                         "NSCBCK is given twice")
        self.assertError(self.bc_errors("NSCBCK     0.0", "NSCBCK"), "NSCBCK has no number")
        # A key of another section
        self.assertError(self.bc_errors("NSCBCK     0.0", "PARTWALL 1"),
                         "unknown key PARTWALL in # BC_OUTFLOW")

    def test_bc_values(self):
        self.assertError(self.bc_errors("KIND       1", "KIND       2"), "KIND 2", "0 or 1")
        self.assertError(self.bc_errors("REFLECT    0\nTHRUST", "REFLECT    3\nTHRUST"), "REFLECT 3")
        self.assertError(self.bc_errors("TYPE       0", "TYPE       4"), "TYPE 4", "0 or 1 or 2")
        self.assertError(self.bc_errors("ORDER      1", "ORDER      3"), "ORDER 3")

    def test_bc_required_keys(self):
        self.assertError(self.bc_errors("NAME       Absorbing_Wall\n", ""), "needs NAME")
        # Subsonic outflow needs the pressure
        self.assertError(self.bc_errors("TYPE       0", "TYPE       1"), "BC_OUTFLOW needs PRESS")
        self.assertError(self.bc_errors("TEMP            348", "!TEMP            348"),
                         "BC_INFLOW_VELTEMP needs TEMP")
        # ... unless the values come from a file
        self.assertEqual(self.bc_errors("TEMP            348", "DISTRIB 1\n!TEMP 348"), [])

    def test_bc_sections(self):
        self.assertError(self.bc_errors("# BC_OUTFLOW", "# BC_OUTFLO"),
                         "unknown section # BC_OUTFLO", "BC_OUTFLOW")
        self.assertError(self.bc_errors("# BC_OUTFLOW", "#BC_OUTFLOW"), "one space")
        self.assertError(self.bc_errors("# BC_OUTFLOW", "# BC_INFLOW"), "old name",
                         "BC_INFLOW_TOTANG")
        self.assertError(self.bc_errors("# END", ""), "no '# END'")
        self.assertError(self.bc_errors("# BC_OUTFLOW", "PATCH 2 2\n# BC_OUTFLOW"),
                         "outside any section")

    def test_section_not_closed(self):
        # The header of the next section closes the one above, which is skipped
        errors = self.bc_errors("REFLECT     0\n#\n", "REFLECT     0\n")
        self.assertError(errors, "# BC_OUTFLOW", "ends the section # BC_INFLOW_VELTEMP")

    def test_patches_covered_once(self):
        self.assertError(self.bc_errors("PATCH      2 2", "PATCH      1 2"),
                         "patch 1 is also in the section")
        self.assertError(self.bc_errors("PATCH      2 2", "PATCH      7 7"),
                         "no boundary condition for patch 2")
        # Patch numbers from the grid
        self.assertEqual(self.bc_errors("PATCH      2 2", "PATCH      2 2", patches=[1] * 6), [])
        self.assertError(self.bc_errors("PATCH      2 2", "PATCH      2 2", patches=[1] * 7),
                         "no boundary condition for patch 7")
        self.assertError(self.bc_errors("PATCH      2 2", "PATCH      2 2",
                                        patches=[1, 1, 0, 1, 1, 1]),
                         "only 5 patches", "removes the 1 without faces")

    def test_periodic_pairs(self):
        self.assertError(self.bc_errors("RELPATCH 4 4", "RELPATCH 5 5"),
                         "patch 5 must point back at patch 3")
        self.assertError(self.bc_errors("RELPATCH 4 4\nAXIS     2", "RELPATCH 4 4\nAXIS     3"),
                         "same axis")
        self.assertError(self.bc_errors("RELPATCH 4 4\nAXIS     2\nANGLE    0.0",
                                        "RELPATCH 4 4\nAXIS     2\nANGLE    5.0"), "add up to 0")
        self.assertError(self.bc_errors("RELPATCH 4 4", "RELPATCH 2 2"),
                         "RELPATCH 2", "not a # BC_PERIODIC patch")

    def test_boundary_gradients_only_at_first_order(self):
        # NSCBC outflow (KIND 1): ORDER 2, or ORDERBF 2 without ORDER, crashes Rocflu
        self.assertError(self.bc_errors("ORDER      1", "ORDER      2"), "ORDER 2", "KIND 1")
        numerics = TEMPLATE.replace("# NUMERICS", "# NUMERICS\nORDERBF 2", 1)
        self.assertError(self.bc_errors("ORDER      1\n", "", text=numerics),
                         "# BC_OUTFLOW", "ORDERBF 2")
        self.assertEqual(self.bc_errors("ORDER      1", "ORDER      1", text=numerics), [])
        # ... unless boundary stencils are 1D
        sdim = numerics.replace("ORDERBF 2", "ORDERBF 2\nSDIMENSBF 1")
        self.assertEqual(self.bc_errors("ORDER      1", "ORDER      2", text=sdim), [])
        # Navier-Stokes needs them everywhere: the inflow's ORDER 2 too
        self.assertError(self.errors({("FLOWMODEL", "MODEL"): 1}), "ORDER 2", "Navier-Stokes")


if __name__ == "__main__":
    unittest.main(verbosity=1)
