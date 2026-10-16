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
from rocpicl_tests.case import edit_inp        # noqa: E402

with open(os.path.join(TESTS_DIR, "templates", "rectshktb.inp")) as f:
    TEMPLATE = f.read()


class CaseValidatorTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir)

    def errors(self, overrides=None, text=None, files=("points.dat",)):
        """Errors for the template with {(SECTION, KEY): value} changes, or
        for the given input text, with the given files in the case directory."""
        if text is None:
            text = edit_inp(TEMPLATE, overrides or {})
        inp = os.path.join(self.dir, "case.inp")
        with open(inp, "w") as f:
            f.write(text)
        for name in files:
            open(os.path.join(self.dir, name), "w").close()
        return case_validator.validate(inp)

    def assertError(self, errors, *words):
        """One of the errors contains all the words."""
        self.assertTrue(any(all(w in e for w in words) for e in errors),
                        "no error with %s in:\n  %s" % (words, "\n  ".join(errors)))

    # -- accepted cases ----------------------------------------------------------

    def test_template_is_fine(self):
        self.assertEqual(self.errors(), [])

    def test_pseudo_turbulence_with_what_it_needs(self):
        self.assertEqual(self.errors({("PICL", "PSEUDOTURB"): 1, ("PICL", "QSFLUCT"): 2,
                                      ("PICL", "HEATTRANSFER"): 2}), [])

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
                           ("ANGULARPERIODIC 0", ("angular",))):
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

    # -- combinations ----------------------------------------------------------------

    def test_fluctuations_need_quasi_steady_force(self):
        self.assertError(self.errors({("PICL", "QSFORCE"): 0, ("PICL", "QSFLUCT"): 1}),
                         "QSFLUCT 1 needs a quasi-steady force")

    def test_pseudo_turbulence_needs(self):
        errors = self.errors({("PICL", "PSEUDOTURB"): 1, ("PICL", "QSFLUCT"): 1})
        self.assertError(errors, "PSEUDOTURB 1 needs QSFLUCT 2")
        self.assertError(errors, "PSEUDOTURB 1 needs HEATTRANSFER")
        self.assertError(self.errors({("PICL", "PSEUDOTURB"): 1, ("PICL", "QSFLUCT"): 2,
                                      ("PICL", "HEATTRANSFER"): 2, ("PICL", "QSFORCE"): 0}),
                         "PSEUDOTURB 1 needs a quasi-steady force")

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

    def test_files_in_case_directory(self):
        self.assertError(self.errors(files=()), "no points.dat")
        self.assertError(self.errors(files=("points.dat", "filein.vtk")), "filein.vtk")
        self.assertError(self.errors(files=("points.dat", "par00002.vtu")), "restart")


if __name__ == "__main__":
    unittest.main(verbosity=1)
