"Test for convergence"

import subprocess
import pytest


class TestHeatConvergence(object):
    """
    Class to bundle all heat convergence tests.
    """

    cwd = "tests/heat"
    solver = "../../src/fenicsR13.py"

    def run_solver(self, inputfile):
        """
        Runs the solver as subprocess with the given input file.
        """
        subprocess.check_call(["python3", self.solver, inputfile], cwd=self.cwd)

    def compare_errors(self, errorsfile, ref_errorsfile):
        """
        Check against reference errors.
        Return exception if diff returns with !=0
        """
        subprocess.check_call([
            "diff", "-u", "--strip-trailing-cr", errorsfile, ref_errorsfile
        ], cwd=self.cwd)

    @pytest.fixture(scope="module", autouse=True)
    def create_meshes(self):
        """
        Creates the test meshes. Executed before any test of the class.
        """
        subprocess.check_call(["python3", "create_meshes.py"], cwd="tests/mesh")

    def test_coeffs_p1p1_stab(self):
        """
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.
        Use P1P1 elements with CIP stabilization.
        """
        self.run_solver("inputs/coeffs_p1p1_stab.yml")
        self.compare_errors("errors.csv", "referrors/coeffs_p1p1_stab.csv")

    def test_coeffs_p2p2_stab(self):
        """
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.
        Use P2P2 elements with CIP stabilization.
        """
        self.run_solver("inputs/coeffs_p2p2_stab.yml")
        self.compare_errors("errors.csv", "referrors/coeffs_p2p2_stab.csv")

    @pytest.mark.skip(reason="Not implemented")
    def test_coeffs_p2p1_nostab(self):
        """
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.
        Use P2P1 elements without stabilization.
        """
        self.run_solver("inputs/coeffs_p1p2_nostab.yml")
        self.compare_errors("errors.csv", "referrors/coeffs_p1p2_nostab.csv")

# TestHeatConvergence().test_coeffs_p2p1_nostab()
TestHeatConvergence().test_coeffs_p2p2_stab()
