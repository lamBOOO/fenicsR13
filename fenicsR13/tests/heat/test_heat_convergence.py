"Test for convergence"

import subprocess
import pytest


class TestHeatConvergence(object):
    """
    Class to bundle all heat convergence tests.
    All tests are compared against reference errors.
    """

    working_dir = "tests/heat"
    solver_path = "../../src/fenicsR13.py"

    def run_solver(self, inputfile):
        """
        Runs the solver as subprocess with the given input file.
        """
        subprocess.check_call(["python3", self.solver_path, inputfile], cwd=self.working_dir)

    def compare_errors(self, errorsfile, ref_errorsfile):
        """
        Check against reference errors.
        Return exception if diff returns with !=0
        """
        subprocess.check_call([
            "diff", "-u", "--strip-trailing-cr", errorsfile, ref_errorsfile
        ], cwd=self.working_dir)

    @pytest.fixture(scope="module", autouse=True)
    def create_meshes(self):
        """
        Creates the test meshes. Executed before any test of the class.
        """
        subprocess.check_call(["python3", "create_meshes.py"], cwd="tests/mesh")

    def test_01_coeffs_p1p1_stab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ============= =======================
        Parameter     Value
        ============= =======================
        :math:`\tau`  :math:`0.1`
        Formulation   Coefficients
        Elements      :math:`P_1P_1`
        Stabilization CIP, :math:`\delta_1=1`
        ============= =======================
        """
        self.run_solver("inputs/01_coeffs_p1p1_stab.yml")
        self.compare_errors("errors.csv", "referrors/01_coeffs_p1p1_stab.csv")

    def test_10_coeffs_p2p2_stab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ============= =======================
        Parameter     Value
        ============= =======================
        :math:`\tau`  :math:`10.0`
        Formulation   Coefficients
        Elements      :math:`P_2P_2`
        Stabilization CIP, :math:`\delta_1=1`
        ============= =======================
        """
        self.run_solver("inputs/10_coeffs_p2p2_stab.yml")
        self.compare_errors("errors.csv", "referrors/10_coeffs_p2p2_stab.csv")

    def test_01_coeffs_p2p2_stab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ============= =======================
        Parameter     Value
        ============= =======================
        :math:`\tau`  :math:`0.1`
        Formulation   Coefficients
        Elements      :math:`P_2P_2`
        Stabilization CIP, :math:`\delta_1=1`
        ============= =======================

        Use P2P2 elements with CIP stabilization.
        """
        self.run_solver("inputs/01_coeffs_p2p2_stab.yml")
        self.compare_errors("errors.csv", "referrors/01_coeffs_p2p2_stab.csv")

    # @pytest.mark.skip(reason="Not implemented")
    def test_01_coeffs_p2p1_nostab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ============= =======================
        Parameter     Value
        ============= =======================
        :math:`\tau`  :math:`0.1`
        Formulation   Coefficients
        Elements      :math:`P_2P_1`
        Stabilization CIP, :math:`\delta_1=1`
        ============= =======================
        """
        self.run_solver("inputs/01_coeffs_p1p2_nostab.yml")
        self.compare_errors("errors.csv", "referrors/01_coeffs_p1p2_nostab.csv")

# TestHeatConvergence().test_coeffs_p2p1_nostab()
# TestHeatConvergence().test_10_coeffs_p2p2_stab()
