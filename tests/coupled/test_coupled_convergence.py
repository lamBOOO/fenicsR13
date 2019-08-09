"""
Module to gather tests for convergence of decoupled coupled system.
"""

import subprocess
import pytest

class TestCoupledConvergence(object):
    """
    Class to bundle all coupled convergence tests.
    All tests are compared against reference errors.
    """

    working_dir = "tests/coupled"
    solver_path = "../../src/fenicsR13.py"

    def run_solver(self, inputfile):
        """
        Runs the solver as subprocess with the given input file.
        """
        subprocess.check_call([
            "python3", self.solver_path, inputfile
        ], cwd=self.working_dir)

    def compare_errors(self, errorsfile, ref_errorsfile):
        """
        Check against reference errors. Compares absolute differences.
        Absolute Error allowed: ``1E-10``
        Return exception if diff returns with !=0
        A comparison for complete equalness can be obtained with:

        .. code-block:: python

            subprocess.check_call([
                "diff", "-u", "--strip-trailing-cr", errorsfile, ref_errorsfile
            ], cwd=self.working_dir)
        """
        subprocess.check_call([
            "numdiff", "-s", "\"\n\r ,\"", "-a", "1E-10",
            errorsfile, ref_errorsfile
        ], cwd=self.working_dir)

    # @pytest.fixture(scope="module", autouse=True)
    @pytest.mark.skip(reason="Not needed because meshes are in repo")
    def create_meshes(self):
        """
        Creates the test meshes. Executed before any test of the class.
        """
        subprocess.check_call(["python3", "create_meshes.py"], cwd="tests/mesh")

    def test_1_coeffs_sources_rot_p1p1p1p1p1_stab(self):
        r"""
        Executes westerkamp2019 coupled system test and check with
        reference errors.

        ========================= ==============================================
        Parameter     Value
        ========================= ==============================================
        :math:`\tau`              :math:`1.0`
        :math:`f_{\mathrm{mass}}` :math:`(1-\frac{5R^2}{18\tau^2})\cos(\phi)`
        :math:`f_{\mathrm{heat}}` :math:`0`
        :math:`v_t^1`             :math:`10.0`
        :math:`v_t^2`             :math:`0.0`
        Elements                  :math:`P_1P_1P_1P_1P_1`
        Stabilization             CIP: :math:`\delta_1,\delta_2=1,\delta_3=0.01`
        ========================= ==============================================
        """
        name = "1_coeffs_sources_rot_p1p1p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_1_coeffs_sources_rot_p2p2p2p2p2_stab(self):
        r"""
        Executes westerkamp2019 coupled system test and check with
        reference errors.

        ========================= ==============================================
        Parameter     Value
        ========================= ==============================================
        :math:`\tau`              :math:`1.0`
        :math:`f_{\mathrm{mass}}` :math:`(1-\frac{5R^2}{18\tau^2})\cos(\phi)`
        :math:`f_{\mathrm{heat}}` :math:`0`
        :math:`v_t^1`             :math:`10.0`
        :math:`v_t^2`             :math:`0.0`
        Elements                  :math:`P_2P_2P_2P_2P_2`
        Stabilization             CIP: :math:`\delta_1,\delta_2=1,\delta_3=0.01`
        ========================= ==============================================
        """
        name = "1_coeffs_sources_rot_p2p2p2p2p2_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)
