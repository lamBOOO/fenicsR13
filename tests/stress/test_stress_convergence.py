"""
Module to gather tests for convergence of decoupled stress system.
"""

import subprocess
import pytest

class TestStressConvergence(object):
    """
    Class to bundle all stress convergence tests.
    All tests are compared against reference errors.
    """

    working_dir = "tests/stress"
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
        print(subprocess.check_output([
            "numdiff", "-s", "\"\n\r ,\"", "-a", "1E-10",
            errorsfile, ref_errorsfile
        ], cwd=self.working_dir))

    # @pytest.fixture(scope="module", autouse=True)
    @pytest.mark.skip(reason="Not needed because meshes are in repo")
    def create_meshes(self):
        """
        Creates the test meshes. Executed before any test of the class.
        """
        subprocess.check_call(["python3", "create_meshes.py"], cwd="tests/mesh")

    def test_01_nosource_rot_p1p1p1_stab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`\tau`              :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_1P_1P_1`
        Stabilization             CIP, :math:`\delta_2=1, \delta_3=0.01`
        ========================= =============================================
        """
        name = "01_nosource_rot_p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_01_source_norot_p1p1p1_stab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`\tau`              :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18\tau^2})\cos(\phi)`
        :math:`v_t^1`             :math:`0`
        Elements                  :math:`P_1P_1P_1`
        Stabilization             CIP, :math:`\delta_2=1, \delta_3=0.01`
        ========================= =============================================
        """
        name = "01_source_norot_p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_01_source_rot_p1p1p1_stab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`\tau`              :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18\tau^2})\cos(\phi)`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_1P_1P_1`
        Stabilization             CIP, :math:`\delta_2=1, \delta_3=0.01`
        ========================= =============================================
        """
        name = "01_source_rot_p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_01_source_rot_p1p2p4_nostab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`\tau`              :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18\tau^2})\cos(\phi)`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_1P_2P_4`
        Stabilization             CIP, :math:`\delta_1=1`
        ========================= =============================================
        """
        name = "01_source_rot_p1p2p4_nostab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_01_source_rot_p2p2p2_stab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`\tau`              :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18\tau^2})\cos(\phi)`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_2P_2P_2`
        Stabilization             CIP, :math:`\delta_2=1, \delta_3=0.01`
        ========================= =============================================
        """
        name = "01_source_rot_p2p2p2_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_10_source_rot_p1p1p1_stab(self):
        r"""
        Executes westerkamp2019 decoupled heat system test and check with
        reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`\tau`              :math:`10.0`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18\tau^2})\cos(\phi)`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_1P_1P_1`
        Stabilization             CIP, :math:`\delta_2=1, \delta_3=0.01`
        ========================= =============================================
        """
        name = "10_source_rot_p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)
