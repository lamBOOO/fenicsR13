"""
Module to gather tests for convergence of decoupled stress system.

This file is executed by ``pytest`` to have good CI.
"""

import subprocess
import pytest

class TestStressConvergence(object):
    """
    Class to bundle all stress convergence tests.

    All tests are compared against reference errors.
    """

    working_dir = "tests/stress"
    solver_path = "fenicsR13"

    def run_solver(self, inputfile):
        """
        Run the solver as subprocess with the given input file.

        Test fails if subprocess return Exception or error.
        """
        subprocess.check_call([
            self.solver_path, inputfile
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
        Create the test meshes. Executed before any test of the class.

        Often not needed if meshes are in Git through LFS for reproducability.
        """
        subprocess.check_call(["python3", "create_meshes.py"], cwd="tests/mesh")

    def test_stress_01_nosource_rot_p1p1p1_stab(self):
        r"""
        Execute decoupled heat system test and check with reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`Kn`                :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_1P_1P_1`
        Stabilization             CIP, :math:`\delta_u=1, \delta_p=0.01`
        ========================= =============================================
        """
        name = "stress_01_nosource_rot_p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_stress_01_source_norot_p1p1p1_stab(self):
        r"""
        Execute decoupled heat system test and check with reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`Kn`                :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18{Kn}^2})\cos(\phi)`
        :math:`v_t^1`             :math:`0`
        Elements                  :math:`P_1P_1P_1`
        Stabilization             CIP, :math:`\delta_u=1, \delta_p=0.01`
        ========================= =============================================
        """
        name = "stress_01_source_norot_p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_stress_01_source_rot_p1p1p1_stab(self):
        r"""
        Execute decoupled heat system test and check with reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`Kn`                :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18{Kn}^2})\cos(\phi)`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_1P_1P_1`
        Stabilization             CIP, :math:`\delta_u=1, \delta_p=0.01`
        ========================= =============================================
        """
        name = "stress_01_source_rot_p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_stress_01_source_rot_p1p2p3_nostab(self):
        r"""
        Execute decoupled heat system test and check with reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`Kn`                :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18{Kn}^2})\cos(\phi)`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_1P_2P_3`
        Stabilization             CIP, :math:`\delta_u=1, \delta_p=0.01`
        ========================= =============================================
        """
        name = "stress_01_source_rot_p1p2p3_nostab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_stress_01_source_rot_p2p2p2_stab(self):
        r"""
        Execute decoupled heat system test and check with reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`Kn`                :math:`0.1`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18{Kn}^2})\cos(\phi)`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_2P_2P_2`
        Stabilization             CIP, :math:`\delta_u=1, \delta_p=0.01`
        ========================= =============================================
        """
        name = "stress_01_source_rot_p2p2p2_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_stress_10_source_rot_p1p1p1_stab(self):
        r"""
        Execute decoupled heat system test and check with reference errors.

        ========================= =============================================
        Parameter     Value
        ========================= =============================================
        :math:`Kn`                :math:`10.0`
        :math:`f_{\mathrm{mass}}` :math:`0.4(1-\frac{5R^2}{18{Kn}^2})\cos(\phi)`
        :math:`v_t^1`             :math:`10.0`
        Elements                  :math:`P_1P_1P_1`
        Stabilization             CIP, :math:`\delta_u=1, \delta_p=0.01`
        ========================= =============================================
        """
        name = "stress_10_source_rot_p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)
