"""
Module to gather tests for convergence of decoupled stress system.

This file is executed by ``pytest`` to have good CI.
"""

import subprocess
import pytest

class TestR13Convergence(object):
    """
    Class to bundle all stress convergence tests.

    All tests are compared against reference errors.
    """

    working_dir = "tests/r13"
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

    def test_r13_1_coeffs_sources_rot_noinflow_p1p1p1p1p1_stab(self):
        r"""
        Execute full linear R13 system test and check with reference errors.

        Test case is similar to  [WES2019]_.

        .. [WES2019] A. Westerkamp and M. Torrilhon. “Finite Element Methods
            for the Linear Regularized 13-Moment Equations Describing Slow
            Rarefied Gas Flows”. In: Journal of Computational Physics 389
            (2019).

        ==================== ===================================================
        Parameter     Value
        ==================== ===================================================
        :math:`Kn`           :math:`1.0`
        :math:`\dot{m}`      :math:`(1-\frac{5R^2}{18{Kn}^2})\cos(\phi)`
        :math:`r`            :math:`0`
        :math:`\theta_w^1`   :math:`1.0`
        :math:`v_t^1`        :math:`10.0`
        :math:`v_n^1`        :math:`0`
        :math:`p_w^1`        :math:`0`
        :math:`\epsilon_w^1` :math:`0`
        :math:`\theta_w^2`   :math:`0.5`
        :math:`v_t^2`        :math:`0.0`
        :math:`v_n^2`        :math:`0`
        :math:`p_w^2`        :math:`0`
        :math:`\epsilon_w^2` :math:`0`
        Elements             :math:`P_1P_1P_1P_1P_1`
        Stabilization        CIP: :math:`\delta_\theta,\delta_u=1,\delta_p=0.01`
        ==================== ===================================================
        """
        name = "r13_1_coeffs_sources_rot_noinflow_p1p1p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_r13_1_coeffs_sources_rot_noinflow_p2p2p2p2p2_stab(self):
        r"""
        Execute full linear R13 system test and check with reference errors.

        Test case is similar to  [WES2019]_.

        ==================== ===================================================
        Parameter     Value
        ==================== ===================================================
        :math:`Kn`           :math:`1.0`
        :math:`\dot{m}`      :math:`(1-\frac{5R^2}{18{Kn}^2})\cos(\phi)`
        :math:`r`            :math:`0`
        :math:`\theta_w^1`   :math:`1.0`
        :math:`v_t^1`        :math:`10.0`
        :math:`v_n^1`        :math:`0`
        :math:`p_w^1`        :math:`0`
        :math:`\epsilon_w^1` :math:`0`
        :math:`\theta_w^2`   :math:`0.5`
        :math:`v_t^2`        :math:`0.0`
        :math:`v_n^2`        :math:`0`
        :math:`p_w^2`        :math:`0`
        :math:`\epsilon_w^2` :math:`0`
        Elements             :math:`P_2P_2P_2P_2P_2`
        Stabilization        CIP: :math:`\delta_\theta,\delta_u=1,\delta_p=0.01`
        ==================== ===================================================
        """
        name = "r13_1_coeffs_sources_rot_noinflow_p2p2p2p2p2_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_r13_1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab(self):
        r"""
        Execute full linear R13 system test and check with reference errors.

        Test case is similar to  [TOR2017]_.

        .. [TOR2017] Torrilhon, M. et al. (2017). “Hierarchical Boltzmann
            simulations and model error estimation”. In: Journal of
            Computational Physics 342 (2017), pp. 66–84.

        ==================== ===================================================
        Parameter     Value
        ==================== ===================================================
        :math:`Kn`           :math:`1.0`
        :math:`\dot{m}`      :math:`0`
        :math:`r`            :math:`0`
        :math:`\theta_w^1`   :math:`1.0`
        :math:`v_t^1`        :math:`0`
        :math:`v_n^1`        :math:`0`
        :math:`p_w^1`        :math:`0`
        :math:`\epsilon_w^1` :math:`10^{-3}`
        :math:`\theta_w^2`   :math:`2.0`
        :math:`v_t^2`        :math:`-1.00 \sin(\phi)`
        :math:`v_n^2`        :math:`+1.00 \cos(\phi)`
        :math:`p_w^2`        :math:`-0.27 \cos(\phi)`
        :math:`\epsilon_w^2` :math:`10^{3}`
        Elements             :math:`P_1P_1P_1P_1P_1`
        Stabilization        CIP: :math:`\delta_\theta,\delta_u=1,\delta_p=0.01`
        ==================== ===================================================
        """
        name = "r13_1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)
