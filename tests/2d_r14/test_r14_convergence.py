"""
Module to gather tests for convergence of the full R14 system.

This file is executed by ``pytest`` to have good CI.

Only the R14 cases live here. The ten r13 cases this package used to run
were duplicates: the seven input files it carried and all ten reference
error tables are byte-identical to the ones in ``tests/2d_r13``, and three
of the ten inputs were missing here altogether, which is why ``working_dir``
used to point at the neighbouring directory. Running them a second time
costs CI time and adds no coverage.

The R14 cases have no reference error tables yet; ``referrors/`` holds only
the copied r13 tables. The solver run is therefore asserted here, which is
what gives the R14 path its coverage, while the comparison against a
reference is reported as an expected failure until the tables exist. See
``test_exact_solution.py`` for the related gap in the exact solutions.
"""

import os
import subprocess
import sys
import pytest


class TestR14Convergence(object):
    """
    Class to bundle all R14 convergence tests.

    All tests are compared against reference errors.
    """

    working_dir = "tests/2d_r14"
    solver_path = "fenicsR13"

    def run_solver(self, inputfile):
        """
        Run the solver as subprocess with the given input file.

        The solver is invoked as a module with the repository root on
        ``PYTHONPATH`` instead of through the ``fenicsR13`` console script.
        The console script runs the copy of the package that is installed in
        the image, which does not contain the r14 formulation; invoking the
        module makes this test exercise the sources in this working tree.

        Test fails if subprocess return Exception or error.
        """
        repo_root = os.path.abspath(
            os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
        )
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(
            [repo_root] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])
        )
        subprocess.check_call([
            sys.executable, "-m", "fenicsR13.fenicsR13", inputfile
        ], cwd=self.working_dir, env=env)

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

    def compare_errors_if_reference_exists(self, name):
        """
        Compare against the reference table, or report the missing table.

        ``referrors/`` contains only the ten r13 tables copied from
        ``tests/2d_r13``; no R14 case has one. Calling ``compare_errors``
        would then fail on a missing file, which reads like a broken test
        rather than a missing reference. The run itself has already been
        asserted by the caller, so the gap is reported as an expected
        failure and turns into a visible pass once a table is added.
        """
        produced = os.path.join(name, "errors.csv")
        assert os.path.exists(os.path.join(self.working_dir, produced)), (
            "the run produced no error table at " + produced
        )
        reference = os.path.join("referrors", name, "errors.csv")
        if not os.path.exists(os.path.join(self.working_dir, reference)):
            pytest.xfail(
                "no reference error table for " + name + "; referrors/ "
                "contains only the r13 tables copied from tests/2d_r13"
            )
        self.compare_errors(name + "/errors.csv", reference)

    def test_r14_1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab(self):
        r"""
        Execute full linear R14 system test and check with reference errors.

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
        Elements             :math:`P_1P_1P_1P_1P_1P_1`
        Stabilization        CIP
        ==================== ===================================================
        """
        name = "r14_1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab"
        self.run_solver("inputs/" + name + ".yml")
        self.compare_errors_if_reference_exists(name)

    def test_r14_1_coeffs_nosources_norot_inflow_p2p2p2p2p2_stab(self):
        r"""
        Execute full linear R14 system test and check with reference errors.

        Same case as the :math:`P_1` variant with second order elements.

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
        Elements             :math:`P_2P_2P_2P_2P_2P_2`
        Stabilization        CIP
        ==================== ===================================================
        """
        name = "r14_1_coeffs_nosources_norot_inflow_p2p2p2p2p2_stab"
        self.run_solver("inputs/" + name + ".yml")
        self.compare_errors_if_reference_exists(name)

    @pytest.mark.xfail(
        strict=True,
        reason="the exact solution this case points at is a byte copy "
               "of the R13 one: no ScriptR class and no sqrt(2/7) "
               "layer, so the error calculation aborts. The R14 "
               "solution for this parameter set is not in the "
               "repository.",
    )
    def test_r14_1_coeffs_sources_rot_noinflow_p1p1p1p1p1_gls(self):
        r"""
        Execute full linear R14 system test and check with reference errors.

        Case with heat and mass sources, rotating inner wall and no inflow.

        ==================== ===================================================
        Parameter     Value
        ==================== ===================================================
        :math:`Kn`           :math:`1.0`
        :math:`\dot{m}`      :math:`1.0 (1 - \frac{5 R^2}{18}) \cos(\phi)`
        :math:`r`            :math:`1.0 (1 - \frac{5 R^2}{18}) \cos(\phi)`
        :math:`\theta_w^1`   :math:`1.0`
        :math:`v_t^1`        :math:`-10`
        :math:`v_n^1`        :math:`0`
        :math:`p_w^1`        :math:`0`
        :math:`\epsilon_w^1` :math:`0`
        :math:`\theta_w^2`   :math:`0.5`
        :math:`v_t^2`        :math:`0`
        :math:`v_n^2`        :math:`0`
        :math:`p_w^2`        :math:`0`
        :math:`\epsilon_w^2` :math:`0`
        Elements             :math:`P_1P_1P_1P_1P_1P_1`
        Stabilization        GLS
        ==================== ===================================================
        """
        name = "r14_1_coeffs_sources_rot_noinflow_p1p1p1p1p1_gls"
        self.run_solver("inputs/" + name + ".yml")
        self.compare_errors_if_reference_exists(name)
