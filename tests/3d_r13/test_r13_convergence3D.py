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

    working_dir = "tests/3d_r13"
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
            "numdiff", "-s", "\"\n\r ,\"", "-a", "1E-2",
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

    def test_r13_kn0_2(self):
        name = "r13_kn0.2"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_r13_kn0_005(self):
        name = "r13_kn0.005"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_r13_kn0_005_shell(self):
        name = "r13_kn0.005_shell"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)
