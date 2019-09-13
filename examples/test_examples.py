"""
Module to perform tests for all example cases.

This file is executed by ``pytest`` to have good CI.
"""

import subprocess
import pytest

class TestExamples(object):
    """
    Class to bundle all examples tests.

    All tests are compared against reference errors.
    """

    solver_path = "../../src/fenicsR13.py"

    def run_solver(self, inputfile, working_dir_):
        """
        Run the solver as subprocess with the given input file.

        Test fails if subprocess return Exception or error.
        """
        subprocess.check_call([
            "python3", self.solver_path, inputfile
        ], cwd=working_dir_)

    # @pytest.fixture(scope="module", autouse=True)
    @pytest.mark.skip(reason="Not needed because meshes are in repo")
    def create_meshes(self, working_dir_):
        """
        Create the test meshes. Executed before any test of the class.

        Often not needed if meshes are in Git through LFS for reproducability.
        """
        subprocess.check_call(["python3", "create_meshes.py"], cwd=working_dir_)

    def test_lid_driven_cavity(self):
        r"""
        Test the lid driven cavity case.
        """
        working_dir = "examples/lid_driven_cavity"
        self.create_meshes(working_dir)
        self.run_solver("input.yml", working_dir)

    def test_channel_flow(self):
        r"""
        Test the channel flow case and generate table with Kn vs. massflow.
        """
        working_dir = "examples/channel_flow"
        self.create_meshes(working_dir)
        self.run_solver("input.yml", working_dir)
        subprocess.check_call(["bash", "postprocessing.sh"], cwd=working_dir)

    def test_knudsen_pump(self):
        r"""
        Test the knudsen pump case.
        """
        working_dir = "examples/knudsen_pump"
        self.create_meshes(working_dir)
        self.run_solver("input.yml", working_dir)
