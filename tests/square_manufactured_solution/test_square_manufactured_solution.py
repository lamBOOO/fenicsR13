"""
This file is executed by ``pytest``.
"""

import subprocess


class TestR13Convergence(object):
    """
    Class to bundle all stress convergence tests.

    All tests are compared against reference errors.
    """

    working_dir = "tests/square_manufactured_solution"
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

        Absolute Error allowed: ``1E-6``
        Return exception if diff returns with !=0
        A comparison for complete equalness can be obtained with:

        .. code-block:: python

            subprocess.check_call([
                "diff", "-u", "--strip-trailing-cr", errorsfile, ref_errorsfile
            ], cwd=self.working_dir)

        """
        print(subprocess.check_output([
            "numdiff", "-s", "\"\n\r ,\"", "-a", "1E-4",
            errorsfile, ref_errorsfile
        ], cwd=self.working_dir))

    def test_2d_heat(self):
        r"""
        Test the 2D Heat convergence of the manufactured solution.
        """
        name = "square_manufactured_solution_2d_heat"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_2d_stress(self):
        r"""
        Test the 2D Stress convergence of the manufactured solution.
        """
        name = "square_manufactured_solution_2d_stress"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_2d_r13(self):
        r"""
        Test the 2D R13 convergence of the manufactured solution.
        """
        name = "square_manufactured_solution_2d_r13"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_3d_heat(self):
        r"""
        Test the 3D Heat convergence of the manufactured solution.
        """
        name = "square_manufactured_solution_3d_heat"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_3d_stress(self):
        r"""
        Test the 3D Stress convergence of the manufactured solution.
        """
        name = "square_manufactured_solution_3d_stress"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)

    def test_3d_r13(self):
        r"""
        Test the 3D R13 convergence of the manufactured solution.
        """
        name = "square_manufactured_solution_3d_r13"
        self.run_solver("inputs/" + name + ".yml")
        errors = name + "/" + "errors.csv"
        referrors = "referrors/" + name + "/errors.csv"
        self.compare_errors(errors, referrors)
