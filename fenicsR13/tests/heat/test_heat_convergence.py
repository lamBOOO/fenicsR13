"Test for convergence"

import subprocess

def test_heat_convergence():
    """
    Executes westerkamp2019 decoupled heat system test and check with
    reference errors
    """

    subprocess.call(["python3", "create_meshes.py"], cwd="tests/mesh")

    subprocess.call(["python3", "../../src/fenicsR13.py"], cwd="tests/heat")

    # Check against reference errors, return exception if diff returns with !=0
    subprocess.check_call([
        "diff", "-u", "--strip-trailing-cr", "errors.csv", "errors_ref.csv"
    ], cwd="tests/heat")
