#!/usr/bin/env python3

# pylint: disable=invalid-name
# pylint: disable=unsubscriptable-object

"""
Program to solve linearized R13 equations.

Different modes allow to solved the ecoupled heat or stress system.
Different meshes can be used to perform a convergence study with given
exact solution.
"""

import sys
import gc
import dolfin as df

import meshes
from input import Input
from solver import Solver
from postprocessor import Postprocessor

def print_information():
    r"""
    Print program name and information.

    That is:

    .. code-block:: bash

        -> Version: v0.4
        -> Contact: Lambert Theisen <lambert.theisen@rwth-aachen.de>
        -> Contact: Prof. Dr. Manuel Torrilhon <mt@mathcces.rwth-aachen.de>
        -> Website: https://git.rwth-aachen.de/lamBOO/fenicsR13
          __            _          ____  _ _____
         / _| ___ _ __ (_) ___ ___|  _ \/ |___ /
        | |_ / _ \ '_ \| |/ __/ __| |_) | | |_ \
        |  _|  __/ | | | | (__\__ \  _ <| |___) |
        |_|  \___|_| |_|_|\___|___/_| \_\_|____/
    """
    print(r"""-> Version: v0.4
-> Contact: Lambert Theisen <lambert.theisen@rwth-aachen.de>
-> Contact: Prof. Dr. Manuel Torrilhon <mt@mathcces.rwth-aachen.de>
-> Website: https://git.rwth-aachen.de/lamBOO/fenicsR13
  __            _          ____  _ _____
 / _| ___ _ __ (_) ___ ___|  _ \/ |___ /
| |_ / _ \ '_ \| |/ __/ __| |_) | | |_ \
|  _|  __/ | | | | (__\__ \  _ <| |___) |
|_|  \___|_| |_|_|\___|___/_| \_\_|____/
""")

def main():
    """
    Execute the main program.

    Searches for an ``"input.yml"`` file in the current directory to use as
    input.

    Usage:

    .. code-block:: bash

        # Usage: <path_to_program> <input_file>
        # Goto case folder:
        cd tests/heat
        python3 ../../src.fenicsR13.py heat_01_coeffs_p1p1_stab.yml

    """
    print_information()

    # Dolfin settings
    df.set_log_level(100) # 1: all logs
    df.parameters["ghost_mode"] = "shared_vertex"

    inputfile = sys.argv[1] if len(sys.argv) == 2 else "input.yml"

    params = Input(inputfile).dict
    mesh_names = params["meshes"]

    convergence_study = params["convergence_study"]["enable"]
    show_plot = params["convergence_study"]["plot"]

    data = []

    for p, mesh_name in enumerate(mesh_names):

        print("Mesh: " + mesh_name)

        mesh_name = mesh_names[p]

        current_mesh = meshes.H5Mesh(mesh_name)
        solver = Solver(params, current_mesh, p)

        solver.assemble()
        solver.solve()
        solver.write()

        if convergence_study:

            errors = solver.calculate_errors()

            data.append({
                "h": current_mesh.mesh.hmax(),
                **errors
            })

            if p == len(mesh_names)-1: # after last mesh
                postp = Postprocessor(data, params["case_name"])
                postp.write_errors()
                postp.plot_errors(show_plot)

        solver = None
        gc.collect()

if __name__ == '__main__':
    main()
