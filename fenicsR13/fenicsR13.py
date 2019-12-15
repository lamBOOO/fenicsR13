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

from fenicsR13.meshes import H5Mesh
from fenicsR13.input import Input
from fenicsR13.solver import Solver
from fenicsR13.postprocessor import Postprocessor

def print_information():
    r"""
    Print program name and information.

    That is:

    .. code-block:: text

        -> Version: v1.1
        -> Contact: Lambert Theisen <lambert.theisen@rwth-aachen.de>
        -> Contact: Manuel Torrilhon <mt@mathcces.rwth-aachen.de>
        -> Repository: <https://git.rwth-aachen.de/lamBOO/fenicsR13>
        -> Documentation: <https://lamboo.pages.rwth-aachen.de/fenicsR13/>
          __            _          ____  _ _____
         / _| ___ _ __ (_) ___ ___|  _ \/ |___ /
        | |_ / _ \ '_ \| |/ __/ __| |_) | | |_ \
        |  _|  __/ | | | | (__\__ \  _ <| |___) |
        |_|  \___|_| |_|_|\___|___/_| \_\_|____/

    """
    print(r"""-> Version: v1.1
-> Contact: Lambert Theisen <lambert.theisen@rwth-aachen.de>
-> Contact: Manuel Torrilhon <mt@mathcces.rwth-aachen.de>
-> Repository: <https://git.rwth-aachen.de/lamBOO/fenicsR13>
-> Documentation: <https://lamboo.pages.rwth-aachen.de/fenicsR13/>
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

        # Install fenicsR13
        pip install .

        # Usage: <path_to_program> <input_file>
        # Goto case folder:
        cd tests/r13
        fenicsR13 inputs/ \
          r13_1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab.yml

    """
    print_information()

    # Dolfin settings
    df.set_log_level(100) # 1: all logs
    df.parameters["ghost_mode"] = "shared_vertex"

    inputfile = sys.argv[1] if len(sys.argv) == 2 else "input.yml"

    input_file = Input(inputfile)
    params = input_file.dict

    # Setup parameter study loop
    if params["parameter_study"]["enable"]:
        parameter_values = params["parameter_study"]["parameter_values"]
        parameter_key = params["parameter_study"]["parameter_key"]
    else:
        parameter_values = [""]
        parameter_key = [""]
    initial_output_folder = params["output_folder"]
    for parameter_value in parameter_values:
        input_file.set_in_input(parameter_key, parameter_value)
        params["output_folder"] = initial_output_folder + str(parameter_value)
        print("Study", parameter_key, ":", parameter_value)

        # Ususal code:
        mesh_names = params["meshes"]

        convergence_study = params["convergence_study"]["enable"]
        show_plot = params["convergence_study"]["plot"]

        data = []

        for p, mesh_name in enumerate(mesh_names):

            print("Mesh: " + mesh_name)

            mesh_name = mesh_names[p]

            current_mesh = H5Mesh(mesh_name)
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
                    postp = Postprocessor(data, params["output_folder"])
                    postp.write_errors()
                    postp.plot_errors(show_plot)

            solver = None
            gc.collect()

if __name__ == '__main__':
    main()
