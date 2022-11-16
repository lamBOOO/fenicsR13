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


def print_information():
    r"""
    Print program name and information.
    """
    print(r"""-> Version: 1.4
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
    df.set_log_level(100)  # 1: all logs
    df.parameters["ghost_mode"] = "shared_vertex"

    inputfile = "input.yml"
    input_file = Input(inputfile)
    params = input_file.dict

    mesh_names = params["meshes"]
    mesh_name = mesh_names[0]
    print("Mesh: " + mesh_name)

    current_mesh = H5Mesh(mesh_name)
    solver = Solver(params, current_mesh, 0)

    solver.assemble()
    solver.solve()
    solver.write()

    solver = None
    gc.collect()



if __name__ == '__main__':
    main()
