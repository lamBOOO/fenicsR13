"""
Program to solve linearized R13 equations.

"""

import dolfin as df
from input import Input
from meshes import H5Mesh
from solver import Solver


def print_information():
    print(r"""    Welcome to
      __            _          ____  _ _____
     / _| ___ _ __ (_) ___ ___|  _ \/ |___ /
    | |_ / _ \ '_ \| |/ __/ __| |_) | | |_ \
    |  _|  __/ | | | | (__\__ \  _ <| |___) |
    |_|  \___|_| |_|_|\___|___/_| \_\_|____/
    -> Version: 1.4
    -> Contact: Lambert Theisen <lambert.theisen@rwth-aachen.de>
    -> Contact: Manuel Torrilhon <mt@mathcces.rwth-aachen.de>
    -> Repository: <https://git.rwth-aachen.de/lamBOO/fenicsR13>
    -> Documentation: <https://lamboo.pages.rwth-aachen.de/fenicsR13/>
""")



# main program instructions
############################
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
