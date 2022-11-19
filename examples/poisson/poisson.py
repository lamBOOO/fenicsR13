"""
Program to solve Poisson example.

"""

import dolfin as df
from input import Input
from meshes import H5Mesh
from solver import Solver


def print_information():
    print(r"""    Welcome to a simple Poisson Solver
    -> Contact: Manuel Torrilhon <mt@acom.rwth-aachen.de>
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
