"""
Program to solve Poisson example.

"""

import os
import dolfin as df
from input import Input
from meshes import H5Mesh
from solver import Solver


def print_information():
    print(r"""    Welcome to a simple Poisson Solver
    -> Contact: Manuel Torrilhon <mt@acom.rwth-aachen.de>
""")


def generate_h5_mesh( geofile, refinement ):
    print("Geometry: " + geofile)
    tmpfile = "tmpmesh"
    refinement = 1

    os.system( "gmsh -setnumber p {} -2 -o {}.msh {}".format(refinement,geofile,geofile) )
    os.system( "dolfin-convert {}.msh {}.xml".format(geofile,tmpfile))
    # Read xml-mesh
    mesh = df.Mesh("{}.xml".format(tmpfile))
    subdomains = df.MeshFunction("size_t", mesh, "{}_physical_region.xml".format(tmpfile))
    boundaries = df.MeshFunction("size_t", mesh, "{}_facet_region.xml".format(tmpfile))
    # Write h5-mesh
    file = df.HDF5File(mesh.mpi_comm(), geofile + ".h5", "w")
    file.write(mesh, "/mesh")
    file.write(subdomains, "/subdomains")
    file.write(boundaries, "/boundaries")
    # Delete xml-mesh
    #os.remove(geofile+".msh")
    os.remove(tmpfile+".xml")
    os.remove(tmpfile+"_physical_region.xml")
    os.remove(tmpfile+"_facet_region.xml")


# main program instructions
############################
print_information()

# Dolfin settings
df.set_log_level(100)  # 1: all logs
df.parameters["ghost_mode"] = "shared_vertex"

inputfile = "input.yml"
input_file = Input(inputfile)
params = input_file.dict

geofile = params["geometry"]
generate_h5_mesh( geofile, 2 )

mesh_name = geofile + ".h5"
print("Mesh: " + mesh_name)
current_mesh = H5Mesh(mesh_name)

solver = Solver(params, current_mesh, 0)

solver.assemble()
solver.solve()
solver.write()
