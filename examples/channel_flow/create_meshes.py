#!/usr/bin/env python3

"Script to generate a set of ring meshes"

import os
import dolfin as df

GMSH_PATH = "gmsh"
GEO_NAME = "channel"

def create_mesh(exponent, remove_intermediate_files=True):
    "Generates a mesh using gmsh"

    mesh_name = "{}{}".format(GEO_NAME, exponent)

    os.system(
        "{} -setnumber p {} -2 -o {}.msh {}.geo".format(
            GMSH_PATH, exponent, mesh_name, GEO_NAME))

    os.system("dolfin-convert {0}.msh {0}.xml".format(mesh_name))

    mesh = df.Mesh("{}.xml".format(mesh_name))
    subdomains = df.MeshFunction(
        "size_t", mesh, "{}_physical_region.xml".format(mesh_name))
    boundaries = df.MeshFunction(
        "size_t", mesh, "{}_facet_region.xml".format(mesh_name))

    file = df.HDF5File(mesh.mpi_comm(), "{}.h5".format(mesh_name), "w")
    file.write(mesh, "/mesh")
    file.write(subdomains, "/subdomains")
    file.write(boundaries, "/boundaries")

    if remove_intermediate_files:
        # Remove intermediate files from..
        # 1) Gmsh:
        os.remove("{}.msh".format(mesh_name))
        # 2) dolfin-convert:
        os.remove("{}.xml".format(mesh_name))
        os.remove("{}_facet_region.xml".format(mesh_name))
        os.remove("{}_physical_region.xml".format(mesh_name))

    return (mesh, subdomains, boundaries)

for p in range(5, 5+1):
    create_mesh(p)
