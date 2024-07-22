#!/usr/bin/env python3

# pylint: disable=invalid-name

"""
Converter from geo-format to a mesh in h5-format.

Installation:

.. code-block:: bash

    pip install .

Usage:

.. code-block:: bash

    Usage: python3 <path_to_geoToH5.py> <geo_file> <h5_file> \
            [<gmsh cli arguments>]
    E.g.: geoToH5 lid.geo lid5.h5 "-setnumber p 5"

3D Meshes:
Append "-3" to the gmsh CLI arguments like:

.. code-block:: bash

    geoToH5 mesh3d.geo mesh3d.h5 "-3"
"""

import os
import sys
import dolfin as df
from datetime import datetime
import random

# Constants
GMSH_PATH = "gmsh"
GEO_NAME = "ring"


def geo_to_h5():
    "Convert given geo-file to a h5-mesh."

    # Check if right arguments passed, print usage if not
    if not 3 <= len(sys.argv) <= 4:
        print("""
Usage: python3 <path_to_geoToH5.py> <geo_file> <h5_file> [<gmsh cli arguments>]
E.g.: geoToH5 lid.geo lid5.h5 "-setnumber p 5"

3D meshes: add "-3" to the gmsh CLI arguments
        """)
        return

    geo_input_file = sys.argv[1]
    h5_output_file = sys.argv[2]
    gmsh_arguments = sys.argv[3] if len(sys.argv) == 4 else ""
    tmp_name = "tmp_" + str(datetime.now().time()) + str(random.random())

    # Create msh-mesh with Gmsh
    os.system(
        "{} -2 {} -o {}.msh {}".format(
            GMSH_PATH, gmsh_arguments, tmp_name, geo_input_file
        )
    )

    # Convert msh-mesh to xml-mesh
    os.system("dolfin-convert {0}.msh {0}.xml".format(tmp_name))

    # Delete msh-mesh
    os.remove("{}.msh".format(tmp_name))

    # Read xml-mesh
    mesh = df.Mesh("{}.xml".format(tmp_name))
    subdomains = df.MeshFunction(
        "size_t", mesh, "{}_physical_region.xml".format(tmp_name))
    boundaries = df.MeshFunction(
        "size_t", mesh, "{}_facet_region.xml".format(tmp_name))

    # Delete xml-mesh
    os.remove("{}.xml".format(tmp_name))
    os.remove("{}_physical_region.xml".format(tmp_name))
    os.remove("{}_facet_region.xml".format(tmp_name))

    # Write h5-mesh
    file = df.HDF5File(mesh.mpi_comm(), h5_output_file, "w")
    file.write(mesh, "/mesh")
    file.write(subdomains, "/subdomains")
    file.write(boundaries, "/boundaries")
