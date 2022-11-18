"""
Module to store the mesh classes.

Currently, the only mesh format is h5.
"""

import os
import dolfin as df


class H5Mesh:
    """
    Mesh class.

    """

    def __init__(self, h5_file):
        """
        Construct the object.

        This includes:

            (#) Read mesh
            (#) Read subdomains
            (#) Read boudnaries

        """
        if not os.path.isfile(h5_file):
            raise Exception(f"{h5_file} not found")

        self.mesh = df.Mesh()

        hdf = df.HDF5File(self.mesh.mpi_comm(), h5_file, "r")

        hdf.read(self.mesh, "/mesh", False)
        dim = self.mesh.topology().dim()

        self.subdomains = df.MeshFunction("size_t", self.mesh, dim)
        hdf.read(self.subdomains, "/subdomains")

        self.boundaries = df.MeshFunction("size_t", self.mesh, dim - 1)
        hdf.read(self.boundaries, "/boundaries")
