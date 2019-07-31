"Class to handle a mesh"

import os
import dolfin as df

class H5Mesh:
    "Mesh class"
    def __init__(self, h5_file):

        if not os.path.isfile(h5_file):
            raise Exception(f"[{self}]: {h5_file} not found")

        self.mesh = df.Mesh()

        hdf = df.HDF5File(self.mesh.mpi_comm(), h5_file, "r")

        hdf.read(self.mesh, "/mesh", False)
        dim = self.mesh.topology().dim()

        self.subdomains = df.MeshFunction("size_t", self.mesh, dim)
        hdf.read(self.subdomains, "/subdomains")

        self.boundaries = df.MeshFunction("size_t", self.mesh, dim - 1)
        hdf.read(self.boundaries, "/boundaries")
