import meshio
import dolfin

msh = meshio.read("ring0.msh")
meshio.write("mesh.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}))




# meshio.write("mesh.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}))
# meshio.write("mf.xdmf", meshio.Mesh(points=msh.points, cells={"line": msh.cells["line"]},
# cell_data={"line": {"name_to_read": msh.cell_data["line"]["gmsh:physical"]}}))

# mesh = Mesh()
# with XDMFFile("msh.xdmf") as infile:
#     infile.read(mesh)
#     mvc = MeshValueCollection("size_t", mesh, 2)
# with XDMFFile("mf.xdmf") as infile:
#     infile.read(mvc, "name_to_read")
#     mf = dolfin.cpp.mesh.MeshFunctionSizet(mesh, mvc)
#     File("dolfinmesh.pvd").write(mesh)
#     File("dolfincellfunc.pvd").write(mf)