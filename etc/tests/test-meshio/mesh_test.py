# import meshio

# geometry = meshio.read("ring0.msh")
# meshio.write("mesh.xdmf", meshio.Mesh(points=geometry.points, cells={"triangle": geometry.cells["triangle"]}))
# meshio.write("mf.xdmf", meshio.Mesh(points=geometry.points, cells={"line": geometry.cells["line"]},cell_data={"line": {"name_to_read": geometry.cell_data["line"]["gmsh:physical"]}}))


import meshio
msh = meshio.read("ring.msh")

meshio.write("mesh.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}))
meshio.write("mf.xdmf", meshio.Mesh(points=msh.points, cells={"line": msh.cells["line"]},
                                    cell_data={"line": {"name_to_read": msh.cell_data["line"]["gmsh:physical"]}}))
meshio.write("cf.xdmf", meshio.Mesh(
    points=msh.points, cells={"line": msh.cells["line"]},
    cell_data={"line": {"name_to_read":
                            msh.cell_data["line"]["gmsh:physical"]}}))


from dolfin import *
mesh = Mesh()
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 2)
with XDMFFile("mf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

mvc = MeshValueCollection("size_t", mesh, 3)
with XDMFFile("cf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
cf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

ds_custom = Measure("ds", domain=mesh, subdomain_data=mf, subdomain_id=3100)
print(assemble(1*ds_custom))