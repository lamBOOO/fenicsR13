# MPI TEST
# comm = d.MPI.comm_world
# rank = d.MPI.rank(comm)
# gmsh_path = "/Applications/gmsh/Gmsh.app/Contents/MacOS/gmsh"
# mesh_name = "ring"
# if rank==0:
# print("genmesh ************************ ")
# os.system("{} -setnumber p {} -2 {}.geo".format(gmsh_path, p,
#                                                 mesh_name))
# os.system("dolfin-convert {0}.msh {0}.xml".format(mesh_name))
# else:
#     print("wait ************************")