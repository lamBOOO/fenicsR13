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



# ##### OLD, NOT WORKING WELL
# R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
# phi = df.Expression("atan2(x[1],x[0])", degree=5)

# x = df.Expression("x[0]", degree=1)
# y = df.Expression("x[1]", degree=1)

# C_0 = df.Constant(-50.80230139855979)
# C_1 = df.Constant(0.6015037593984962)
# C_2 = df.Constant(0)
# C_3 = df.Constant(-444.7738727200452)
# C_4 = df.Constant(-0.12443443849461801)
# C_5 = df.Constant(39.38867688999618)
# C_6 = df.Constant(-0.6917293233082705)
# C_7 = df.Constant(0)
# C_8 = df.Constant(0)
# C_9 = df.Constant(0)
# C_10 = df.Constant(0)
# C_11 = df.Constant(0)
# C_12 = df.Constant(2.255312046238658E-11)
# C_13 = df.Constant(407.2248457002586)
# C_14 = df.Constant(-104.89346597195336)
# C_15 = df.Constant(4.870715709115059E-7)

# lambda_1 = df.Constant(df.sqrt(5/9))
# lambda_2 = df.Constant(df.sqrt(5/6))
# lambda_3 = df.Constant(df.sqrt(3/2))

# gamma_0 = 1
# gamma = 1

# A_1 = df.Constant(0.4)

# # # x = df.SpatialCoordinate(mesh_)
# # class PressureExact(df.UserExpression):
# #     def eval(self, values, x):
# #         C_0 = df.Constant(-50.80230139855979)
# #         values[0] = ufl.atan_2(x[1], x[0])*C_0
# # p_e_i = df.interpolate(PressureExact(), space)
# # # p_e_i = df.interpolate(df.Expression("2*bessel", degree=2, bessel=PressureExact()), space)
# # p_e_i.rename("p_e", "p_e")
# # file_p = df.File(output_folder + "p_e.pvd")
# # file_p.write(p_e_i)

# d_0 = C_9 + C_2*ufl.bessel_K(0, R*lambda_2/tau) + C_8*ufl.bessel_I(0, R*lambda_2/tau)
# d = - (10*A_1 * R**2)/(27*tau) + (4*C_4*R)/(tau) - (2*C_5*tau)/R + C_14*ufl.bessel_K(1, R*lambda_2/tau) + C_15* ufl.bessel_I(1, R*lambda_2/tau)
# p = d_0 + ufl.cos(phi) * d
# # test = ufl.bessel_K(0,1)
# # p = df.Expression("d_0 + cos(phi) * d", degree=2, R=R, phi=phi, d_0=d_0, d=d)

# # p_e_i = df.project(p, space)

# # p_e_i = df.interpolate(phi, space)
# # p_e_i = df.project(test, space)
# # p_e_i = df.project(phi, space)

# x = df.SpatialCoordinate(mesh_)
# # p_e_i = df.project(ufl.atan_2(x[1], x[0]), space)

