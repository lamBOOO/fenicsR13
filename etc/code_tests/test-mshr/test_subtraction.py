"""Simple mesh test: Obtains same mesh as in Westerkamp2019"""

def main():
    "Generates the mesh"

    import mshr as m
    import dolfin as d
    import matplotlib.pyplot as plt

    d.set_log_level(13) # PROGRESS

    r_1 = 0.5 # inner
    r_2 = 2.0 # outer
    res = 10 # resolution

    circle_inner = m.Circle(d.Point(0.0, 0.0), r_1)
    circle_outer = m.Circle(d.Point(0.0, 0.0), r_2)

    domain = circle_outer - circle_inner

    domain.set_subdomain(1, circle_inner)
    domain.set_subdomain(2, circle_outer)

    mesh = m.generate_mesh(domain, res)

    print("max edge length:", mesh.hmax())

    mesh_file_pvd = d.File("mesh.pvd")
    mesh_file_pvd.write(mesh)

    plt.figure()
    d.plot(mesh, title="Mesh")
    plt.show()

if __name__ == "__main__":
    main()
