

from dolfin import * # pylint: disable=W0614

if not has_cgal():
    print ("DOLFIN must be compiled with CGAL to run this demo.")
    exit(0)

# Create empty Mesh
mesh = Mesh()

# Create list of polygonal domain vertices
domain_vertices = [Point(0.0, 0.0),
                   Point(10.0, 0.0),
                   Point(10.0, 2.0),
                   Point(8.0, 2.0),
                   Point(7.5, 1.0),
                   Point(2.5, 1.0),
                   Point(2.0, 4.0),
                   Point(0.0, 4.0),
                   Point(0.0, 0.0)]

# Generate mesh and plot
PolygonalMeshGenerator.generate(mesh, domain_vertices, 0.25);
plot(mesh, interactive=True)

# Generate 3D mesh from OFF file input (tetrahedron)
PolyhedralMeshGenerator.generate(mesh, "../tetrahedron.off", 0.05)
plot(mesh, interactive=True)

# Generate 3D mesh from OFF file input (cube)
PolyhedralMeshGenerator.generate(mesh, "../cube.off", 0.05)
plot(mesh, interactive=True)