from ufl import *

degree = 4
family = "N1curl"
cell_type = tetrahedron
element = FiniteElement(family, cell_type, degree)
coord_element = VectorElement("Lagrange", cell_type, 1)
mesh = Mesh(coord_element)
V = FunctionSpace(mesh, element)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(curl(u), curl(v))*dx