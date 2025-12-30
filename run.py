from ngsolve.webgui import Draw
from netgen import gui
from triangulation import *
from FEMSpaces import *
from problems import *
from solvers import *

print("Resolucion por elementos finitos del problema de Navier-Stokes")
print("--------------------------------------------------------------")
print("\n")

mesh = rectangle_mesh(0,1,0,1,0.05)
et_dirichlet = "top"
uD = CF((1,0))
X = taylor_hood(mesh)
gf_NS = picard_iter(mesh,X,et_dirichlet,uD,nu=0.01)
uh = gf_NS.components[0]
Draw(uh)

input("presione enter para terminar...")