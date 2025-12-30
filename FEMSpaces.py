from ngsolve import *
from triangulation import *

def taylor_hood(mesh,order=2):
    V = VectorH1(mesh, order=order, dirichlet=etiquetas(mesh))
    Q = H1(mesh, order=1)
    N = NumberSpace(mesh)
    return FESpace([V,Q,N])

def mini_elements(mesh):
    V = VectorH1(mesh, order=1, dirichlet=etiquetas(mesh))
    V.SetOrder(TRIG,3)
    V.Update()
    Q = H1(mesh, order=1)
    N = NumberSpace(mesh)
    return FESpace([V,Q,N])