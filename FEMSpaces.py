from ngsolve import *
from triangulation import *

def taylor_hood(mesh,etiquetas,order=2):
    """
    Docstring for taylor_hood
    
    :param mesh: Mesh 
    :param order: Polinomial order
    """
    V = VectorH1(mesh, order=order, dirichlet=etiquetas)
    Q = H1(mesh, order=1)
    N = NumberSpace(mesh)
    return FESpace([V,Q,N])

def mini_elements(mesh,etiquetas):
    """
    Docstring for mini_elements
    
    :param mesh: Mesh
    """
    V = VectorH1(mesh, order=1, dirichlet=etiquetas)
    V.SetOrder(TRIG,3)
    V.Update()
    Q = H1(mesh, order=1)
    N = NumberSpace(mesh)
    return FESpace([V,Q,N])