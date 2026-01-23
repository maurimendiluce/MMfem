from netgen import gui
from ngsolve.webgui import Draw
from triangulation import *
from FEMSpaces import *
from problems import *
from solvers import *
from adapt import *

print("Método adaptativo: problema de Navier-Stokes con dato no regular")
print("----------------------------------------------------------------")

example = int(input("Convexo [0] | No convexo [1]:"))
if example == 0:
    mesh = rectangle_mesh(0,1,0,1,0.5)
    cant_ref = int(input("Cantidad de refinamientos: "))
    et_dirichlet = "top"
    uD = CF((1,0))
    datos = adaptative_method(cant_ref,mesh,et_dirichlet,uD,theta=0.7)
    order(datos)
    results(datos)

if example == 1:
    mesh = L_mesh(0.4)
    cant_ref = int(input("Cantidad de refinamientos: "))
    et_dirichlet = "flow"
    uD = CF((-1,0))
    datos = adaptative_method(cant_ref,mesh,et_dirichlet,uD,theta=0.5)
    order(datos)
    results(datos)

input("presione enter para terminar...")