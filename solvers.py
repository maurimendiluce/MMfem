from ngsolve import *
from problems import *

def picard_iter(mesh,X,etiq_dirich,uD,nu=1,termino_conv="convectivo"):

    #elegimos el tipo de definicion para el termino convectivo
    #opciones = {
    #    "1": "Convectivo puro",
    #    "2": "Convectivo con div",
    #    "3": "Skew-simmetric"
    #}

    #print("Opciones disponibles:")
    #for k, v in opciones.items():
    #    print(f"{k}) {v}")

    #opcion = input("Elegí la forma del término convectivo: ")

    #if opcion not in opciones:
    #    print("Opción inválida")
    #    return

    #print(f"Elegiste: {opciones[opcion]}")
    #print("Continuando ejecución...")


    #if opciones[opcion] == "Convectivo puro":
    #    elegido = "convect"
    #if opciones[opcion] == "Convectivo con div":
    #    elegido = "div"
    #if opciones[opcion] == "Skew-simmetric":
    #    elegido = "skew"
    
    tol = 1e-10
    iter = 50

    gf_stokes = Stokes(mesh,X,etiq_dirich,uD,nu=nu)
    u_old = gf_stokes.components[0]

    gf_picard = Picard(mesh,X,u_old,etiq_dirich,uD,nu=nu,lin=termino_conv)
    i=0
    while i<iter:
        error = sqrt(Integrate(InnerProduct(Grad(gf_picard.components[0])-Grad(u_old),Grad(gf_picard.components[0])-Grad(u_old)),mesh))/sqrt(Integrate(InnerProduct(Grad(gf_picard.components[0]),Grad(gf_picard.components[0])),mesh))
        print("Iterancion:",i, "| error=",error, "| Tolerancia=", tol)
        if error<tol:
            return gf_picard
        u_old = gf_picard.components[0]
        gf_picard = Picard(mesh,X,u_old,etiq_dirich,uD,nu=nu,lin=termino_conv)
        i += 1
    

    return gf_picard

def newton_iter(mesh,X,etiq_dirich,uD,nu=1):
    
    tol = 1e-10
    iter = 50

    gf_stokes = Stokes(mesh,X,etiq_dirich,uD,nu=nu)
    u_old = gf_stokes.components[0]

    gf_newton = Newton(mesh,X,u_old,etiq_dirich,uD,nu=nu)
    i=0
    while i<iter:
        error = sqrt(Integrate(InnerProduct(Grad(gf_newton.components[0])-Grad(u_old),Grad(gf_newton.components[0])-Grad(u_old)),mesh))/sqrt(Integrate(InnerProduct(Grad(gf_newton.components[0]),Grad(gf_newton.components[0])),mesh))
        print("Iterancion:",i, "| error=",error, "| Tolerancia=", tol)
        if error<tol:
            return gf_newton
        u_old = gf_newton.components[0]
        gf_newton = Newton(mesh,X,u_old,etiq_dirich,uD,nu=nu)
        i += 1
    

    return gf_newton