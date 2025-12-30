from netgen import gui
import netgen.geom2d as geom2d
from ngsolve import *
from ngsolve.webgui import Draw

#geometrias
def rectangle_mesh(a,b,c,d,h):

    """
    Genera una malla para el rectangulo [a,b]x[c,d]

    Parameters
    ----------
    a,b,c,d = vertices
    h = tamaño de la malla

    Returns
    -------
    mesh: objeto de netgen
    Examples
    --------
    >>> mesh = rectangle_mesh(0, 1, 0, 1, 0.1)
    """
    geo = geom2d.SplineGeometry()
    p1 = geo.AppendPoint (a,c)
    p2 = geo.AppendPoint (b,c)
    p3 = geo.AppendPoint (b,d)
    p4 = geo.AppendPoint (a,d)

    geo.Append (["line", p1, p2],bc = "floor")
    geo.Append (["line", p2, p3],bc = "rigth")
    geo.Append (["line", p3, p4],bc = "top")
    geo.Append (["line", p4, p1],bc = "left")

    mesh = Mesh(geo.GenerateMesh(maxh=h))

    return mesh

def etiquetas(mesh):
    """
    Devuelve las etiquetas de borde del dominio.

    Parameters
    ----------
    mesh : netgen mesh
        Dominio discretizado.

    Returns
    -------
    labels : str
        Etiquetas de borde.
    
    Examples
    --------
    >>> mesh = rectangle_mesh(0, 1, 0, 1, 0.1)
    >>> etiquetas(mesh)
    'left|right|top|bottom'
    """
    bc = mesh.GetBoundaries()
    s = "|".join(bc)
    return s 

#armado de espacios
def espacios(mesh):
    #elegimos el tipo de espacios para hacer Elementos Finitos
    opciones = {
        "1": "Taylor-Hood P2P1",
        "2": "Mini-elements",
    }

    print("Opciones disponibles:")
    for k, v in opciones.items():
        print(f"{k}) {v}")

    opcion = input("Elegí el método a usar: ")

    if opcion not in opciones:
        print("Opción inválida")
        return

    print(f"Elegiste: {opciones[opcion]}")
    if opciones[opcion] == "Taylor-Hood P2P1":
        V = VectorH1(mesh, order=2, dirichlet=etiquetas(mesh))
        Q = H1(mesh, order=1)
        N = NumberSpace(mesh)
        X = FESpace([V,Q,N])
    if opciones[opcion] == "Mini-elements":
        V = VectorH1(mesh, order=1, dirichlet=etiquetas(mesh))
        V.SetOrder(TRIG,3)
        V.Update()
        Q = H1(mesh, order=1)
        N = NumberSpace(mesh)
        X = FESpace([V,Q,N])

    return X

#problemas
def Stokes(mesh,X,etiq_dirich,uD,nu=1):

    (u,p,lam), (v,q,mu) = X.TnT()
    stokes = BilinearForm(X)
    stokes +=  (nu*InnerProduct(grad(u), grad(v)) - div(u)*q + div(v)*p-lam*q-mu*p)*dx
    stokes.Assemble()

    F = LinearForm(X)
    

    gf = GridFunction(X)

    #u0 = CF((1,0))
    gf.components[0].Set(uD,definedon=mesh.Boundaries(etiq_dirich))

    res = F.vec.CreateVector()
    res.data = F.vec - stokes.mat * gf.vec
    inv = stokes.mat.Inverse(freedofs=X.FreeDofs())
    gf.vec.data += inv * res

    
    return gf

def Picard(mesh,X,u_old,etiq_dirich,uD,nu=1,lin="convect"):
    
    (u,p,lam), (v,q,mu) = X.TnT()
    NSlin = BilinearForm(X)
    if lin == "convect":
        conv = InnerProduct(Grad(u)*u_old,v)*dx
    elif lin == "div":
        conv = InnerProduct(Grad(u)*u_old,v)*dx+ 0.5*InnerProduct(div(u)*u_old,v)*dx
    elif lin == "skew":
        conv = 0.5*InnerProduct(Grad(u)*u_old,v)*dx-0.5*InnerProduct(Grad(u_old)*v,u)*dx

    stokes = (nu*InnerProduct(Grad(u),Grad(v))-div(u)*q-div(v)*p-lam*q-mu*p)*dx
    NSlin += stokes + conv

    NSlin.Assemble()

    F = LinearForm(X)

    gf = GridFunction(X)

    gf.components[0].Set(uD,definedon=mesh.Boundaries(etiq_dirich))

    res = F.vec.CreateVector()
    res.data = F.vec - NSlin.mat * gf.vec
    inv = NSlin.mat.Inverse(freedofs=X.FreeDofs())
    gf.vec.data += inv * res

    return gf

def Newton(mesh,X,u_old,etiq_dirich,uD,nu=1):
    
    (u,p,lam), (v,q,mu) = X.TnT()
    NSlin = BilinearForm(X)
    conv = InnerProduct(Grad(u)*u_old,v)*dx+InnerProduct(Grad(u_old)*u,v)*dx-InnerProduct(Grad(u_old)*u_old,v)*dx

    stokes = (nu*InnerProduct(Grad(u),Grad(v))-div(u)*q-div(v)*p-lam*q-mu*p)*dx
    NSlin += stokes + conv

    NSlin.Assemble()

    F = LinearForm(X)

    gf = GridFunction(X)

    gf.components[0].Set(uD,definedon=mesh.Boundaries(etiq_dirich))

    res = F.vec.CreateVector()
    res.data = F.vec - NSlin.mat * gf.vec
    inv = NSlin.mat.Inverse(freedofs=X.FreeDofs())
    gf.vec.data += inv * res

    return gf


#solvers
def NavierStokes_Picard_iter(mesh,X,etiq_dirich,uD,nu=1):

    #elegimos el tipo de definicion para el termino convectivo
    opciones = {
        "1": "Convectivo puro",
        "2": "Convectivo con div",
        "3": "Skew-simmetric"
    }

    print("Opciones disponibles:")
    for k, v in opciones.items():
        print(f"{k}) {v}")

    opcion = input("Elegí la forma del término convectivo: ")

    if opcion not in opciones:
        print("Opción inválida")
        return

    print(f"Elegiste: {opciones[opcion]}")
    print("Continuando ejecución...")


    if opciones[opcion] == "Convectivo puro":
        elegido = "convect"
    if opciones[opcion] == "Convectivo con div":
        elegido = "div"
    if opciones[opcion] == "Skew-simmetric":
        elegido = "skew"
    
    tol = 1e-10
    iter = 50

    gf_stokes = Stokes(mesh,X,etiq_dirich,uD,nu=nu)
    u_old = gf_stokes.components[0]

    gf_picard = Picard(mesh,X,u_old,etiq_dirich,uD,nu=nu,lin=elegido)
    i=0
    while i<iter:
        error = sqrt(Integrate(InnerProduct(Grad(gf_picard.components[0])-Grad(u_old),Grad(gf_picard.components[0])-Grad(u_old)),mesh))/sqrt(Integrate(InnerProduct(Grad(gf_picard.components[0]),Grad(gf_picard.components[0])),mesh))
        print("Iterancion:",i, "| error=",error, "| Tolerancia=", tol)
        if error<tol:
            return gf_picard
        u_old = gf_picard.components[0]
        gf_picard = Picard(mesh,X,u_old,etiq_dirich,uD,nu=nu,lin=elegido)
        i += 1
    

    return gf_picard

def NavierStokes_Newton_iter(mesh,X,etiq_dirich,uD,nu=1):
    
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
        gf_picard = Newton(mesh,X,u_old,etiq_dirich,uD,nu=nu)
        i += 1
    

    return gf_picard


__all__ = ["rectangle_mesh","etiquetas"]