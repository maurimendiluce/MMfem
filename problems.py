from ngsolve import *
from triangulation import *
from FEMSpaces import *

def Stokes(mesh,X,etiq_dirich,uD,nu=1):
    """
    Docstring for Stokes
    
    :param mesh: Description
    :param X: Description
    :param etiq_dirich: Description
    :param uD: Description
    :param nu: Description
    """
    (u,p,lam), (v,q,mu) = X.TnT()
    stokes = BilinearForm(X)
    stokes +=  (nu*InnerProduct(grad(u), grad(v)) - div(u)*q + div(v)*p-lam*q-mu*p)*dx
    stokes.Assemble()

    F = LinearForm(X)
    

    gf = GridFunction(X)

    gf.components[0].Set(uD,definedon=mesh.Boundaries(etiq_dirich))

    res = F.vec.CreateVector()
    res.data = F.vec - stokes.mat * gf.vec
    inv = stokes.mat.Inverse(freedofs=X.FreeDofs())
    gf.vec.data += inv * res

    
    return gf

def Picard(mesh,X,u_old,etiq_dirich,uD,nu=1,lin="convectivo"):
    """
    Docstring for Picard
    
    :param mesh: Description
    :param X: Description
    :param u_old: Description
    :param etiq_dirich: Description
    :param uD: Description
    :param nu: Description
    :param lin: Description
    """
    (u,p,lam), (v,q,mu) = X.TnT()
    NSlin = BilinearForm(X)
    if lin == "convectivo":
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
    """
    Docstring for Newton
    
    :param mesh: Description
    :param X: Description
    :param u_old: Description
    :param etiq_dirich: Description
    :param uD: Description
    :param nu: Description
    """
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