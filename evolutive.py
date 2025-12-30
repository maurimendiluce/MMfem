from ngsolve import *

def picard_evolutive(mesh,X,u_old,etiq_dirich,uD,dt,nu=1):
    
    (u,p,lam), (v,q,mu) = X.TnT()
    NSlin = BilinearForm(X)
    NSlin += (nu*InnerProduct(Grad(u),Grad(v))+InnerProduct(Grad(u)*u_old,v)-div(u)*q-div(v)*p-lam*q-mu*p)*dx

def implicit(mesh,X,etiq_dirich,uD,nu=1):

    u0 = uD
    (u,p,lam), (v,q,mu) = X.TnT()
