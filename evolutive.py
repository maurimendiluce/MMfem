from ngsolve import *
from netgen.geom2d import unit_square
from ngsolve.webgui import Draw


def picard_evolutive(dt=0.01,nu=0.1):
    
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.05))
    V = VectorH1(mesh, order=2, dirichlet="left|bottom|top")
    Q = H1(mesh, order=1)

    X = V * Q
    (u,p),(v,q) = X.TnT()
    t = Parameter(0.0)
    uin = CoefficientFunction((sin(2*pi*t), 0))

    gfu = GridFunction(X)        # (u^{n+1}, p^{n+1})
    gfu_old = GridFunction(X)    # (u^n, p^n)
    gfu_old.components[0].Set((0,0))


    u_old, p_old = gfu_old.components
    a = BilinearForm(X, symmetric=False)

    # masa
    a += (1/dt) * InnerProduct(u, v) * dx

    # viscosidad
    a += nu * InnerProduct(grad(u), grad(v)) * dx

    # presión
    a += -div(v)*p * dx - div(u)*q * dx

    # convectivo semi–implícito
    a += InnerProduct(grad(u) * u_old, v) * dx

    f = CoefficientFunction((0,0))

    L = LinearForm(X)
    L += (1/dt) * InnerProduct(u_old, v) * dx
    L += InnerProduct(f, v) * dx

    a.Assemble()

    inv = a.mat.Inverse(X.FreeDofs())

    Draw(gfu.components[0], mesh, "velocity")

    tcur = 0.0
    T = 2.0

    while tcur < T:
        print(f"t = {tcur:.3f}")

        # actualizar tiempo
        t.Set(tcur)

        # imponer Dirichlet fuerte
        gfu.components[0].Set(uin, definedon=mesh.Boundaries("left"))

        # ensamblar RHS (u_old cambió)
        L.Assemble()

        # resolver
        gfu.vec.data = inv * L.vec

        # actualizar estado previo
        gfu_old.vec.data = gfu.vec

        Redraw()
        tcur += dt


def implicit(mesh,X,etiq_dirich,uD,nu=1):

    u0 = uD
    (u,p,lam), (v,q,mu) = X.TnT()

