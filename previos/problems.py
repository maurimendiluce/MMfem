from ngsolve import *

class Stokes:

    def __init__(self, mesh, X, nu=1.0):
        """
        Parameters
        ----------
        mesh : ngsolve.Mesh
        X : ngsolve.FESpace
            Espacio mixto velocidad-presión.
        nu : float
            Viscosidad cinemática.
        """
        self.mesh = mesh
        self.X = X
        self.nu = nu
    
    def solve_stokes(self, etiq_dirich, uD):
        """
        Resuelve el problema de Stokes estacionario.

        Returns
        -------
        gf : ngsolve.GridFunction
            Solución (u, p, lambda).
        """
        (u,p,lam), (v,q,mu) = self.X.TnT()
        stokes = BilinearForm(self.X)
        stokes +=  (self.nu*InnerProduct(grad(u), grad(v)) - div(u)*q + div(v)*p-lam*q-mu*p)*dx
        stokes.Assemble()

        F = LinearForm(self.X)

        gf = GridFunction(self.X)

        gf.components[0].Set(uD,definedon=self.mesh.Boundaries(etiq_dirich))

        res = F.vec.CreateVector()
        res.data = F.vec - stokes.mat * gf.vec
        inv = stokes.mat.Inverse(freedofs=self.X.FreeDofs())
        gf.vec.data += inv * res

        return gf