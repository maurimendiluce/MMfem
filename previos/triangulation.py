from netgen import gui
import netgen.geom2d as geom2d
from ngsolve import Mesh

class RectangleMesh:
    """
    Malla de un rectángulo [a,b] x [c,d].
    """

    def __init__(self, a, b, c, d, h):
        """
        Parameters
        ----------
        a, b, c, d : float
            Vértices del rectángulo.
        h : float
            Tamaño característico de la malla.
        """
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.h = h

        self.mesh = None

    def generate(self):
        """
        Genera la malla del rectángulo.

        Returns
        -------
        mesh : ngsolve.Mesh
        """
        geo = geom2d.SplineGeometry()

        p1 = geo.AppendPoint(self.a, self.c)
        p2 = geo.AppendPoint(self.b, self.c)
        p3 = geo.AppendPoint(self.b, self.d)
        p4 = geo.AppendPoint(self.a, self.d)

        geo.Append(["line", p1, p2], bc="floor")
        geo.Append(["line", p2, p3], bc="right")
        geo.Append(["line", p3, p4], bc="top")
        geo.Append(["line", p4, p1], bc="left")

        self.mesh = Mesh(geo.GenerateMesh(maxh=self.h))
        return self.mesh

    def etiquetas(self):
        """
        Devuelve las etiquetas de borde del dominio.

        Returns
        -------
        labels : str
            Etiquetas de borde separadas por '|'.
        """
        if self.mesh is None:
            raise RuntimeError("Primero debe generarse la malla con generate().")

        return "|".join(self.mesh.GetBoundaries())