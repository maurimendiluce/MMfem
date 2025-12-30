import netgen.geom2d as geom2d
from ngsolve import Mesh

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