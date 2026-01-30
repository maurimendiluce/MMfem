"""
Mesh generation and boundary management for 2D finite element problems.

This module provides utilities for creating triangular meshes for common
geometries used in finite element simulations, together with helper functions
for boundary label management compatible with NGSolve.

Main features
-------------
- Rectangular domains with customizable boundary names
- L-shaped (backward-facing step) benchmark domain
- Utilities for extracting and composing boundary condition labels
"""

from typing import Dict, List, Union, Iterable
import netgen.geom2d as geom2d
from ngsolve import Mesh
#from ngsolve.webgui import Draw

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

_DEFAULT_RECT_BOUNDARIES: Dict[str, str] = {
    "bottom": "bottom",
    "right": "right",
    "top": "top",
    "left": "left",
}

_ALLOWED_RECT_KEYS = set(_DEFAULT_RECT_BOUNDARIES.keys())


def _validate_boundary_names(boundary_names: Dict[str, str]) -> None:
    """Validate boundary name keys for rectangular domains."""
    unknown = set(boundary_names) - _ALLOWED_RECT_KEYS
    if unknown:
        raise ValueError(
            f"Unknown boundary keys {unknown}. "
            f"Allowed keys are {_ALLOWED_RECT_KEYS}."
        )


def _append_polygon(
    geo: geom2d.SplineGeometry,
    points: List[int],
    boundaries: Iterable[str],
) -> None:
    """
    Append polygon edges to a SplineGeometry.

    Parameters
    ----------
    geo : SplineGeometry
        Netgen geometry object.
    points : list of int
        Point indices defining the polygon (counterclockwise).
    boundaries : iterable of str
        Boundary labels for each edge.
    """
    for (p0, p1), bc in zip(
        zip(points, points[1:] + points[:1]), boundaries
    ):
        geo.Append(["line", p0, p1], bc=bc)


# -----------------------------------------------------------------------------
# Mesh generators
# -----------------------------------------------------------------------------

def rectangle_mesh(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    h: float,
    boundary_names: Dict[str, str] | None = None,
) -> Mesh:
    """
    Generate a triangular mesh for a rectangular domain
    [x_min, x_max] × [y_min, y_max].

    Parameters
    ----------
    x_min, x_max : float
        Left and right x-coordinates.
    y_min, y_max : float
        Bottom and top y-coordinates.
    h : float
        Maximum mesh size (element diameter).
    boundary_names : dict, optional
        Mapping of boundary positions to labels.
        Allowed keys: {"bottom", "right", "top", "left"}.

    Returns
    -------
    Mesh
        NGSolve mesh with labeled boundaries.

    Examples
    --------
    >>> mesh = rectangle_mesh(0, 1, 0, 1, h=0.1)
    >>> mesh.GetBoundaries()
    ('bottom', 'right', 'top', 'left')
    """
    if boundary_names is None:
        boundary_names = _DEFAULT_RECT_BOUNDARIES.copy()
    else:
        _validate_boundary_names(boundary_names)
        boundary_names = _DEFAULT_RECT_BOUNDARIES | boundary_names

    geo = geom2d.SplineGeometry()

    points = [
        geo.AppendPoint(x_min, y_min),
        geo.AppendPoint(x_max, y_min),
        geo.AppendPoint(x_max, y_max),
        geo.AppendPoint(x_min, y_max),
    ]

    boundaries = [
        boundary_names["bottom"],
        boundary_names["right"],
        boundary_names["top"],
        boundary_names["left"],
    ]

    _append_polygon(geo, points, boundaries)

    return Mesh(geo.GenerateMesh(maxh=h))


def l_shaped_mesh(h: float) -> Mesh:
    """
    Generate a triangular mesh for an L-shaped (backward-facing step) domain.

    The domain is:
        [0,10] × [-0.5,1] \\ [0,4] × [-0.5,0]

    This geometry is commonly used as a CFD benchmark and exhibits
    a re-entrant corner singularity.

    Parameters
    ----------
    h : float
        Maximum mesh size.

    Returns
    -------
    Mesh
        NGSolve mesh with labeled boundaries:
        - "side"      : wall boundaries
        - "flow_in"   : inlet (x = 0)
        - "flow_out"  : outlet (x = 10)
    """
    geo = geom2d.SplineGeometry()

    points = [
        geo.AppendPoint(0, 0),
        geo.AppendPoint(4, 0),
        geo.AppendPoint(4, -0.5),
        geo.AppendPoint(10, -0.5),
        geo.AppendPoint(10, 1),
        geo.AppendPoint(0, 1),
    ]

    boundaries = [
        "side",
        "side",
        "side",
        "flow_out",
        "side",
        "flow_in",
    ]

    _append_polygon(geo, points, boundaries)

    return Mesh(geo.GenerateMesh(maxh=h))


# Backwards compatibility
L_mesh = l_shaped_mesh


# -----------------------------------------------------------------------------
# Boundary utilities
# -----------------------------------------------------------------------------

def get_boundary_labels(mesh: Mesh) -> List[str]:
    """
    Return all boundary labels of a mesh.

    Parameters
    ----------
    mesh : Mesh
        NGSolve mesh object.

    Returns
    -------
    list of str
        Boundary labels.
    """
    return list(mesh.GetBoundaries())


def boundary_string(mesh: Mesh) -> str:
    """
    Return boundary labels as a pipe-separated string.

    This format is required by NGSolve for boundary condition specification.
    """
    return "|".join(get_boundary_labels(mesh))


def get_dirichlet_boundaries(
    mesh: Mesh,
    neumann_labels: Union[str, List[str]],
) -> str:
    """
    Compute Dirichlet boundary labels by excluding Neumann boundaries.

    Parameters
    ----------
    mesh : Mesh
        NGSolve mesh object.
    neumann_labels : str or list of str
        Boundary labels where Neumann conditions apply.

    Returns
    -------
    str
        Pipe-separated Dirichlet boundary labels.

    Raises
    ------
    ValueError
        If any Neumann label is not found in the mesh.
    """
    if isinstance(neumann_labels, str):
        neumann_labels = [neumann_labels]

    labels = set(mesh.GetBoundaries())
    neumann_labels = set(neumann_labels)

    missing = neumann_labels - labels
    if missing:
        raise ValueError(
            f"Neumann labels {missing} not found. "
            f"Available boundaries: {labels}"
        )

    return "|".join(labels - neumann_labels)


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

__all__ = [
    "rectangle_mesh",
    "l_shaped_mesh",
    "L_mesh",
    "get_boundary_labels",
    "boundary_string",
    "get_dirichlet_boundaries",
]
