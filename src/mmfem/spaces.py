"""
Finite element spaces for incompressible flow problems.

This module provides functions to construct finite element spaces suitable
for solving the Navier-Stokes equations, including inf-sup stable pairs.

Functions:
    taylor_hood: Taylor-Hood elements (P2-P1 velocity-pressure pair)
    mini_elements: MINI elements (P1-bubble enriched velocity, P1 pressure)
"""

from typing import Union
from ngsolve import Mesh, FESpace, VectorH1, H1, NumberSpace, TRIG


def taylor_hood(
    mesh: Mesh,
    dirichlet_boundaries: str,
    velocity_order: int = 2
) -> FESpace:
    """
    Construct a Taylor-Hood finite element space for Navier-Stokes problems.

    The Taylor-Hood element is a popular inf-sup stable pair consisting of:
    - P_k vector-valued finite elements for velocity (default k=2)
    - P_{k-1} scalar finite elements for pressure (k-1 typically = 1)
    - NumberSpace for pressure mean value constraint

    Parameters
    ----------
    mesh : Mesh
        Computational domain mesh
    dirichlet_boundaries : str
        Pipe-separated string of boundary labels where Dirichlet conditions
        are applied (e.g., "left|right|top")
    velocity_order : int, optional
        Polynomial order for velocity space (default: 2).
        Pressure order will be velocity_order - 1

    Returns
    -------
    FESpace
        Product finite element space [V × Q × N] where:
        - V: Vector H1 space for velocity (order=velocity_order)
        - Q: Scalar H1 space for pressure (order=velocity_order-1)
        - N: NumberSpace for mean value constraint

    Examples
    --------
    >>> from mmfem.mesh import rectangle_mesh
    >>> from mmfem.spaces import taylor_hood
    >>> mesh = rectangle_mesh(0, 1, 0, 1, h=0.1)
    >>> X = taylor_hood(mesh, "left|right|top|floor")
    >>> print(f"Total DOFs: {X.ndof}")

    >>> # Higher order Taylor-Hood (P3-P2)
    >>> X_high = taylor_hood(mesh, "left|right|top|floor", velocity_order=3)

    Notes
    -----
    The Taylor-Hood P2-P1 element is inf-sup stable and widely used for
    incompressible flow simulations. The NumberSpace constraint ensures
    uniqueness of the pressure solution up to a constant.

    References
    ----------
    Taylor, C., & Hood, P. (1973). *A numerical solution of the Navier-Stokes equations using the finite element technique.*
    """
    # Velocity space (vector-valued H1)
    V = VectorH1(mesh, order=velocity_order, dirichlet=dirichlet_boundaries)
    
    # Pressure space (scalar H1, one order lower)
    pressure_order = max(1, velocity_order - 1)
    Q = H1(mesh, order=pressure_order)
    
    # NumberSpace for pressure mean value constraint
    N = NumberSpace(mesh)
    
    return FESpace([V, Q, N])


def mini_elements(
    mesh: Mesh,
    dirichlet_boundaries: str
) -> FESpace:
    """
    Construct MINI finite element space for Navier-Stokes problems.

    The MINI element consists of:
    - P1 vector-valued finite elements enriched with cubic bubbles for velocity
    - P1 scalar finite elements for pressure
    - NumberSpace for pressure mean value constraint

    The bubble enrichment is achieved by setting order=3 on triangular elements,
    which adds internal degrees of freedom that vanish on element boundaries.

    Parameters
    ----------
    mesh : Mesh
        Computational domain mesh
    dirichlet_boundaries : str
        Pipe-separated string of boundary labels where Dirichlet conditions
        are applied

    Returns
    -------
    FESpace
        Product finite element space [V × Q × N] where:
        - V: Vector H1 space with P1 + bubble enrichment
        - Q: Scalar H1 space for pressure (order=1)
        - N: NumberSpace for mean value constraint

    Examples
    --------
    >>> from mmfem.mesh import rectangle_mesh
    >>> from mmfem.spaces import mini_elements
    >>> mesh = rectangle_mesh(0, 1, 0, 1, h=0.1)
    >>> X = mini_elements(mesh, "left|right|top|floor")
    >>> print(f"Total DOFs: {X.ndof}")

    Notes
    -----
    MINI elements provide a simple, low-order inf-sup stable discretization.
    The bubble functions improve the approximation quality while maintaining
    local support. They are particularly useful for:
    - Problems requiring local mass conservation
    - Situations where Taylor-Hood elements are too expensive

    The bubble enrichment increases the number of degrees of freedom per
    element but maintains sparsity of the system matrix.

    References
    ----------
    Arnold, D., Brezzi, F., & Fortin, M. (1984). *A stable finite element for the Stokes equations.*
    """
    # Velocity space (initially P1)
    V = VectorH1(mesh, order=1, dirichlet=dirichlet_boundaries)
    
    # Enrich with cubic bubbles on triangular elements
    V.SetOrder(TRIG, 3)
    V.Update()
    
    # Pressure space (P1)
    Q = H1(mesh, order=1)
    
    # NumberSpace for pressure mean value constraint
    N = NumberSpace(mesh)
    
    return FESpace([V, Q, N])


def stabilization_p1p1(mesh: Mesh,
    dirichlet_boundaries: str
) -> FESpace:
    """
    Docstring for stabilization_p1p1
    
    :param mesh: Description
    :type mesh: Mesh
    :param dirichlet_boundaries: Description
    :type dirichlet_boundaries: str
    :return: Description
    :rtype: Any
    """
    # Velocity space (vector-valued H1)
    V = VectorH1(mesh, order=1, dirichlet=dirichlet_boundaries)
    
    # Pressure space (scalar H1, one order lower)
    Q = H1(mesh, order=1)
    
    # NumberSpace for pressure mean value constraint
    N = NumberSpace(mesh)
    
    return FESpace([V, Q, N])

# Optional: Function to get space information
def get_space_info(fespace: FESpace) -> dict:
    """
    Extract information about a finite element space.

    Parameters
    ----------
    fespace : FESpace
        Finite element space (typically a product space)

    Returns
    -------
    dict
        Dictionary containing:
        - 'total_dofs': Total degrees of freedom
        - 'velocity_dofs': Velocity DOFs (if applicable)
        - 'pressure_dofs': Pressure DOFs (if applicable)
        - 'constraint_dofs': Constraint DOFs (if applicable)

    Examples
    --------
    >>> X = taylor_hood(mesh, "left|right|top|floor")
    >>> info = get_space_info(X)
    >>> print(info)
    """
    info = {'total_dofs': fespace.ndof}
    
    # Try to extract component information
    try:
        components = fespace.components
        if len(components) >= 3:
            info['velocity_dofs'] = components[0].ndof
            info['pressure_dofs'] = components[1].ndof
            info['constraint_dofs'] = components[2].ndof
    except AttributeError:
        pass
    
    return info