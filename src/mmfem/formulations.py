"""
Variational formulations for Stokes and Navier-Stokes problems.

This module contains functions to assemble and solve variational formulations of incompressible flow problems using finite element methods.

Functions:

    stokes_problem: Solve the linear Stokes problem
    picard_step: Perform one Picard linearization step
    newton_step: Perform one Newton linearization step

"""

from typing import Callable, Literal
from ngsolve import (
    Mesh, FESpace, GridFunction, BilinearForm, LinearForm, CF,
    grad, div, InnerProduct, dx, Grad
)


ConvectionType = Literal["standard", "divergence", "skew_symmetric"]


def stokes_problem(
    mesh: Mesh,
    fespace: FESpace,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    viscosity: float = 1.0
) -> GridFunction:
    """
    Solve the steady Stokes problem.

    The Stokes equations model slow, viscous, incompressible flow:

    .. math::

        -\\nu \\Delta u + \\nabla p = 0 \\quad \\text{in } \\Omega \\\\
        \\nabla \\cdot u = 0 \\quad \\text{in } \\Omega \\\\
        u = u_D \\quad \\text{on } \\Gamma_D

    Parameters
    ----------
    mesh : Mesh
        Computational domain
    fespace : FESpace
        Product finite element space [V × Q × N]
    dirichlet_boundaries : str
        Boundary labels where Dirichlet conditions apply
    velocity_bc : CF
        Dirichlet boundary data for velocity
    viscosity : float, optional
        Kinematic viscosity ν (default: 1.0)

    Returns
    -------
    GridFunction
        Solution containing velocity, pressure, and Lagrange multiplier.

    Examples
    --------
    >>> from ngsolve import CF
    >>> from mmfem import rectangle_mesh, taylor_hood, stokes_problem
    >>> mesh = rectangle_mesh(0, 1, 0, 1, h=0.05)
    >>> X = taylor_hood(mesh, "left|right|floor|top")
    >>> u_bc = CF((1, 0))
    >>> solution = stokes_problem(mesh, X, "top", u_bc, viscosity=0.01)

    Notes
    -----
    The weak formulation finds :math:`(u, p, \\lambda) \\in V \\times Q \\times \\mathbb{R}`
    such that

    .. math::

        \\nu (\\nabla u, \\nabla v)
        - (p, \\nabla \\cdot v)
        - (\\nabla \\cdot u, q)
        - \\lambda q
        - \\mu p = 0

    for all test functions :math:`(v, q, \\mu)`.

    The Lagrange multiplier enforces the constraint

    .. math::

        \\int_\\Omega p \\, dx = 0
    """
    # Extract test and trial functions
    (u, p, lam), (v, q, mu) = fespace.TnT()

    # Bilinear form: a(u,p,λ; v,q,μ)
    a = BilinearForm(fespace)
    a += viscosity * InnerProduct(grad(u), grad(v)) * dx  # Viscous term
    a += -div(u) * q * dx                                   # Continuity (∇·u = 0)
    a += -div(v) * p * dx                                   # Pressure gradient
    a += -lam * q * dx                                      # Mean value constraint
    a += -mu * p * dx                                       # Mean value constraint
    a.Assemble()

    # Linear form: L(v,q,μ) = 0 (no body forces)
    f = LinearForm(fespace)
    f.Assemble()

    # Solution vector
    solution = GridFunction(fespace)

    # Apply Dirichlet boundary conditions
    solution.components[0].Set(
        velocity_bc,
        definedon=mesh.Boundaries(dirichlet_boundaries)
    )

    # Solve the linear system
    residual = f.vec.CreateVector()
    residual.data = f.vec - a.mat * solution.vec
    
    inverse = a.mat.Inverse(freedofs=fespace.FreeDofs())
    solution.vec.data += inverse * residual

    return solution


def picard_step(
    mesh: Mesh,
    fespace: FESpace,
    velocity_prev: GridFunction,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    viscosity: float = 1.0,
    convection_form: ConvectionType = "standard"
) -> GridFunction:
    """
    Perform one Picard iteration step for the Navier-Stokes equations.

    The Picard method linearizes the nonlinear convection term by using
    the velocity from the previous iteration, solving:

    .. math::

        -\\nu \\Delta u + (u_{old} \\cdot \\nabla)u + \\nabla p = 0
        \\nabla \\cdot u = 0

    Parameters
    ----------
    mesh : Mesh
        Computational domain
    fespace : FESpace
        Product finite element space [V × Q × N]
    velocity_prev : GridFunction
        Velocity from previous iteration (for linearization)
    dirichlet_boundaries : str
        Boundary labels where Dirichlet conditions apply
    velocity_bc : CF
        Dirichlet boundary data for velocity
    viscosity : float, optional
        Kinematic viscosity ν (default: 1.0)
    convection_form : {'standard', 'divergence', 'skew_symmetric'}, optional
        Form of the convection term:
        - 'standard': (u_old·∇)u
        - 'divergence': (u_old·∇)u + 0.5(∇·u_old)u
        - 'skew_symmetric': 0.5[(u_old·∇)u - (∇u_old)·u]

    Returns
    -------
    GridFunction
        Updated solution (u, p, λ)

    Examples
    --------
    >>> # Picard iteration
    >>> u_old = stokes_solution.components[0]
    >>> for i in range(max_iter):
    ...     solution = picard_step(mesh, X, u_old, "top", u_bc, viscosity=0.01)
    ...     u_new = solution.components[0]
    ...     error = compute_error(u_new, u_old)
    ...     if error < tolerance:
    ...         break
    ...     u_old = u_new

    Notes
    -----
    The Picard method (also called fixed-point iteration) is:
    - Linearly convergent
    - Robust for low Reynolds numbers
    - Requires good initial guess for high Reynolds numbers

    Different convection forms can improve stability:
    - 'divergence': Better for non-solenoidal velocity approximations
    - 'skew_symmetric': Improved energy conservation properties
    """
    # Extract test and trial functions
    (u, p, lam), (v, q, mu) = fespace.TnT()

    # Bilinear form
    a = BilinearForm(fespace)

    # Stokes part (linear)
    a += viscosity * InnerProduct(Grad(u), Grad(v)) * dx
    a += -div(u) * q * dx
    a += -div(v) * p * dx
    a += -lam * q * dx
    a += -mu * p * dx

    # Convection term (linearized with u_old)
    if convection_form == "standard":
        convection = InnerProduct(Grad(u) * velocity_prev, v) * dx
    elif convection_form == "divergence":
        convection = (
            InnerProduct(Grad(u) * velocity_prev, v) * dx +
            0.5 * InnerProduct(div(velocity_prev) * u, v) * dx
        )
    elif convection_form == "skew_symmetric":
        convection = (
            0.5 * InnerProduct(Grad(u) * velocity_prev, v) * dx -
            0.5 * InnerProduct(Grad(velocity_prev) * v, u) * dx
        )
    else:
        raise ValueError(
            f"Invalid convection_form '{convection_form}'. "
            f"Must be 'standard', 'divergence', or 'skew_symmetric'."
        )

    a += convection
    a.Assemble()

    # Linear form (no body forces)
    f = LinearForm(fespace)
    f.Assemble()

    # Solution
    solution = GridFunction(fespace)
    solution.components[0].Set(
        velocity_bc,
        definedon=mesh.Boundaries(dirichlet_boundaries)
    )

    # Solve
    residual = f.vec.CreateVector()
    residual.data = f.vec - a.mat * solution.vec
    inverse = a.mat.Inverse(freedofs=fespace.FreeDofs())
    solution.vec.data += inverse * residual

    return solution


def newton_step(
    mesh: Mesh,
    fespace: FESpace,
    velocity_prev: GridFunction,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    viscosity: float = 1.0
) -> GridFunction:
    """
    Perform one Newton iteration step for the Navier-Stokes equations.

    The Newton method uses the full Jacobian of the nonlinear system,
    solving for the update math:`\delta u`:

    .. math::

        J(u_{old})[\\delta u] = -R(u_{old})
    
    where J is the Jacobian and R is the residual.

    Parameters
    ----------
    mesh : Mesh
        Computational domain
    fespace : FESpace
        Product finite element space [V × Q × N]
    velocity_prev : GridFunction
        Velocity from previous iteration
    dirichlet_boundaries : str
        Boundary labels where Dirichlet conditions apply
    velocity_bc : CF
        Dirichlet boundary data for velocity
    viscosity : float, optional
        Kinematic viscosity ν (default: 1.0)

    Returns
    -------
    GridFunction
        Updated solution (u, p, λ)

    Examples
    --------
    >>> # Newton iteration
    >>> u_old = stokes_solution.components[0]
    >>> for i in range(max_iter):
    ...     solution = newton_step(mesh, X, u_old, "top", u_bc, viscosity=0.01)
    ...     u_new = solution.components[0]
    ...     error = compute_error(u_new, u_old)
    ...     if error < tolerance:
    ...         break
    ...     u_old = u_new

    Notes
    -----
    The Newton method:
    - Has quadratic convergence near the solution
    - Requires good initial guess (typically from Stokes or Picard)
    - Can fail to converge for high Reynolds numbers without continuation

    The linearization includes both convection terms:

    ..math::

        (u_{old} \\cdot \\nabla)u + (u\\cdot \\nabla)u_{old}
    
    which results from differentiating :math`(u\\cdot \\nabla)u`.
    """
    # Extract test and trial functions
    (u, p, lam), (v, q, mu) = fespace.TnT()

    # Bilinear form (Jacobian)
    a = BilinearForm(fespace)

    # Stokes part
    a += viscosity * InnerProduct(Grad(u), Grad(v)) * dx
    a += -div(u) * q * dx
    a += -div(v) * p * dx
    a += -lam * q * dx
    a += -mu * p * dx

    # Linearized convection term (derivative of (u·∇)u)
    # ∂/∂u [(u·∇)u] · δu = (u_old·∇)δu + (δu·∇)u_old
    convection = (
        InnerProduct(Grad(u) * velocity_prev, v) * dx +
        InnerProduct(Grad(velocity_prev) * u, v) * dx -
        InnerProduct(Grad(velocity_prev) * velocity_prev, v) * dx
    )
    
    a += convection
    a.Assemble()

    # Linear form
    f = LinearForm(fespace)
    f.Assemble()

    # Solution
    solution = GridFunction(fespace)
    solution.components[0].Set(
        velocity_bc,
        definedon=mesh.Boundaries(dirichlet_boundaries)
    )

    # Solve
    residual = f.vec.CreateVector()
    residual.data = f.vec - a.mat * solution.vec
    inverse = a.mat.Inverse(freedofs=fespace.FreeDofs())
    solution.vec.data += inverse * residual

    return solution