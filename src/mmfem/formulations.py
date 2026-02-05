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
    grad, div, InnerProduct, dx, Grad,SymbolicBFI, Integrate,L2
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

def stokesp1p1_problem(
    mesh: Mesh,
    fespace: FESpace,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    viscosity: float = 1.0
) -> GridFunction:
    """
    Docstring for stokes_problem
    
    :param mesh: Description
    :type mesh: Mesh
    :param fespace: Description
    :type fespace: FESpace
    :param dirichlet_boundaries: Description
    :type dirichlet_boundaries: str
    :param velocity_bc: Description
    :type velocity_bc: CF
    :param viscosity: Description
    :type viscosity: float
    :return: Description
    :rtype: Any
    """
    # Extract test and trial functions
    (u, p, lam), (v, q, mu) = fespace.TnT()

    fes_P0 = L2(mesh, order=0)
    Π = GridFunction(fes_P0)
    #Π_p = Integrate(p, mesh) / Integrate(1, mesh)
    #Π_q = Integrate(q, mesh) / Integrate(1, mesh)

    # Bilinear form: a(u,p,λ; v,q,μ)
    a = BilinearForm(fespace)
    a += viscosity * InnerProduct(grad(u), grad(v)) * dx  # Viscous term
    a += -div(u) * q * dx                                   # Continuity (∇·u = 0)
    a += -div(v) * p * dx                                   # Pressure gradient
    a += - (p-Π) * (q-Π) * dx                             
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

def newton_step_p1p1(
    mesh: Mesh,
    fespace: FESpace,
    velocity_prev: GridFunction,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    viscosity: float = 1.0
) -> GridFunction:
    """
    Docstring for newton_step
    
    :param mesh: Description
    :type mesh: Mesh
    :param fespace: Description
    :type fespace: FESpace
    :param velocity_prev: Description
    :type velocity_prev: GridFunction
    :param dirichlet_boundaries: Description
    :type dirichlet_boundaries: str
    :param velocity_bc: Description
    :type velocity_bc: CF
    :param viscosity: Description
    :type viscosity: float
    :return: Description
    :rtype: Any
    """
    # Extract test and trial functions
    (u, p, lam), (v, q, mu) = fespace.TnT()

    #Π_p = Integrate(p, mesh) / Integrate(1, mesh)
    #Π_q = Integrate(q, mesh) / Integrate(1, mesh)
    fes_P0 = L2(mesh, order=0)
    Π = GridFunction(fes_P0)

    # Bilinear form (Jacobian)
    a = BilinearForm(fespace)

    # Stokes part
    a += viscosity * InnerProduct(Grad(u), Grad(v)) * dx
    a += -div(u) * q * dx
    a += -div(v) * p * dx
    a += - (p-Π) * (q-Π) * dx                   
    a += -lam * q * dx
    a += -mu * p * dx

    # Linearized convection term (derivative of (u·∇)u)
    # ∂/∂u [(u·∇)u] · δu = (u_old·∇)δu + (δu·∇)u_old
    #convection = (
    #    InnerProduct(Grad(u) * velocity_prev, v) * dx +
    #    InnerProduct(Grad(velocity_prev) * u, v) * dx -
    #    InnerProduct(Grad(velocity_prev) * velocity_prev, v) * dx
    #)
    convection = (
        0.5*(InnerProduct(Grad(u)*velocity_prev,v)*dx-InnerProduct(Grad(u)*v,velocity_prev)*dx) +
        0.5*(InnerProduct(Grad(velocity_prev)*u,v)*dx-InnerProduct(Grad(velocity_prev)*v,u)*dx) -
        0.5*(InnerProduct(Grad(velocity_prev)*velocity_prev,v)*dx-InnerProduct(Grad(velocity_prev)*v,velocity_prev)*dx)
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

def unsteady_stokes_step(
    mesh: Mesh,
    fespace: FESpace,
    velocity_prev: GridFunction,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    dt: float,
    viscosity: float = 1.0,
    theta: float = 1.0
) -> GridFunction:
    """
    Perform one time step for the unsteady Stokes problem.

    Solves the time-dependent Stokes equations using theta-method:

    .. math::

        \\frac{u^{n+1} - u^n}{\\Delta t} - \\nu \\Delta u^{n+\\theta} + \\nabla p^{n+1} = 0 \\\\
        \\nabla \\cdot u^{n+1} = 0

    where :math:`u^{n+\\theta} = \\theta u^{n+1} + (1-\\theta) u^n`.

    Parameters
    ----------
    mesh : Mesh
        Computational domain
    fespace : FESpace
        Product finite element space [V × Q × N]
    velocity_prev : GridFunction
        Velocity at previous time step u^n
    dirichlet_boundaries : str
        Boundary labels where Dirichlet conditions apply
    velocity_bc : CF
        Dirichlet boundary data for velocity
    dt : float
        Time step size Δt
    viscosity : float, optional
        Kinematic viscosity ν (default: 1.0)
    theta : float, optional
        Time discretization parameter (default: 1.0)
        - theta = 0.0: Forward Euler (explicit)
        - theta = 0.5: Crank-Nicolson (second order)
        - theta = 1.0: Backward Euler (implicit, unconditionally stable)

    Returns
    -------
    GridFunction
        Solution at new time step (u^{n+1}, p^{n+1}, λ^{n+1})

    Examples
    --------
    >>> # Time stepping loop
    >>> u_n = initial_velocity
    >>> for n in range(num_steps):
    ...     solution = unsteady_stokes_step(mesh, X, u_n, "walls", u_bc, dt=0.01)
    ...     u_n = solution.components[0]

    Notes
    -----
    The weak formulation finds :math:`(u^{n+1}, p^{n+1}, \\lambda^{n+1})`
    such that

    .. math::

        \\frac{1}{\\Delta t}(u^{n+1}, v) + \\theta \\nu (\\nabla u^{n+1}, \\nabla v)
        - (p^{n+1}, \\nabla \\cdot v) = \\frac{1}{\\Delta t}(u^n, v) 
        + (1-\\theta) \\nu (\\nabla u^n, \\nabla v)

    for all test functions v, plus the incompressibility constraint.
    """
    # Extract test and trial functions
    (u, p, lam), (v, q, mu) = fespace.TnT()

    # Bilinear form: implicit terms
    a = BilinearForm(fespace)
    
    # Time derivative (implicit part: u^{n+1}/dt)
    a += (1.0 / dt) * InnerProduct(u, v) * dx
    
    # Viscous term (theta-method: theta * viscosity * grad(u^{n+1}))
    a += theta * viscosity * InnerProduct(grad(u), grad(v)) * dx
    
    # Pressure and incompressibility
    a += -div(u) * q * dx
    a += -div(v) * p * dx
    a += -lam * q * dx
    a += -mu * p * dx
    
    a.Assemble()

    # Linear form: explicit terms (right-hand side)
    f = LinearForm(fespace)
    
    # Time derivative (explicit part: u^n/dt)
    f += (1.0 / dt) * InnerProduct(velocity_prev, v) * dx
    
    # Viscous term (explicit part: (1-theta) * viscosity * grad(u^n))
    if theta < 1.0:
        f += (1.0 - theta) * viscosity * InnerProduct(grad(velocity_prev), grad(v)) * dx
    
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


def unsteady_navier_stokes_semiimplicit_step(
    mesh: Mesh,
    fespace: FESpace,
    velocity_prev: GridFunction,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    dt: float,
    viscosity: float = 1.0,
    convection_form: ConvectionType = "standard"
) -> GridFunction:
    """
    Perform one semi-implicit time step for Navier-Stokes.

    Uses backward Euler for diffusion and explicit treatment of convection:

    .. math::

        \\frac{u^{n+1} - u^n}{\\Delta t} - \\nu \\Delta u^{n+1} 
        + (u^n \\cdot \\nabla)u^n + \\nabla p^{n+1} = 0

    This is unconditionally stable for diffusion but has CFL restriction
    for convection.

    Parameters
    ----------
    mesh : Mesh
        Computational domain
    fespace : FESpace
        Product finite element space [V × Q × N]
    velocity_prev : GridFunction
        Velocity at previous time step u^n
    dirichlet_boundaries : str
        Boundary labels where Dirichlet conditions apply
    velocity_bc : CF
        Dirichlet boundary data for velocity
    dt : float
        Time step size Δt
    viscosity : float, optional
        Kinematic viscosity ν (default: 1.0)
    convection_form : {'standard', 'divergence', 'skew_symmetric'}, optional
        Form of the convection term (default: 'standard')

    Returns
    -------
    GridFunction
        Solution at new time step

    Notes
    -----
    Semi-implicit scheme:
    - Implicit: time derivative + viscous term → unconditionally stable
    - Explicit: convection term → CFL condition required
    - Linear system to solve at each time step (no iterations needed)
    - Efficient for moderate Reynolds numbers
    """
    # Extract test and trial functions
    (u, p, lam), (v, q, mu) = fespace.TnT()

    # Bilinear form: implicit terms
    a = BilinearForm(fespace)
    
    # Time derivative
    a += (1.0 / dt) * InnerProduct(u, v) * dx
    
    # Viscous term (implicit)
    a += viscosity * InnerProduct(grad(u), grad(v)) * dx
    
    # Pressure and incompressibility
    a += -div(u) * q * dx
    a += -div(v) * p * dx
    a += -lam * q * dx
    a += -mu * p * dx
    
    a.Assemble()

    # Linear form: explicit terms
    f = LinearForm(fespace)
    
    # Time derivative (u^n term)
    f += (1.0 / dt) * InnerProduct(velocity_prev, v) * dx
    
    # Convection term (explicit: evaluated at u^n)
    if convection_form == "standard":
        convection_rhs = -InnerProduct(Grad(velocity_prev) * velocity_prev, v) * dx
    elif convection_form == "divergence":
        convection_rhs = (
            -InnerProduct(Grad(velocity_prev) * velocity_prev, v) * dx -
            0.5 * InnerProduct(div(velocity_prev) * velocity_prev, v) * dx
        )
    elif convection_form == "skew_symmetric":
        convection_rhs = (
            -0.5 * InnerProduct(Grad(velocity_prev) * velocity_prev, v) * dx +
            0.5 * InnerProduct(Grad(velocity_prev) * v, velocity_prev) * dx
        )
    else:
        raise ValueError(
            f"Invalid convection_form '{convection_form}'. "
            f"Must be 'standard', 'divergence', or 'skew_symmetric'."
        )
    
    f += convection_rhs
    f.Assemble()

    # Solution vector
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


def unsteady_navier_stokes_implicit_step(
    mesh: Mesh,
    fespace: FESpace,
    velocity_prev: GridFunction,
    velocity_guess: GridFunction,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    dt: float,
    viscosity: float = 1.0,
    method: Literal["picard", "newton"] = "picard"
) -> GridFunction:
    """
    Perform one fully implicit time step for Navier-Stokes using linearization.

    Solves the nonlinear system at each time step:

    .. math::

        \\frac{u^{n+1} - u^n}{\\Delta t} - \\nu \\Delta u^{n+1} 
        + (u^{n+1} \\cdot \\nabla)u^{n+1} + \\nabla p^{n+1} = 0

    using either Picard or Newton linearization.

    Parameters
    ----------
    mesh : Mesh
        Computational domain
    fespace : FESpace
        Product finite element space [V × Q × N]
    velocity_prev : GridFunction
        Velocity at previous time step u^n
    velocity_guess : GridFunction
        Initial guess for u^{n+1} (for linearization)
    dirichlet_boundaries : str
        Boundary labels where Dirichlet conditions apply
    velocity_bc : CF
        Dirichlet boundary data for velocity
    dt : float
        Time step size Δt
    viscosity : float, optional
        Kinematic viscosity ν (default: 1.0)
    method : {'picard', 'newton'}, optional
        Linearization method (default: 'picard')

    Returns
    -------
    GridFunction
        Updated solution (one linearization step)

    Notes
    -----
    This function performs ONE linearization step. For full convergence,
    call iteratively until ||u^{k+1} - u^k|| < tolerance.
    
    Use `time_stepping_implicit` solver for automatic iteration management.
    """
    # Extract test and trial functions
    (u, p, lam), (v, q, mu) = fespace.TnT()

    # Bilinear form
    a = BilinearForm(fespace)
    
    # Time derivative
    a += (1.0 / dt) * InnerProduct(u, v) * dx
    
    # Viscous term
    a += viscosity * InnerProduct(Grad(u), Grad(v)) * dx
    
    # Pressure and incompressibility
    a += -div(u) * q * dx
    a += -div(v) * p * dx
    a += -lam * q * dx
    a += -mu * p * dx

    # Convection term (linearized)
    if method == "picard":
        # Picard: (u_guess · ∇)u
        convection = InnerProduct(Grad(u) * velocity_guess, v) * dx
    elif method == "newton":
        # Newton: (u_guess · ∇)u + (u · ∇)u_guess - (u_guess · ∇)u_guess
        convection = (
            InnerProduct(Grad(u) * velocity_guess, v) * dx +
            InnerProduct(Grad(velocity_guess) * u, v) * dx -
            InnerProduct(Grad(velocity_guess) * velocity_guess, v) * dx
        )
    else:
        raise ValueError(f"Invalid method '{method}'. Must be 'picard' or 'newton'.")
    
    a += convection
    a.Assemble()

    # Linear form: time derivative term from u^n
    f = LinearForm(fespace)
    f += (1.0 / dt) * InnerProduct(velocity_prev, v) * dx
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