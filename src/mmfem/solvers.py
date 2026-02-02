"""
Iterative solvers for nonlinear Navier-Stokes equations.

This module provides high-level iterative solvers that combine the variational
formulations with convergence criteria and iteration control.

Functions:
    picard_iteration: Fixed-point iteration for Navier-Stokes
    newton_iteration: Newton-Raphson iteration for Navier-Stokes
    compute_velocity_error: Compute relative H1 error between velocities
"""

from typing import Tuple, Optional, Callable, Literal
import numpy as np
from ngsolve import (
    Mesh, FESpace, GridFunction, CF, Integrate, sqrt, InnerProduct, Grad, VTKOutput
)

from .formulations import stokes_problem, picard_step, newton_step, ConvectionType


def compute_velocity_error(
    velocity_new: GridFunction,
    velocity_old: GridFunction,
    mesh: Mesh
) -> float:
    """
    Compute relative H1 error between two velocity fields.

    The relative error is defined as:

    .. math::

        error = || \\nabla (u_{new} - u_{old}) ||_{L^2} / || \\nabla u_{new} ||_{L^2}

    Parameters
    ----------

    velocity_new : GridFunction
        New velocity field
    velocity_old : GridFunction
        Previous velocity field
    mesh : Mesh
        Computational domain

    Returns
    -------
    
    float
        Relative H1 error

    Examples
    --------
    >>> error = compute_velocity_error(u_new, u_old, mesh)
    >>> print(f"Relative error: {error:.2e}")
    """
    # Gradient difference
    grad_diff = Grad(velocity_new) - Grad(velocity_old)
    
    # L2 norm of gradient difference
    numerator = sqrt(Integrate(
        InnerProduct(grad_diff, grad_diff),
        mesh
    ))
    
    # L2 norm of new gradient
    denominator = sqrt(Integrate(
        InnerProduct(Grad(velocity_new), Grad(velocity_new)),
        mesh
    ))
    
    # Avoid division by zero
    if denominator < 1e-15:
        return float('inf')
    
    return numerator / denominator


def picard_iteration(
    mesh: Mesh,
    fespace: FESpace,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    viscosity: float = 1.0,
    convection_form: ConvectionType = "standard",
    tolerance: float = 1e-10,
    max_iterations: int = 100,
    verbose: bool = True,
    callback: Optional[Callable] = None
) -> Tuple[GridFunction, dict]:
    """
    Solve Navier-Stokes equations using Picard (fixed-point) iteration.

    The method iteratively solves linearized problems until convergence:
        1. Start with Stokes solution (u^0)
        2. For k = 0, 1, 2, ...:
           Solve linearized problem with u^k to get u^{k+1}
           Check convergence: ||u^{k+1} - u^k|| < tolerance
           If converged, return u^{k+1}

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
    convection_form : {'standard', 'divergence', 'skew_symmetric'}, optional
        Form of the convection term (default: 'standard')
    tolerance : float, optional
        Convergence tolerance for relative H1 error (default: 1e-10)
    max_iterations : int, optional
        Maximum number of iterations (default: 100)
    verbose : bool, optional
        Print iteration information (default: True)
    callback : callable, optional
        Function called after each iteration with signature:
        callback(iteration: int, solution: GridFunction, error: float)

    Returns
    -------
    solution : GridFunction
        Converged solution containing (velocity, pressure, multiplier)
    info : dict
        Dictionary with convergence information:
        - 'converged': bool, whether the method converged
        - 'iterations': int, number of iterations performed
        - 'errors': list of float, error at each iteration
        - 'final_error': float, final error value

    Raises
    ------
    RuntimeError
        If maximum iterations reached without convergence

    Examples
    --------
    >>> from mmfem import rectangle_mesh, taylor_hood, picard_iteration
    >>> from ngsolve import CF
    >>> 
    >>> mesh = rectangle_mesh(0, 1, 0, 1, h=0.05)
    >>> X = taylor_hood(mesh, "left|right|bottom|top")
    >>> u_bc = CF((1, 0))
    >>> 
    >>> solution, info = picard_iteration(
    ...     mesh, X, "top", u_bc,
    ...     viscosity=0.01,
    ...     tolerance=1e-8
    ... )
    >>> 
    >>> print(f"Converged in {info['iterations']} iterations")
    >>> u, p, _ = solution.components

    Notes
    -----
    Picard iteration is robust but only linearly convergent. It works well for:
    - Low to moderate Reynolds numbers
    - Problems with good initial guesses

    For high Reynolds numbers, consider:
    - Using continuation methods (gradually decrease viscosity)
    - Switching to Newton iteration after Picard converges
    """
    if verbose:
        print("=" * 70)
        print("Picard Iteration for Navier-Stokes Equations")
        print("=" * 70)
        print(f"Viscosity (ν): {viscosity}")
        print(f"Convection form: {convection_form}")
        print(f"Tolerance: {tolerance:.2e}")
        print(f"Max iterations: {max_iterations}")
        print("-" * 70)

    # Step 1: Solve Stokes problem for initial guess
    if verbose:
        print("Solving initial Stokes problem...")
    
    stokes_solution = stokes_problem(
        mesh, fespace, dirichlet_boundaries, velocity_bc, viscosity
    )
    velocity_old = stokes_solution.components[0]

    # Initialize tracking
    errors = []
    iteration = 0

    # Step 2: Picard iteration
    while iteration < max_iterations:
        # Solve linearized problem
        solution = picard_step(
            mesh, fespace, velocity_old,
            dirichlet_boundaries, velocity_bc,
            viscosity, convection_form
        )
        
        velocity_new = solution.components[0]

        # Compute error
        error = compute_velocity_error(velocity_new, velocity_old, mesh)
        errors.append(error)

        if verbose:
            print(f"Iteration {iteration:3d} | Error: {error:.6e} | Tolerance: {tolerance:.2e}")

        # Call callback if provided
        if callback is not None:
            callback(iteration, solution, error)

        # Check convergence
        if error < tolerance:
            info = {
                'converged': True,
                'iterations': iteration,
                'errors': errors,
                'final_error': error
            }
            if verbose:
                print("-" * 70)
                print(f"✓ Converged in {iteration} iterations")
                print("=" * 70)
            return solution, info

        # Update for next iteration
        velocity_old = velocity_new
        iteration += 1

    # Max iterations reached
    info = {
        'converged': False,
        'iterations': max_iterations,
        'errors': errors,
        'final_error': errors[-1] if errors else float('inf')
    }
    
    if verbose:
        print("-" * 70)
        print(f"✗ Did not converge in {max_iterations} iterations")
        print(f"Final error: {errors[-1]:.6e}")
        print("=" * 70)
    
    raise RuntimeError(
        f"Picard iteration did not converge in {max_iterations} iterations. "
        f"Final error: {errors[-1]:.6e}"
    )


def newton_iteration(
    mesh: Mesh,
    fespace: FESpace,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    viscosity: float = 1.0,
    tolerance: float = 1e-10,
    max_iterations: int = 50,
    verbose: bool = True,
    use_stokes_initial: bool = True,
    callback: Optional[Callable] = None
) -> Tuple[GridFunction, dict]:
    """
    Solve Navier-Stokes equations using Newton-Raphson iteration.

    The Newton method has quadratic convergence but requires a good initial guess.
    It solves the nonlinear system by iteratively solving linearized problems
    using the full Jacobian.

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
    tolerance : float, optional
        Convergence tolerance for relative H1 error (default: 1e-10)
    max_iterations : int, optional
        Maximum number of iterations (default: 50)
    verbose : bool, optional
        Print iteration information (default: True)
    use_stokes_initial : bool, optional
        Use Stokes solution as initial guess (default: True).
        If False, uses zero initial guess.
    callback : callable, optional
        Function called after each iteration

    Returns
    -------
    solution : GridFunction
        Converged solution containing (velocity, pressure, multiplier)
    info : dict
        Dictionary with convergence information

    Examples
    --------
    >>> solution, info = newton_iteration(
    ...     mesh, X, "top", u_bc,
    ...     viscosity=0.01,
    ...     tolerance=1e-10
    ... )

    Notes
    -----
    Newton iteration:
    - Has quadratic convergence near the solution
    - May diverge without good initial guess
    - Typically requires fewer iterations than Picard

    Best practices:
    1. Use Stokes solution as initial guess (default)
    2. For high Reynolds numbers, first converge with Picard
    3. Consider continuation methods for difficult problems
    """
    if verbose:
        print("=" * 70)
        print("Newton Iteration for Navier-Stokes Equations")
        print("=" * 70)
        print(f"Viscosity (ν): {viscosity}")
        print(f"Tolerance: {tolerance:.2e}")
        print(f"Max iterations: {max_iterations}")
        print("-" * 70)

    # Initial guess
    if use_stokes_initial:
        if verbose:
            print("Solving initial Stokes problem...")
        initial_solution = stokes_problem(
            mesh, fespace, dirichlet_boundaries, velocity_bc, viscosity
        )
    else:
        initial_solution = GridFunction(fespace)
        initial_solution.components[0].Set(
            velocity_bc,
            definedon=mesh.Boundaries(dirichlet_boundaries)
        )
    
    velocity_old = initial_solution.components[0]

    # Newton iteration
    errors = []
    iteration = 0

    while iteration < max_iterations:
        # Newton step
        solution = newton_step(
            mesh, fespace, velocity_old,
            dirichlet_boundaries, velocity_bc,
            viscosity
        )
        
        velocity_new = solution.components[0]

        # Compute error
        error = compute_velocity_error(velocity_new, velocity_old, mesh)
        errors.append(error)

        if verbose:
            print(f"Iteration {iteration:3d} | Error: {error:.6e} | Tolerance: {tolerance:.2e}")

        # Callback
        if callback is not None:
            callback(iteration, solution, error)

        # Check convergence
        if error < tolerance:
            info = {
                'converged': True,
                'iterations': iteration,
                'errors': errors,
                'final_error': error
            }
            if verbose:
                print("-" * 70)
                print(f"✓ Converged in {iteration} iterations")
                print("=" * 70)
            return solution, info

        # Update
        velocity_old = velocity_new
        iteration += 1

    # Max iterations reached
    info = {
        'converged': False,
        'iterations': max_iterations,
        'errors': errors,
        'final_error': errors[-1] if errors else float('inf')
    }
    
    if verbose:
        print("-" * 70)
        print(f"✗ Did not converge in {max_iterations} iterations")
        print(f"Final error: {errors[-1]:.6e}")
        print("=" * 70)
    
    raise RuntimeError(
        f"Newton iteration did not converge in {max_iterations} iterations. "
        f"Final error: {errors[-1]:.6e}"
    )

def time_stepping_semiimplicit(
    mesh: Mesh,
    fespace: FESpace,
    initial_velocity: GridFunction,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    dt: float,
    T_final: float,
    viscosity: float = 1.0,
    convection_form: ConvectionType = "standard",
    save_frequency: int = 1,
    verbose: bool = True,
    callback: Optional[Callable] = None
) -> Tuple[list, list]:
    """
    Solve time-dependent Navier-Stokes using semi-implicit time stepping.

    The semi-implicit scheme treats diffusion implicitly (unconditionally stable)
    and convection explicitly (requires CFL condition).

    Parameters
    ----------
    mesh : Mesh
        Computational domain
    fespace : FESpace
        Product finite element space [V × Q × N]
    initial_velocity : GridFunction
        Initial velocity field u^0
    dirichlet_boundaries : str
        Boundary labels where Dirichlet conditions apply
    velocity_bc : CF
        Dirichlet boundary data for velocity
    dt : float
        Time step size Δt
    T_final : float
        Final simulation time
    viscosity : float, optional
        Kinematic viscosity ν (default: 1.0)
    convection_form : {'standard', 'divergence', 'skew_symmetric'}, optional
        Form of the convection term (default: 'standard')
    save_frequency : int, optional
        Save solution every N time steps (default: 1)
    verbose : bool, optional
        Print progress information (default: True)
    callback : callable, optional
        Function called after each saved time step with signature:
        callback(time: float, timestep: int, solution: GridFunction)

    Returns
    -------
    solutions : list of GridFunction
        List of solutions at saved time steps
    times : list of float
        List of corresponding time values

    Examples
    --------
    >>> from mmfem import rectangle_mesh, taylor_hood, stokes_problem
    >>> from ngsolve import CF
    >>> 
    >>> # Setup
    >>> mesh = rectangle_mesh(0, 1, 0, 1, h=0.05)
    >>> X = taylor_hood(mesh, "left|right|bottom|top")
    >>> 
    >>> # Initial condition (Stokes solution)
    >>> u_bc = CF((1, 0))
    >>> stokes_sol = stokes_problem(mesh, X, "top", u_bc, viscosity=0.01)
    >>> u0 = stokes_sol.components[0]
    >>> 
    >>> # Time stepping
    >>> solutions, times = time_stepping_semiimplicit(
    ...     mesh, X, u0, "top", u_bc,
    ...     dt=0.01, T_final=1.0,
    ...     viscosity=0.01
    ... )

    Notes
    -----
    CFL condition for stability of convection term:
    
    .. math::
    
        \\Delta t \\leq \\frac{Ch}{||u||_{\\infty}}
    
    where C ≈ 1 and h is the mesh size.
    """
    from .formulations import unsteady_navier_stokes_semiimplicit_step

    if verbose:
        print("=" * 70)
        print("Semi-Implicit Time Stepping for Navier-Stokes")
        print("=" * 70)
        print(f"Time step (Δt): {dt}")
        print(f"Final time (T): {T_final}")
        print(f"Viscosity (ν): {viscosity}")
        print(f"Convection form: {convection_form}")
        print(f"Number of steps: {int(T_final / dt)}")
        print("-" * 70)

    # Initialize
    velocity_n = initial_velocity
    time = 0.0
    timestep = 0
    
    solutions = []
    times = []
    
    # Save initial condition
    solution_0 = GridFunction(fespace)
    solution_0.components[0].vec.data = velocity_n.vec
    solutions.append(solution_0)
    times.append(0.0)
    # Time stepping loop
    while time < T_final - 1e-10:
        # Adjust last time step
        if time + dt > T_final:
            dt = T_final - time
        
        # Perform one time step
        solution = unsteady_navier_stokes_semiimplicit_step(
            mesh, fespace, velocity_n,
            dirichlet_boundaries, velocity_bc,
            dt, viscosity, convection_form
        )
        # Update
        velocity_n = solution.components[0]
        time += dt
        timestep += 1
        # Verbose output
        if verbose and timestep % max(1, save_frequency) == 0:
            print(f"Step {timestep:4d} | Time: {time:6.4f} | dt: {dt:.6f}")
        
        # Save solution
        if timestep % save_frequency == 0:
            solutions.append(solution)
            times.append(time)
            
            # Call callback
            if callback is not None:
                callback(time, timestep, solution)
    
    if verbose:
        print("-" * 70)
        print(f"✓ Simulation complete: {timestep} time steps")
        print("=" * 70)
    
    return solutions, times


def time_stepping_implicit(
    mesh: Mesh,
    fespace: FESpace,
    initial_velocity: GridFunction,
    dirichlet_boundaries: str,
    velocity_bc: CF,
    dt: float,
    T_final: float,
    viscosity: float = 1.0,
    method: Literal["picard", "newton"] = "picard",
    tolerance: float = 1e-8,
    max_nonlinear_iter: int = 20,
    save_frequency: int = 1,
    verbose: bool = True,
    callback: Optional[Callable] = None
) -> Tuple[list, list]:
    """
    Solve time-dependent Navier-Stokes using fully implicit time stepping.

    The fully implicit scheme treats both diffusion and convection implicitly,
    requiring nonlinear iterations at each time step but allowing larger time steps.

    Parameters
    ----------
    mesh : Mesh
        Computational domain
    fespace : FESpace
        Product finite element space [V × Q × N]
    initial_velocity : GridFunction
        Initial velocity field u^0
    dirichlet_boundaries : str
        Boundary labels where Dirichlet conditions apply
    velocity_bc : CF
        Dirichlet boundary data for velocity
    dt : float
        Time step size Δt
    T_final : float
        Final simulation time
    viscosity : float, optional
        Kinematic viscosity ν (default: 1.0)
    method : {'picard', 'newton'}, optional
        Linearization method for nonlinear iterations (default: 'picard')
    tolerance : float, optional
        Convergence tolerance for nonlinear iterations (default: 1e-8)
    max_nonlinear_iter : int, optional
        Maximum nonlinear iterations per time step (default: 20)
    save_frequency : int, optional
        Save solution every N time steps (default: 1)
    verbose : bool, optional
        Print progress information (default: True)
    callback : callable, optional
        Function called after each saved time step

    Returns
    -------
    solutions : list of GridFunction
        List of solutions at saved time steps
    times : list of float
        List of corresponding time values

    Examples
    --------
    >>> # Fully implicit time stepping (allows larger time steps)
    >>> solutions, times = time_stepping_implicit(
    ...     mesh, X, u0, "top", u_bc,
    ...     dt=0.05,  # Larger time step than semi-implicit
    ...     T_final=1.0,
    ...     viscosity=0.01,
    ...     method='picard'
    ... )

    Notes
    -----
    Advantages of fully implicit scheme:
    - Unconditionally stable (no CFL restriction)
    - Can use larger time steps
    - Better for stiff problems

    Disadvantages:
    - Requires nonlinear iterations at each time step
    - More expensive per time step than semi-implicit
    """
    from .formulations import unsteady_navier_stokes_implicit_step

    if verbose:
        print("=" * 70)
        print("Fully Implicit Time Stepping for Navier-Stokes")
        print("=" * 70)
        print(f"Time step (Δt): {dt}")
        print(f"Final time (T): {T_final}")
        print(f"Viscosity (ν): {viscosity}")
        print(f"Method: {method}")
        print(f"Nonlinear tolerance: {tolerance:.2e}")
        print(f"Max nonlinear iterations: {max_nonlinear_iter}")
        print(f"Number of steps: {int(T_final / dt)}")
        print("-" * 70)

    # Initialize
    velocity_n = initial_velocity
    time = 0.0
    timestep = 0
    
    solutions = []
    times = []
    
    # Save initial condition
    solution_0 = GridFunction(fespace)
    solution_0.components[0].vec.data = velocity_n.vec
    solutions.append(solution_0)
    times.append(0.0)

    # Time stepping loop
    while time < T_final - 1e-10:
        # Adjust last time step
        if time + dt > T_final:
            dt = T_final - time
        
        # Nonlinear iteration for this time step
        velocity_guess = velocity_n  # Use previous time step as initial guess
        
        for nl_iter in range(max_nonlinear_iter):
            # Perform one linearization step
            solution = unsteady_navier_stokes_implicit_step(
                mesh, fespace, velocity_n, velocity_guess,
                dirichlet_boundaries, velocity_bc,
                dt, viscosity, method
            )
            
            velocity_new = solution.components[0]
            
            # Compute nonlinear iteration error
            error = compute_velocity_error(velocity_new, velocity_guess, mesh)
            
            if verbose and timestep % save_frequency == 0:
                if nl_iter == 0:
                    print(f"Step {timestep:4d} | Time: {time:6.4f}")
                print(f"  NL iter {nl_iter:2d} | Error: {error:.6e}")
            
            # Check convergence
            if error < tolerance:
                break
            
            velocity_guess = velocity_new
        
        if nl_iter == max_nonlinear_iter - 1 and error >= tolerance:
            print(f"WARNING: Nonlinear iteration did not converge at time {time:.4f}")
            print(f"         Final error: {error:.6e}")
        
        # Update for next time step
        velocity_n = velocity_new
        time += dt
        timestep += 1
        
        # Save solution
        if timestep % save_frequency == 0:
            solutions.append(solution)
            times.append(time)
            
            # Call callback
            if callback is not None:
                callback(time, timestep, solution)
    
    if verbose:
        print("-" * 70)
        print(f"✓ Simulation complete: {timestep} time steps")
        print("=" * 70)
    
    return solutions, times