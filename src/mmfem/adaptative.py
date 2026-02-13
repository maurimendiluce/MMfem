"""
Adaptive mesh refinement for Navier-Stokes problems.

This module provides a posteriori error estimation and adaptive mesh refinement
strategies for improving the accuracy of finite element solutions.

Functions:
    estimate_error: Compute a posteriori error estimates
    mark_elements: Mark elements for refinement using Dörfler strategy
    adaptive_solve: Solve with adaptive mesh refinement
    compute_convergence_rates: Analyze convergence rates
    generate_convergence_plot: Visualize convergence behavior
    export_results_latex: Export results as LaTeX table
"""

from typing import Dict, List, Tuple, Optional, Literal
import numpy as np
import matplotlib.pyplot as plt
from ngsolve import (
    Mesh, FESpace, GridFunction, H1, specialcf, grad, dx, Integrate,
    sqrt, InnerProduct, Trace, BND, Draw
)

from .spaces import taylor_hood
from .solvers import newton_iteration


DomainType = Literal["convex","non-convex"]


def estimate_error(
    velocity: GridFunction,
    pressure: GridFunction,
    mesh: Mesh,
    dirichlet_labels: str,
    domain_type: DomainType = "convex"
) -> Tuple[np.ndarray, float]:
    """
    Compute a posteriori error estimates for Navier-Stokes solution.

    Parameters
    ----------
    velocity : GridFunction
        Computed velocity field (vector-valued)
    pressure : GridFunction
        Computed pressure field (scalar)
    mesh : Mesh
        Computational mesh
    domain_type : {'general', 'convex'}, optional
        Type of domain for error estimate scaling (default: 'general')

    Returns
    -------
    eta_local : ndarray
        Element-wise error estimates
    eta_global : float
        Global error estimate

    Examples
    --------
    >>> eta_local, eta_global = estimate_error(u, p, mesh)
    >>> print(f"Global error: {eta_global:.6e}")
    """
    h = specialcf.mesh_size
    n = specialcf.normal(mesh.dim)
    
    V1 = H1(mesh, order=2, autoupdate=True, dirichlet=dirichlet_labels)
    u_x = GridFunction(V1)
    u_y = GridFunction(V1)
    u_x.Set(velocity[0])
    u_y.Set(velocity[1])
    
    laplacian_x = Trace(u_x.Operator("hesse"))
    laplacian_y = Trace(u_y.Operator("hesse"))
    
    p_x = grad(pressure)[0]
    p_y = grad(pressure)[1]
    
    convect_x = u_x * grad(u_x)[0] + u_y * grad(u_x)[1]
    convect_y = u_x * grad(u_y)[0] + u_y * grad(u_y)[1]
    
    residual_x = laplacian_x - p_x - convect_x
    residual_y = laplacian_y - p_y - convect_y
    
    div_u = grad(u_x)[0] + grad(u_y)[1]
    
    jump_x = (grad(u_x) - grad(u_x).Other()) * n
    jump_y = (grad(u_y) - grad(u_y).Other()) * n
    
    if domain_type == "convex":
        alpha, beta, gamma = 8, 4, 5
    elif domain_type == "non-convex":
        alpha, beta, gamma = 20/3, 0, 11/3
    else:
        raise ValueError(f"Invalid domain_type '{domain_type}'")
    
    eta_residual = (h**alpha) * (residual_x**4 + residual_y**4) * dx
    
    if beta > 0:
        eta_divergence = (h**beta) * (div_u**4) * dx
    else:
        eta_divergence = div_u**4 * dx
    
    eta_jump = (
        (h**gamma) * (jump_x**4) * dx(element_vb=BND) +
        (h**gamma) * (jump_y**4) * dx(element_vb=BND)
    )
    
    eta_total = eta_residual + eta_divergence + eta_jump
    eta_local = Integrate(eta_total, mesh, element_wise=True)
    eta_global = np.power(np.sum(eta_local), 0.25)
    
    return eta_local, eta_global


def mark_elements(
    mesh: Mesh,
    eta_local: np.ndarray,
    strategy: Literal["maximum", "dorfler"] = "maximum",
    theta: float = 0.7
) -> int:
    """
    Mark elements for refinement based on error estimates.

    Parameters
    ----------
    mesh : Mesh
        Computational mesh
    eta_local : ndarray
        Element-wise error estimates
    strategy : {'maximum', 'dorfler'}, optional
        Marking strategy (default: 'maximum')
    theta : float, optional
        Marking parameter (default: 0.7)

    Returns
    -------
    n_marked : int
        Number of elements marked

    Examples
    --------
    >>> n_marked = mark_elements(mesh, eta_local, theta=0.7)
    >>> mesh.Refine()
    """
    if not 0 < theta < 1:
        raise ValueError(f"theta must be in (0, 1), got {theta}")
    
    n_marked = 0
    
    if strategy == "maximum":
        eta_max = np.max(eta_local)
        threshold = theta * eta_max
        
        for el in mesh.Elements():
            if eta_local[el.nr] > threshold:
                mesh.SetRefinementFlag(el, True)
                n_marked += 1
    
    elif strategy == "dorfler":
        total_error = np.sum(eta_local)
        target_error = theta * total_error
        
        sorted_indices = np.argsort(eta_local)[::-1]
        
        cumulative_error = 0.0
        for idx in sorted_indices:
            if cumulative_error >= target_error:
                break
            
            mesh.SetRefinementFlag(mesh.Elements()[idx], True)
            cumulative_error += eta_local[idx]
            n_marked += 1
    
    else:
        raise ValueError(f"Invalid strategy '{strategy}'")
    
    return n_marked


def adaptive_solve(
    mesh: Mesh,
    dirichlet_boundaries: str,
    dirichlet_labels: str,
    velocity_bc,
    viscosity: float = 1.0,
    n_refinements: int = 5,
    theta: float = 0.7,
    marking_strategy: Literal["maximum", "dorfler"] = "maximum",
    domain_type: DomainType = "convex",
    tolerance: float = 1e-8,
    verbose: bool = True
) -> Tuple[GridFunction, Dict]:
    """
    Solve Navier-Stokes with adaptive mesh refinement.

    Parameters
    ----------
    mesh : Mesh
        Initial coarse mesh
    dirichlet_boundaries : str
        Boundary labels for Dirichlet conditions
    velocity_bc : CF
        Velocity boundary condition
    viscosity : float, optional
        Kinematic viscosity (default: 1.0)
    n_refinements : int, optional
        Number of refinement cycles (default: 5)
    theta : float, optional
        Marking parameter (default: 0.7)
    marking_strategy : {'maximum', 'dorfler'}, optional
        Marking strategy (default: 'maximum')
    domain_type : {'general', 'convex'}, optional
        Domain type (default: 'general')
    tolerance : float, optional
        Solver tolerance (default: 1e-8)
    verbose : bool, optional
        Print progress (default: True)

    Returns
    -------
    solution : GridFunction
        Final solution
    history : dict
        Refinement history

    Examples
    --------
    >>> solution, history = adaptive_solve(mesh, "left|right|floor", u_bc)
    """
    if verbose:
        print("="*80)
        print("ADAPTIVE MESH REFINEMENT FOR NAVIER-STOKES")
        print("="*80)
        print(f"Initial mesh: {mesh.ne} elements, {mesh.nv} vertices")
        print(f"Viscosity: {viscosity}")
        print(f"Refinement cycles: {n_refinements}")
        print("-"*80)
    
    history = {
        'n_vertices': [],
        'n_elements': [],
        'n_dofs': [],
        'eta_global': [],
        'error_l4': [],
        'n_marked': []
    }
    
    solution_old = None
    
    for cycle in range(n_refinements):
        if verbose:
            print(f"\nCycle {cycle + 1}/{n_refinements}")
            print("-"*80)
        
        X = taylor_hood(mesh, dirichlet_labels)
        
        if verbose:
            print(f"Mesh: {mesh.ne} elements, {mesh.nv} vertices, {X.ndof} DOFs")
            print("Solving...")
        
        solution, solve_info = newton_iteration(
            mesh, X, dirichlet_boundaries, velocity_bc,
            viscosity=viscosity,
            tolerance=tolerance,
            verbose=False
        )
        
        if verbose:
            print(f"  Converged in {solve_info['iterations']} iterations")
        
        velocity = solution.components[0]
        pressure = solution.components[1]
        
        if verbose:
            print("Computing error...")
        
        eta_local, eta_global = estimate_error(
            velocity, pressure, mesh,dirichlet_labels, domain_type=domain_type
        )
        
        if verbose:
            print(f"  Global error: {eta_global:.6e}")
        
        if solution_old is not None:
            velocity_old = solution_old.components[0]
            error_diff = velocity - velocity_old
            error_l4 = np.power(
                Integrate(InnerProduct(error_diff, error_diff)**2, mesh),
                0.25
            )
            history['error_l4'].append(error_l4)
            
            if verbose:
                print(f"  L4 error: {error_l4:.6e}")
        else:
            history['error_l4'].append(np.nan)
        
        history['n_vertices'].append(mesh.nv)
        history['n_elements'].append(mesh.ne)
        history['n_dofs'].append(X.ndof)
        history['eta_global'].append(eta_global)
        
        if cycle < n_refinements - 1:
            if verbose:
                print("Marking elements...")
            
            n_marked = mark_elements(
                mesh, eta_local,
                strategy=marking_strategy,
                theta=theta
            )
            
            history['n_marked'].append(n_marked)
            
            if verbose:
                print(f"  Marked {n_marked} elements ({100*n_marked/mesh.ne:.1f}%)")
                print("Refining...")
            
            solution_old = solution
            mesh.Refine()
        else:
            history['n_marked'].append(0)
    
    history['convergence_rate_eta'] = compute_convergence_rates(
        history['n_elements'], history['eta_global']
    )
    history['convergence_rate_error'] = compute_convergence_rates(
        history['n_elements'][1:], history['error_l4'][1:]
    )
    
    if verbose:
        print("\n" + "="*80)
        print("COMPLETE")
        print("="*80)
        print(f"Final mesh: {mesh.ne} elements")
        print(f"Final error: {history['eta_global'][-1]:.6e}")
        if len(history['convergence_rate_eta']) > 0:
            avg_rate = np.mean(history['convergence_rate_eta'])
            print(f"Avg rate: {avg_rate:.2f}")
        print("="*80)
    
    return solution, history


def compute_convergence_rates(
    n_elements: List[float],
    errors: List[float]
) -> List[float]:
    """Compute convergence rates."""
    rates = []
    
    valid_data = [(ne, e) for ne, e in zip(n_elements, errors) 
                  if not (np.isnan(e) or np.isinf(e)) and e > 0]
    
    if len(valid_data) < 2:
        return rates
    
    for i in range(1, len(valid_data)):
        ne_prev, e_prev = valid_data[i-1]
        ne_curr, e_curr = valid_data[i]
        
        log_ne_diff = np.log(ne_curr) - np.log(ne_prev)
        log_e_diff = np.log(e_curr) - np.log(e_prev)
        
        if abs(log_ne_diff) > 1e-10:
            rate = -log_e_diff / log_ne_diff
            rates.append(rate)
    
    return rates


def generate_convergence_plot(
    history: Dict,
    filename: Optional[str] = None,
    show: bool = True
) -> None:
    """Generate convergence plot."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    n_elements = np.array(history['n_elements'])
    eta_global = np.array(history['eta_global'])
    error_l4 = np.array(history['error_l4'])
    
    ax.loglog(n_elements, eta_global, 'o-', 
              linewidth=2, markersize=8, label='Error estimate (η)')
    
    valid_errors = ~np.isnan(error_l4)
    if np.any(valid_errors):
        ax.loglog(n_elements[valid_errors], error_l4[valid_errors], 's-',
                  linewidth=2, markersize=8, label='L4 error')
    
    if len(history['convergence_rate_eta']) > 0:
        avg_rate_eta = np.mean(history['convergence_rate_eta'])
        x_ref = np.array([n_elements[0], n_elements[-1]])
        y_ref = eta_global[0] * (x_ref / n_elements[0])**(-avg_rate_eta)
        ax.loglog(x_ref, y_ref, 'k--', alpha=0.5, 
                  label=f'Rate ≈ {avg_rate_eta:.2f}')
    
    ax.set_xlabel('Number of Elements', fontsize=12)
    ax.set_ylabel('Error', fontsize=12)
    ax.set_title('Adaptive Refinement Convergence', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"Plot saved to '{filename}'")
    
    if show:
        plt.show()


def export_results_latex(
    history: Dict,
    filename: str = "adaptive_results.tex"
) -> None:
    """Export results as LaTeX table."""
    n_elements = np.array(history['n_elements'])
    eta_global = np.array(history['eta_global'])
    error_l4 = np.array(history['error_l4'])
    rate_eta = history['convergence_rate_eta']
    rate_error = history['convergence_rate_error']
    
    with open(filename, 'w') as f:
        f.write(r"\begin{tabular}{rccccc}" + "\n")
        f.write(r"\hline" + "\n")
        f.write(r"$n_e$ & $\eta$ & rate($\eta$) & $L^4$ error & rate(error) \\" + "\n")
        f.write(r"\hline" + "\n")
        
        for i in range(len(n_elements)):
            ne = int(n_elements[i])
            eta = eta_global[i]
            err = error_l4[i]
            
            r_eta = rate_eta[i-1] if i > 0 and i-1 < len(rate_eta) else 0
            r_err = rate_error[i-1] if i > 0 and i-1 < len(rate_error) else 0
            
            f.write(f"{ne} & {eta:.3e} & ")
            f.write(f"{r_eta:.2f} & " if i > 0 else "-- & ")
            
            if not np.isnan(err):
                f.write(f"{err:.3e} & ")
                f.write(f"{r_err:.2f} \\\\\n" if i > 1 else "-- \\\\\n")
            else:
                f.write("-- & -- \\\\\n")
        
        f.write(r"\hline" + "\n")
        f.write(r"\end{tabular}" + "\n")
    
    print(f"Results exported to '{filename}'")