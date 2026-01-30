"""
MMfem: A Python library for solving Navier-Stokes problems using the Finite Element Method.

This package provides tools for:
- Mesh generation and domain discretization
- Finite element space construction (Taylor-Hood, MINI elements)
- Numerical solvers (Picard, Newton iterations)
- Adaptive mesh refinement
- Time-dependent problems

Main modules:
    - mesh: Mesh generation and boundary handling
    - spaces: Finite element spaces (Taylor-Hood, MINI)
    - formulations: Variational formulations (Stokes, Navier-Stokes)
    - solvers: Iterative solvers (Picard, Newton)
    - adaptive: Adaptive mesh refinement strategies
    - time_dependent: Time-stepping schemes for evolutionary problems

Example:
    >>> from mmfem import mesh, spaces, solvers
    >>> domain = mesh.rectangle_mesh(0, 1, 0, 1, h=0.05)
    >>> fem_space = spaces.taylor_hood(domain, "top")
    >>> solution = solvers.picard_iteration(domain, fem_space, ...)

Authors:
    Mauricio Mendiluce

Version:
    0.1.0
"""

__version__ = "0.1.0"
__author__ = "Mauricio Mendiluce"
__email__ = "mmendiluce@dm.uba.ar"

# Import main components for easier access
from .mesh import rectangle_mesh, L_mesh, get_boundary_labels, get_dirichlet_boundaries
from .spaces import taylor_hood, mini_elements, get_space_info
from .formulations import stokes_problem, picard_step, newton_step
from .solvers import picard_iteration, newton_iteration

__all__ = [
    # Mesh functions
    "rectangle_mesh",
    "L_mesh",
    "get_boundary_labels",
    "get_dirichlet_boundaries",
    
    # FEM spaces
    "taylor_hood",
    "mini_elements",
    "get_space_info",
    
    # Problem formulations
    "stokes_problem",
    "picard_step",
    "newton_step",
    
    # Solvers
    "picard_iteration",
    "newton_iteration",
]