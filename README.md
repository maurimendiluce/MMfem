# MMfem: Finite Element Method for Navier-Stokes Equations

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![NGSolve](https://img.shields.io/badge/NGSolve-required-orange)](https://ngsolve.org/)

A Python library for solving incompressible Navier-Stokes problems using the Finite Element Method (FEM). Built on top of NGSolve, MMfem provides an intuitive API for setting up and solving fluid dynamics problems with various numerical methods.

## Features

- **Multiple FEM Spaces**: Taylor-Hood (P2-P1) and MINI elements
- **Flexible Solvers**: Picard and Newton iterations with customizable convergence criteria
- **Mesh Generation**: Built-in support for rectangular and L-shaped domains
- **Adaptive Refinement**: A posteriori error estimation and automatic mesh refinement
- **Time-Dependent Problems**: Support for evolutionary Navier-Stokes equations (coming soon)
- **Comprehensive Documentation**: Detailed docstrings and examples

## Installation

### Prerequisites

MMfem requires Python 3.10 or later and NGSolve. Install NGSolve first:

```bash
# Using pip
pip install ngsolve

# Or using conda
conda install -c conda-forge ngsolve
```

### Installing MMfem

Clone the repository and install:

```bash
git clone https://github.com/maurimendiluce/MMfem.git
cd MMfem
pip install -e .
```

## Quick Start

Solve the classic lid-driven cavity problem with just a few lines:

```python
from ngsolve import CF
from mmfem import rectangle_mesh, taylor_hood, picard_iteration

# Create mesh
mesh = rectangle_mesh(0, 1, 0, 1, h=0.05)

# Define finite element space
X = taylor_hood(mesh, "left|right|bottom|top")

# Boundary condition: horizontal velocity at top
u_bc = CF((1, 0))

# Solve using Picard iteration
solution, info = picard_iteration(
    mesh, X, 
    dirichlet_boundaries="top",
    velocity_bc=u_bc,
    viscosity=0.01
)

# Extract velocity and pressure
velocity, pressure, _ = solution.components
print(f"Converged in {info['iterations']} iterations")
```

## Usage Examples

### Lid-Driven Cavity

The classic benchmark problem in computational fluid dynamics:

```python
from netgen import gui
from ngsolve import CF, Draw
from mmfem import rectangle_mesh, taylor_hood, picard_iteration

# Setup
mesh = rectangle_mesh(0, 1, 0, 1, h=0.02)
X = taylor_hood(mesh, "left|right|floor")
u_bc = CF((1, 0))  # Unit velocity at top

# Solve
solution, info = picard_iteration(
    mesh, X,
    dirichlet_boundaries="top",
    velocity_bc=u_bc,
    viscosity=0.01,
    tolerance=1e-8
)

# Visualize
velocity = solution.components[0]
Draw(velocity, mesh, "velocity")
```

### Custom Convection Forms

Choose different formulations of the convection term:

```python
# Standard convection: (u·∇)u
solution, _ = picard_iteration(
    mesh, X, "top", u_bc,
    viscosity=0.01,
    convection_form="standard"
)

# Divergence form: (u·∇)u + 0.5(∇·u)u
solution, _ = picard_iteration(
    mesh, X, "top", u_bc,
    viscosity=0.01,
    convection_form="divergence"
)

# Skew-symmetric form: 0.5[(u·∇)u - (∇u)·u]
solution, _ = picard_iteration(
    mesh, X, "top", u_bc,
    viscosity=0.01,
    convection_form="skew_symmetric"
)
```

### Newton Iteration

For faster convergence with a good initial guess:

```python
from mmfem import newton_iteration

solution, info = newton_iteration(
    mesh, X,
    dirichlet_boundaries="top",
    velocity_bc=u_bc,
    viscosity=0.01,
    tolerance=1e-10,
    max_iterations=50
)

print(f"Newton converged in {info['iterations']} iterations")
print(f"Final error: {info['final_error']:.2e}")
```
### Monitoring Convergence

Use callbacks to track convergence:

```python
errors_history = []

def monitor_convergence(iteration, solution, error):
    errors_history.append(error)
    if iteration % 5 == 0:
        print(f"  Iteration {iteration}: error = {error:.2e}")

solution, info = picard_iteration(
    mesh, X, "top", u_bc,
    viscosity=0.01,
    callback=monitor_convergence
)

# Plot convergence
import matplotlib.pyplot as plt
plt.semilogy(errors_history)
plt.xlabel('Iteration')
plt.ylabel('Relative Error')
plt.title('Convergence History')
plt.grid(True)
plt.show()
```

### Adaptive Mesh Refinement

Automatically refine the mesh based on error estimates:

*coming soon*


## API Reference

You can find the documentation [here](docs/build/index.html)
