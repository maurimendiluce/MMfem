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


## Mathematical Background

### Navier-Stokes Equations

MMfem solves the incompressible Navier-Stokes equations:

```
∂u/∂t - ν∆u + (u·∇)u + ∇p = f    in Ω
                     ∇·u = 0    in Ω
                       u = u_D  on ∂Ω_D
          ν∂u/∂n - pn = g      on ∂Ω_N
```

where:
- **u**: velocity field (vector)
- **p**: pressure (scalar)
- **ν**: kinematic viscosity
- **f**: body force
- **Ω**: computational domain

### Weak Formulation

Find (u, p) ∈ V × Q such that:

```
ν(∇u, ∇v) + ((u·∇)u, v) - (p, ∇·v) = (f, v)    ∀v ∈ V
                            (∇·u, q) = 0         ∀q ∈ Q
```

### Discretization Methods

#### Taylor-Hood Elements
- **Velocity**: P_k vector elements (default k=2)
- **Pressure**: P_{k-1} scalar elements
- **Stability**: Inf-sup stable for k ≥ 2

#### MINI Elements
- **Velocity**: P1 + cubic bubble
- **Pressure**: P1
- **Stability**: Inf-sup stable

### Linearization Schemes

#### Picard Iteration
Linearize by evaluating convection at previous iterate:
```
ν∆u^{n+1} - (u^n·∇)u^{n+1} + ∇p^{n+1} = f
```
- **Convergence**: Linear
- **Robustness**: High for low Reynolds numbers

#### Newton Iteration
Use full Jacobian of nonlinear system:
```
J(u^n)[δu] = -R(u^n)
```
- **Convergence**: Quadratic
- **Robustness**: Requires good initial guess

## Project Structure

```
MMfem/
├── mmfem/                   # Main package
│   ├── __init__.py         # Package initialization
│   ├── mesh.py             # Mesh generation
│   ├── spaces.py           # FEM spaces
│   ├── formulations.py     # Variational formulations
│   ├── solvers.py          # Iterative solvers
├── examples/               # Example scripts
│   ├── lid_driven_cavity.py
├── tests/                  # Unit tests
│   ├── test_mesh.py
├── docs/                   # Documentation
│   ├── index.md
├── README.md              # This file
├── requirements.txt       # Dependencies
└── LICENSE                # License file
```

## Contributing

Contributions are welcome! Please follow these guidelines:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feature/amazing-feature`)
5. **Open** a Pull Request

### Code Style

- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guide
- Use type hints for function signatures
- Write comprehensive docstrings (NumPy style)
- Add unit tests for new features

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use MMfem in your research, please cite:

```bibtex
@software{mmfem2026,
  author = {Mendiluce, Mauricio},
  title = {MMfem: Finite Element Method for Navier-Stokes Equations},
  year = {2026},
  url = {https://github.com/maurimendiluce/MMfem}
}
```

## Acknowledgments

- Built on top of [NGSolve](https://ngsolve.org/)
- Inspired by classical FEM literature and best practices
- Thanks to the NGSolve community for excellent documentation

## References

1. Elman, H., Silvester, D., & Wathen, A. (2014). *Finite Elements and Fast Iterative Solvers*. Oxford University Press.
2. Girault, V., & Raviart, P. A. (1986). *Finite Element Methods for Navier-Stokes Equations*. Springer.
3. Braess, D. (2007). *Finite Elements: Theory, Fast Solvers, and Applications*. Cambridge University Press.

## Contact

- GitHub: [@maurimendiluce](https://github.com/maurimendiluce)
- Email: mmendiluce@dm.uba.ar