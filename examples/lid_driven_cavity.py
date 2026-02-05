"""
Lid-Driven Cavity Flow Example

This script demonstrates how to solve the classic lid-driven cavity problem
using MMfem. The problem consists of a square cavity with a moving top wall.

Problem setup:
    - Domain: [0,1] × [0,1]
    - Boundary conditions:
        * Top wall (y=1): u = (1, 0) (moving lid)
        * Other walls: u = (0, 0) (no-slip)
    - Reynolds number: Re = UL/ν = 100 (ν = 0.01)

This is a standard benchmark in computational fluid dynamics.
"""
from netgen import gui
from ngsolve import CF, Draw
from mmfem import rectangle_mesh, taylor_hood,stabilization_p1p1, picard_iteration,newton_iteration_p1p1
import matplotlib.pyplot as plt


def main():
    """Run the lid-driven cavity simulation."""
    
    print("="*70)
    print("LID-DRIVEN CAVITY FLOW SIMULATION")
    print("="*70)
    
    # Problem parameters
    viscosity = 0.1  # Kinematic viscosity (Re = 100)
    h = 0.02          # Mesh size
    tolerance = 1e-8  # Convergence tolerance
    
    print(f"\nProblem Parameters:")
    print(f"  Viscosity (ν): {viscosity}")
    print(f"  Reynolds number: {1.0 / viscosity}")
    print(f"  Mesh size: {h}")
    print(f"  Tolerance: {tolerance:.2e}")
    
    # Step 1: Create mesh
    print("\n" + "-"*70)
    print("Step 1: Generating mesh...")
    mesh = rectangle_mesh(0, 1, 0, 1, h=h)
    print(f"  Number of elements: {mesh.ne}")
    print(f"  Number of vertices: {mesh.nv}")
    
    # Step 2: Define finite element space
    print("\n" + "-"*70)
    print("Step 2: Setting up finite element space...")
    #X = taylor_hood(mesh, "left|right|bottom|top")
    X = stabilization_p1p1(mesh, "left|right|bottom|top")
    print(f"  Total DOFs: {X.ndof}")
    print(f"  Velocity DOFs: {X.components[0].ndof}")
    print(f"  Pressure DOFs: {X.components[1].ndof}")
    
    # Step 3: Define boundary conditions
    print("\n" + "-"*70)
    print("Step 3: Defining boundary conditions...")
    # Moving lid: horizontal velocity = 1
    u_bc = CF((1, 0))
    print(f"  Dirichlet BC on 'top': u = (1, 0)")
    print(f"  No-slip BC on other walls: u = (0, 0)")
    
    # Step 4: Solve using Picard iteration
    print("\n" + "-"*70)
    print("Step 4: Solving Navier-Stokes equations...")
    
    solution, info = newton_iteration_p1p1(
        mesh=mesh,
        fespace=X,
        dirichlet_boundaries="top",
        velocity_bc=u_bc,
        viscosity=viscosity,
        tolerance=tolerance,
        max_iterations=100,
        verbose=True
    )
    
    # Step 5: Extract results
    print("\n" + "-"*70)
    print("Step 5: Extracting results...")
    velocity, pressure, _ = solution.components
    
    # Compute some statistics
    from ngsolve import Integrate, sqrt, InnerProduct
    
    velocity_magnitude = sqrt(InnerProduct(velocity, velocity))
    max_velocity = velocity_magnitude(mesh(0.5, 1.0))  # At center of top wall
    avg_velocity = Integrate(velocity_magnitude, mesh) / Integrate(1, mesh)
    
    print(f"  Maximum velocity magnitude: {max_velocity:.6f}")
    print(f"  Average velocity magnitude: {avg_velocity:.6f}")
    
    # Step 6: Visualize results
    print("\n" + "-"*70)
    print("Step 6: Visualizing results...")
    
    # Visualize using NGSolve's Draw (opens GUI)
    Draw(velocity, mesh, "velocity")
    Draw(pressure, mesh, "pressure")
    
    # Plot convergence history
    if info['errors']:
        plt.figure(figsize=(10, 6))
        plt.semilogy(info['errors'], 'o-', linewidth=2, markersize=6)
        plt.xlabel('Iteration', fontsize=12)
        plt.ylabel('Relative H1 Error', fontsize=12)
        plt.title('Picard Iteration Convergence', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('convergence_history.png', dpi=150)
        print("  Convergence plot saved to 'convergence_history.png'")
    
    # Final summary
    print("\n" + "="*70)
    print("SIMULATION COMPLETE")
    print("="*70)
    print(f"✓ Converged in {info['iterations']} iterations")
    print(f"✓ Final error: {info['final_error']:.2e}")
    print(f"✓ Solution ready for visualization")
    print("="*70)
    
    # Keep GUI open
    input("\nPress Enter to close visualization and exit...")


if __name__ == "__main__":
    main()