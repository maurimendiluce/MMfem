"""
Unsteady Lid-Driven Cavity Flow Example

This script demonstrates how to solve the time-dependent lid-driven cavity 
problem using MMfem. The flow starts from rest and evolves towards the 
steady-state solution.

Problem setup:
    - Domain: [0,1] × [0,1]
    - Boundary conditions:
        * Top wall (y=1): u = (1, 0) (moving lid)
        * Other walls: u = (0, 0) (no-slip)
    - Initial condition: u = (0, 0) (flow at rest)
    - Reynolds number: Re = UL/ν = 100 (ν = 0.01)

This demonstrates the transient dynamics of the classic benchmark problem.
"""
from ngsolve import CF, GridFunction, VTKOutput
from mmfem import (
    rectangle_mesh, taylor_hood,
    time_stepping_semiimplicit, time_stepping_implicit,
    export_time_series
)
from mmfem.formulations import stokes_problem
import matplotlib.pyplot as plt
import numpy as np


def main():
    """Run the unsteady lid-driven cavity simulation."""
    
    print("="*70)
    print("UNSTEADY LID-DRIVEN CAVITY FLOW SIMULATION")
    print("="*70)
    
    # Problem parameters
    viscosity = 0.001   # Kinematic viscosity (Re = 1000)
    h = 0.05           # Mesh size
    dt = 0.01          # Time step
    T_final = 10.0      # Final time
    
    print(f"\nProblem Parameters:")
    print(f"  Viscosity (ν): {viscosity}")
    print(f"  Reynolds number: {1.0 / viscosity}")
    print(f"  Mesh size: {h}")
    print(f"  Time step (Δt): {dt}")
    print(f"  Final time: {T_final}")
    print(f"  Number of time steps: {int(T_final / dt)}")
    
    # Step 1: Create mesh
    print("\n" + "-"*70)
    print("Step 1: Generating mesh...")
    mesh = rectangle_mesh(0, 1, 0, 1, h=h)
    print(f"  Number of elements: {mesh.ne}")
    print(f"  Number of vertices: {mesh.nv}")
    
    # Step 2: Define finite element space
    print("\n" + "-"*70)
    print("Step 2: Setting up finite element space...")
    X = taylor_hood(mesh, "left|right|bottom|top")
    print(f"  Total DOFs: {X.ndof}")
    print(f"  Velocity DOFs: {X.components[0].ndof}")
    print(f"  Pressure DOFs: {X.components[1].ndof}")
    
    # Step 3: Define boundary conditions
    print("\n" + "-"*70)
    print("Step 3: Defining boundary conditions...")
    u_bc = CF((1, 0))  # Moving lid
    print(f"  Dirichlet BC on 'top': u = (1, 0)")
    print(f"  No-slip BC on other walls: u = (0, 0)")
    
    # Step 4: Set initial condition
    print("\n" + "-"*70)
    print("Step 4: Setting initial condition...")
    print("  Computing Stokes solution as initial condition...")
    
    # Use Stokes solution as initial condition
    stokes_sol = stokes_problem(mesh, X, "top", u_bc, viscosity)
    u0 = stokes_sol.components[0]
    
    print(f"  Initial condition ready")
    
    # Step 5: Time stepping
    print("\n" + "-"*70)
    print("Step 5: Time stepping simulation...")
    print("\nChoose time stepping method:")
    print("  1. Semi-implicit (fast, requires small time step)")
    print("  2. Fully implicit (stable, allows larger time step)")
    
    choice = input("Enter choice (1 or 2, default=1): ").strip() or "1"
    
    # Storage for kinetic energy
    kinetic_energies = []
    
    def save_energy(time, timestep, solution):
        """Callback to compute and save kinetic energy."""
        from ngsolve import Integrate, InnerProduct
        velocity = solution.components[0]
        KE = 0.5 * Integrate(InnerProduct(velocity, velocity), mesh)
        kinetic_energies.append((time, KE))
    
    if choice == "1":
        print("\nRunning semi-implicit time stepping...")
        solutions, times = time_stepping_semiimplicit(
            mesh=mesh,
            fespace=X,
            initial_velocity=u0,
            dirichlet_boundaries="top",
            velocity_bc=u_bc,
            dt=dt,
            T_final=T_final,
            viscosity=viscosity,
            convection_form="standard",
            save_frequency=10,  # Save every 10 steps
            verbose=True,
            callback=save_energy
        )
    else:
        print("\nRunning fully implicit time stepping...")
        solutions, times = time_stepping_implicit(
            mesh=mesh,
            fespace=X,
            initial_velocity=u0,
            dirichlet_boundaries="top",
            velocity_bc=u_bc,
            dt=0.05,  # Can use larger time step
            T_final=T_final,
            viscosity=viscosity,
            method='picard',
            tolerance=1e-8,
            max_nonlinear_iter=20,
            save_frequency=5,  # Save every 5 steps
            verbose=True,
            callback=save_energy
        )
    
    # Step 6: Post-processing
    print("\n" + "-"*70)
    print("Step 6: Post-processing and visualization...")
    
    # Plot kinetic energy evolution
    if kinetic_energies:
        ke_times, ke_values = zip(*kinetic_energies)
        
        plt.figure(figsize=(10, 6))
        plt.plot(ke_times, ke_values, 'b-', linewidth=2)
        plt.xlabel('Time', fontsize=12)
        plt.ylabel('Kinetic Energy', fontsize=12)
        plt.title('Evolution of Kinetic Energy', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('kinetic_energy.png', dpi=150)
        print("  Kinetic energy plot saved to 'kinetic_energy.png'")
    
    # Export to VTK for visualization in ParaView
    print("\n  Exporting solutions to VTK format...")
    
    pvd_file = export_time_series(
        mesh=mesh,
        solutions=solutions,
        times=times,
        field_names=["velocity", "pressure"],
        output_name="cavity_animation",
        output_dir="vtk_output",
        subdivision=2,
        include_pressure=True,
        verbose=True
    )
    
    # Compute and display statistics at different times
    from ngsolve import Integrate, sqrt, InnerProduct
    
    print("\n  Velocity statistics at selected times:")
    print(f"  {'Time':>8s} | {'Max |u|':>10s} | {'Avg |u|':>10s}")
    print("  " + "-"*35)
    
    for i, (t, sol) in enumerate(zip(times, solutions)):
        if i % max(1, len(times)//5) == 0:  # Show ~5 snapshots
            velocity = sol.components[0]
            velocity_magnitude = sqrt(InnerProduct(velocity, velocity))
            
            # Note: max velocity computation needs evaluation at points
            avg_velocity = Integrate(velocity_magnitude, mesh) / Integrate(1, mesh)
            
            print(f"  {t:8.4f} | {'--':>10s} | {avg_velocity:10.6f}")
    
    # Final summary
    print("\n" + "="*70)
    print("SIMULATION COMPLETE")
    print("="*70)
    print(f"✓ Computed {len(solutions)} time snapshots")
    print(f"✓ Final time reached: {times[-1]:.4f}")
    print(f"✓ Solutions exported to VTK format")
    print(f"✓ View results in ParaView: paraview unsteady_cavity.vtu")
    print("="*70)


if __name__ == "__main__":
    main()