"""
Simple Unsteady Navier-Stokes Example

This example demonstrates the basic usage of the time-dependent Navier-Stokes
solvers without all the bells and whistles. It's a minimal working example.

Problem: Flow in a rectangular domain with a parabolic inflow profile
"""
from ngsolve import CF, x, y, Integrate, InnerProduct, sqrt
from mmfem import (
    rectangle_mesh, taylor_hood,
    stokes_problem,
    time_stepping_semiimplicit, time_stepping_implicit,
    export_time_series
)
from mmfem.formulations import (
    unsteady_navier_stokes_semiimplicit_step,
    unsteady_navier_stokes_implicit_step
)
from mmfem.solvers import compute_velocity_error


def example_semi_implicit():
    """Example using semi-implicit time stepping."""
    print("\n" + "="*60)
    print("SEMI-IMPLICIT TIME STEPPING EXAMPLE")
    print("="*60)
    
    # Setup
    mesh = rectangle_mesh(0, 2, 0, 1, h=0.1)
    X = taylor_hood(mesh, "left|right|bottom|top")
    
    # Parabolic inflow at left boundary
    u_in = CF((4 * y * (1 - y), 0))
    
    # Initial condition: Stokes solution
    print("\nComputing initial Stokes solution...")
    stokes_sol = stokes_problem(mesh, X, "left", u_in, viscosity=0.01)
    velocity_n = stokes_sol.components[0]
    
    # Time parameters
    dt = 0.01
    n_steps = 50
    
    print(f"\nTime stepping:")
    print(f"  Time step: {dt}")
    print(f"  Number of steps: {n_steps}")
    print(f"  Final time: {dt * n_steps}")
    
    # Time loop
    print("\nProgress:")
    for n in range(n_steps):
        # Take one time step
        solution = unsteady_navier_stokes_semiimplicit_step(
            mesh=mesh,
            fespace=X,
            velocity_prev=velocity_n,
            dirichlet_boundaries="left",
            velocity_bc=u_in,
            dt=dt,
            viscosity=0.01
        )
        
        velocity_new = solution.components[0]
        
        # Compute change from previous step
        if n > 0:
            error = compute_velocity_error(velocity_new, velocity_n, mesh)
        else:
            error = 0.0
        
        # Update for next step
        velocity_n = velocity_new
        
        if (n + 1) % 10 == 0:
            time = (n + 1) * dt
            ke = 0.5 * Integrate(InnerProduct(velocity_n, velocity_n), mesh)
            print(f"  Step {n+1:3d} | Time: {time:.3f} | KE: {ke:.6f} | Change: {error:.2e}")
    
    print("\n✓ Semi-implicit simulation complete")
    return velocity_n


def example_fully_implicit():
    """Example using fully implicit time stepping with Picard iterations."""
    print("\n" + "="*60)
    print("FULLY IMPLICIT TIME STEPPING EXAMPLE")
    print("="*60)
    
    # Setup
    mesh = rectangle_mesh(0, 2, 0, 1, h=0.1)
    X = taylor_hood(mesh, "left|right|bottom|top")
    
    # Parabolic inflow
    u_in = CF((4 * y * (1 - y), 0))
    
    # Initial condition
    print("\nComputing initial Stokes solution...")
    stokes_sol = stokes_problem(mesh, X, "left", u_in, viscosity=0.01)
    velocity_n = stokes_sol.components[0]
    
    # Time parameters
    dt = 0.05  # Larger time step than semi-implicit!
    n_steps = 20
    
    print(f"\nTime stepping:")
    print(f"  Time step: {dt} (larger than semi-implicit!)")
    print(f"  Number of steps: {n_steps}")
    print(f"  Final time: {dt * n_steps}")
    
    # Time loop
    print("\nProgress:")
    for n in range(n_steps):
        # Nonlinear iteration at this time step
        velocity_guess = velocity_n  # Use previous time as initial guess
        
        max_nl_iter = 10
        nl_tol = 1e-8
        
        for nl_iter in range(max_nl_iter):
            # One Picard iteration
            solution = unsteady_navier_stokes_implicit_step(
                mesh=mesh,
                fespace=X,
                velocity_prev=velocity_n,
                velocity_guess=velocity_guess,
                dirichlet_boundaries="left",
                velocity_bc=u_in,
                dt=dt,
                viscosity=0.01,
                method='picard'
            )
            
            velocity_new = solution.components[0]
            
            # Check nonlinear convergence
            nl_error = compute_velocity_error(velocity_new, velocity_guess, mesh)
            
            if nl_error < nl_tol:
                break
            
            velocity_guess = velocity_new
        
        # Update for next time step
        velocity_n = velocity_new
        time = (n + 1) * dt
        ke = 0.5 * Integrate(InnerProduct(velocity_n, velocity_n), mesh)
        
        print(f"  Step {n+1:3d} | Time: {time:.3f} | NL iters: {nl_iter+1:2d} | KE: {ke:.6f}")
    
    print("\n✓ Fully implicit simulation complete")
    return velocity_n


def main():
    """Run both examples."""
    print("\n" + "="*60)
    print("SIMPLE UNSTEADY NAVIER-STOKES EXAMPLES")
    print("="*60)
    print("\nThis script demonstrates two approaches:")
    print("  1. Semi-implicit: Fast but requires small time steps")
    print("  2. Fully implicit: More stable, allows larger time steps")
    
    # Run semi-implicit example
    u_semi = example_semi_implicit()
    
    # Run fully implicit example
    u_implicit = example_fully_implicit()
    
    print("\n" + "="*60)
    print("BOTH SIMULATIONS COMPLETE")
    print("="*60)
    print("\nKey observations:")
    print("  • Semi-implicit used dt=0.01 (50 steps to reach t=0.5)")
    print("  • Fully implicit used dt=0.05 (20 steps to reach t=1.0)")
    print("  • Fully implicit requires nonlinear iterations per step")
    print("  • Trade-off: fewer steps vs. more work per step")
    
    # Optional: Export for visualization
    print("\n" + "="*60)
    print("OPTIONAL: Export to VTK for ParaView")
    print("="*60)
    print("\nWould you like to export results for visualization?")
    export_choice = input("Export to VTK? (y/n, default=n): ").strip().lower()
    
    if export_choice == 'y':
        print("\nExporting semi-implicit results...")
        # Note: We'd need to store solutions from the examples above
        # For now, just show how it would work:
        print("  (Re-run semi-implicit to get solutions)")
        
        mesh = rectangle_mesh(0, 2, 0, 1, h=0.1)
        X = taylor_hood(mesh, "left|right|bottom|top")
        u_in = CF((4 * y * (1 - y), 0))
        
        stokes_sol = stokes_problem(mesh, X, "left", u_in, viscosity=0.01)
        u0 = stokes_sol.components[0]
        
        solutions_export, times_export = time_stepping_semiimplicit(
            mesh=mesh,
            fespace=X,
            initial_velocity=u0,
            dirichlet_boundaries="left",
            velocity_bc=u_in,
            dt=0.05,
            T_final=0.5,
            viscosity=0.01,
            save_frequency=2,
            verbose=False
        )
        
        pvd_file = export_time_series(
            mesh=mesh,
            solutions=solutions_export,
            times=times_export,
            field_names=["velocity", "pressure"],
            output_name="simple_example",
            verbose=True
        )
        
        print(f"\n✓ Export complete! Open with: paraview {pvd_file}")
    
    print("="*60)


if __name__ == "__main__":
    main()