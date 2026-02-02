"""
VTK export utilities for time-dependent simulations.

This module provides helper functions to export time series data
to VTK format for visualization in ParaView.
"""
import os
from typing import List, Union
from ngsolve import Mesh, GridFunction, VTKOutput


def export_time_series(
    mesh: Mesh,
    solutions: List[GridFunction],
    times: List[float],
    field_names: List[str],
    output_name: str = "simulation",
    output_dir: str = "vtk_output",
    subdivision: int = 2,
    include_pressure: bool = True,
    verbose: bool = True
) -> str:
    """
    Export a time series of solutions to VTK format with PVD animation file.
    
    Creates individual .vtu files for each time step and a .pvd file that
    ParaView can use to display the animation.
    
    Parameters
    ----------
    mesh : Mesh
        The computational mesh
    solutions : List[GridFunction]
        List of solution GridFunctions at different time steps.
        Each GridFunction should have components [velocity, pressure, lambda].
    times : List[float]
        Time value for each solution
    field_names : List[str]
        Names of fields to export (e.g., ["velocity", "pressure"])
    output_name : str, optional
        Base name for output files (default: "simulation")
    output_dir : str, optional
        Directory to store VTK files (default: "vtk_output")
    subdivision : int, optional
        Subdivision level for output (default: 2)
    include_pressure : bool, optional
        Whether to include pressure field (default: True)
    verbose : bool, optional
        Print progress information (default: True)
    
    Returns
    -------
    str
        Path to the .pvd animation file
    
    Examples
    --------
    >>> solutions, times = time_stepping_semiimplicit(...)
    >>> pvd_file = export_time_series(
    ...     mesh, solutions, times,
    ...     field_names=["velocity", "pressure"],
    ...     output_name="cavity_flow"
    ... )
    >>> print(f"Open with: paraview {pvd_file}")
    
    Notes
    -----
    The .pvd file format is XML-based and tells ParaView which .vtu files
    correspond to which time steps. This enables:
    - Animation playback
    - Time slider control
    - Temporal queries
    
    ParaView usage:
    1. Open the .pvd file in ParaView
    2. Click "Apply" in the Properties panel
    3. Use the animation controls to play/pause
    4. Use the time slider to navigate
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    if verbose:
        print(f"\nExporting {len(solutions)} time steps to VTK...")
    
    # Export each time step
    for i, (sol, t) in enumerate(zip(solutions, times)):
        # Extract fields
        coefs = []
        names = []
        
        # Velocity (always included)
        velocity = sol.components[0]
        coefs.append(velocity)
        names.append("velocity")
        
        # Pressure (optional)
        if include_pressure and len(sol.components) > 1:
            pressure = sol.components[1]
            coefs.append(pressure)
            names.append("pressure")
        
        # Create VTK file for this time step
        vtk_file = f"{output_dir}/{output_name}_{i:04d}"
        vtk = VTKOutput(
            ma=mesh,
            coefs=coefs,
            names=names,
            filename=vtk_file,
            subdivision=subdivision
        )
        vtk.Do()
        
        if verbose and (i + 1) % max(1, len(solutions) // 10) == 0:
            print(f"  Progress: {i+1}/{len(solutions)} files exported")
    
    # Create PVD collection file
    pvd_filename = f"{output_name}.pvd"
    pvd_content = '<?xml version="1.0"?>\n'
    pvd_content += '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
    pvd_content += '  <Collection>\n'
    
    for i, t in enumerate(times):
        vtu_file = f"{output_dir}/{output_name}_{i:04d}.vtu"
        pvd_content += f'    <DataSet timestep="{t}" file="{vtu_file}"/>\n'
    
    pvd_content += '  </Collection>\n'
    pvd_content += '</VTKFile>\n'
    
    with open(pvd_filename, "w") as f:
        f.write(pvd_content)
    
    if verbose:
        print(f"\n✓ Export complete:")
        print(f"  • {len(solutions)} VTK files in '{output_dir}/'")
        print(f"  • Animation file: '{pvd_filename}'")
        print(f"\n📽️  To view animation:")
        print(f"     paraview {pvd_filename}")
    
    return pvd_filename


def export_single_timestep(
    mesh: Mesh,
    solution: GridFunction,
    filename: str,
    subdivision: int = 2,
    include_pressure: bool = True
) -> None:
    """
    Export a single solution to VTK format.
    
    Parameters
    ----------
    mesh : Mesh
        The computational mesh
    solution : GridFunction
        Solution with components [velocity, pressure, lambda]
    filename : str
        Output filename (without .vtu extension)
    subdivision : int, optional
        Subdivision level (default: 2)
    include_pressure : bool, optional
        Include pressure field (default: True)
    
    Examples
    --------
    >>> export_single_timestep(mesh, final_solution, "final_state")
    """
    coefs = [solution.components[0]]  # velocity
    names = ["velocity"]
    
    if include_pressure and len(solution.components) > 1:
        coefs.append(solution.components[1])
        names.append("pressure")
    
    vtk = VTKOutput(
        ma=mesh,
        coefs=coefs,
        names=names,
        filename=filename,
        subdivision=subdivision
    )
    vtk.Do()


def export_comparison(
    mesh: Mesh,
    solutions_dict: dict,
    times: List[float],
    output_name: str = "comparison",
    subdivision: int = 2
) -> List[str]:
    """
    Export multiple solution series for comparison.
    
    Useful for comparing different methods or parameters.
    
    Parameters
    ----------
    mesh : Mesh
        The computational mesh
    solutions_dict : dict
        Dictionary mapping method names to solution lists
        Example: {"semi_implicit": solutions1, "fully_implicit": solutions2}
    times : List[float]
        Time values (assumed same for all methods)
    output_name : str, optional
        Base name for outputs (default: "comparison")
    subdivision : int, optional
        Subdivision level (default: 2)
    
    Returns
    -------
    List[str]
        List of .pvd filenames created
    
    Examples
    --------
    >>> pvd_files = export_comparison(
    ...     mesh,
    ...     {
    ...         "semi_implicit": solutions1,
    ...         "fully_implicit": solutions2
    ...     },
    ...     times
    ... )
    >>> for pvd in pvd_files:
    ...     print(f"Created: {pvd}")
    """
    pvd_files = []
    
    for method_name, solutions in solutions_dict.items():
        pvd_file = export_time_series(
            mesh=mesh,
            solutions=solutions,
            times=times,
            field_names=["velocity", "pressure"],
            output_name=f"{output_name}_{method_name}",
            output_dir=f"vtk_{method_name}",
            subdivision=subdivision,
            verbose=True
        )
        pvd_files.append(pvd_file)
    
    print(f"\n✓ Comparison export complete: {len(pvd_files)} methods")
    return pvd_files


# Example usage
if __name__ == "__main__":
    print("VTK Export Utilities")
    print("=" * 50)
    print("\nUsage example:")
    print("""
    from vtk_utils import export_time_series
    
    # After running time stepping:
    solutions, times = time_stepping_semiimplicit(...)
    
    # Export for ParaView animation:
    pvd_file = export_time_series(
        mesh, solutions, times,
        field_names=["velocity", "pressure"],
        output_name="my_simulation"
    )
    
    # Open in ParaView:
    # $ paraview my_simulation.pvd
    """)