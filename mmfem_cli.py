"""
MMfem Interactive Command Line Interface

Interactive tool for solving Navier-Stokes problems without writing code.
Supports both steady-state and time-dependent (unsteady) problems.
"""

import sys
from pathlib import Path
from typing import Optional
from ngsolve import CF, Draw


def clear_screen():
    """Clear the terminal screen."""
    import os
    os.system('cls' if os.name == 'nt' else 'clear')


def print_header():
    """Print application header."""
    print("="*80)
    print(" "*25 + "MMfem - Interactive CLI")
    print(" "*15 + "Navier-Stokes Finite Element Solver")
    print(" "*20 + "(Steady & Unsteady Problems)")
    print("="*80)


def print_menu(title: str, options: dict):
    """Print a menu with options."""
    print(f"\n{title}")
    print("-"*80)
    for key, description in options.items():
        print(f"  [{key}] {description}")
    print("-"*80)


def get_choice(prompt: str, valid_choices: list) -> str:
    """Get validated user choice."""
    while True:
        choice = input(f"\n{prompt}: ").strip().lower()
        if choice in valid_choices:
            return choice
        print(f"❌ Opción inválida. Elige una de: {', '.join(valid_choices)}")


def get_float(prompt: str, default: Optional[float] = None, min_val: Optional[float] = None) -> float:
    """Get validated float input."""
    while True:
        default_text = f" (default: {default})" if default is not None else ""
        user_input = input(f"\n{prompt}{default_text}: ").strip()
        
        if not user_input and default is not None:
            return default
        
        try:
            value = float(user_input)
            if min_val is not None and value < min_val:
                print(f"❌ El valor debe ser mayor o igual a {min_val}")
                continue
            return value
        except ValueError:
            print("❌ Por favor ingresa un número válido")


def get_int(prompt: str, default: Optional[int] = None, min_val: Optional[int] = None) -> int:
    """Get validated integer input."""
    while True:
        default_text = f" (default: {default})" if default is not None else ""
        user_input = input(f"\n{prompt}{default_text}: ").strip()
        
        if not user_input and default is not None:
            return default
        
        try:
            value = int(user_input)
            if min_val is not None and value < min_val:
                print(f"❌ El valor debe ser mayor o igual a {min_val}")
                continue
            return value
        except ValueError:
            print("❌ Por favor ingresa un número entero válido")


def setup_problem():
    """Interactive problem setup."""
    clear_screen()
    print_header()
    
    print("\n" + "="*80)
    print("CONFIGURACIÓN DEL PROBLEMA")
    print("="*80)
    
    # 0. Problem type (steady vs unsteady)
    print_menu("0. Tipo de problema:", {
        "1": "Estacionario (steady-state)",
        "2": "Evolutivo (time-dependent)"
    })
    problem_type = get_choice("Opción", ["1", "2"])
    is_unsteady = (problem_type == "2")
    
    # 1. Domain selection
    print_menu("1. Selecciona el dominio:", {
        "1": "Rectángulo (cavity flow)",
        "2": "Dominio en L (con singularidad)"
    })
    domain_type = get_choice("Opción", ["1", "2"])
    
    # 2. Mesh size
    print("\n2. Tamaño de malla")
    if domain_type == "1":
        h = get_float("Tamaño de malla h", default=0.05, min_val=0.01)
    else:
        h = get_float("Tamaño de malla h", default=0.3, min_val=0.05)
    
    # 3. Viscosity
    print("\n3. Viscosidad")
    viscosity = get_float("Viscosidad (ν)", default=0.01, min_val=0.0001)
    reynolds = 1.0 / viscosity
    print(f"   → Número de Reynolds: Re ≈ {reynolds:.0f}")
    
    # 4. Time parameters (only for unsteady)
    if is_unsteady:
        print("\n4. Parámetros temporales")
        dt = get_float("Paso de tiempo (Δt)", default=0.01, min_val=0.0001)
        T_final = get_float("Tiempo final (T)", default=1.0, min_val=0.01)
        n_steps = int(T_final / dt)
        print(f"   → Número de pasos de tiempo: {n_steps}")
        
        save_freq = get_int("Guardar cada N pasos", default=10, min_val=1)
    
    # 5. Solver selection
    if is_unsteady:
        print_menu("5. Selecciona el esquema temporal:", {
            "1": "Semi-implícito (rápido, requiere dt pequeño)",
            "2": "Totalmente implícito con Picard (estable, dt grande)",
            "3": "Totalmente implícito con Newton (convergencia rápida)"
        })
        solver_type = get_choice("Opción", ["1", "2", "3"])
        
        if solver_type in ["2", "3"]:
            # Nonlinear iteration parameters for implicit methods
            nl_tol = get_float("Tolerancia no lineal", default=1e-8, min_val=1e-15)
            max_nl_iter = get_int("Máx. iteraciones no lineales por paso", default=20, min_val=1)
        else:
            nl_tol = None
            max_nl_iter = None
        
        # Convection form (only for semi-implicit)
        if solver_type == "1":
            print_menu("Forma del término convectivo:", {
                "1": "Standard: (u·∇)u",
                "2": "Divergence: (u·∇)u + 0.5(∇·u)u",
                "3": "Skew-symmetric: 0.5[(u·∇)u - (∇u)·u]"
            })
            conv_choice = get_choice("Opción", ["1", "2", "3"])
            conv_forms = {"1": "standard", "2": "divergence", "3": "skew_symmetric"}
            convection_form = conv_forms[conv_choice]
        else:
            convection_form = "standard"
    else:
        # Steady-state solver
        print_menu("5. Selecciona el método de solución:", {
            "1": "Picard (robusto, convergencia lineal)",
            "2": "Newton (rápido, convergencia cuadrática)"
        })
        solver_type = get_choice("Opción", ["1", "2"])
        
        # Solver parameters
        tolerance = get_float("Tolerancia de convergencia", default=1e-8, min_val=1e-15)
        max_iter = get_int("Máximo de iteraciones", default=100 if solver_type == "1" else 50, min_val=1)
        
        if solver_type == "1":
            print_menu("Forma del término convectivo:", {
                "1": "Standard: (u·∇)u",
                "2": "Divergence: (u·∇)u + 0.5(∇·u)u",
                "3": "Skew-symmetric: 0.5[(u·∇)u - (∇u)·u]"
            })
            conv_choice = get_choice("Opción", ["1", "2", "3"])
            conv_forms = {"1": "standard", "2": "divergence", "3": "skew_symmetric"}
            convection_form = conv_forms[conv_choice]
        else:
            convection_form = None
    
    # 6. Visualization/Export options
    print_menu("6. Opciones de visualización/exportación:", {
        "1": "Solo visualizar (ventana NGSolve)",
        "2": "Exportar para ParaView (.pvd)" + (" + visualizar" if not is_unsteady else ""),
        "3": "Ambos" if is_unsteady else "Sin visualización"
    })
    viz_choice = get_choice("Opción", ["1", "2", "3"])
    
    # Summary
    print("\n" + "="*80)
    print("RESUMEN DE CONFIGURACIÓN")
    print("="*80)
    print(f"  Tipo: {'Evolutivo (time-dependent)' if is_unsteady else 'Estacionario (steady-state)'}")
    print(f"  Dominio: {'Rectángulo' if domain_type == '1' else 'L-shaped'}")
    print(f"  Tamaño de malla: h = {h}")
    print(f"  Viscosidad: ν = {viscosity} (Re ≈ {reynolds:.0f})")
    
    if is_unsteady:
        print(f"  Paso de tiempo: Δt = {dt}")
        print(f"  Tiempo final: T = {T_final}")
        print(f"  Pasos de tiempo: {n_steps}")
        print(f"  Guardar cada: {save_freq} pasos")
        
        if solver_type == "1":
            print(f"  Esquema: Semi-implícito (IMEX)")
            print(f"  Convección: {convection_form}")
        elif solver_type == "2":
            print(f"  Esquema: Totalmente implícito (Picard)")
            print(f"  Tolerancia NL: {nl_tol:.2e}")
            print(f"  Max iter NL: {max_nl_iter}")
        else:
            print(f"  Esquema: Totalmente implícito (Newton)")
            print(f"  Tolerancia NL: {nl_tol:.2e}")
            print(f"  Max iter NL: {max_nl_iter}")
    else:
        if solver_type == "1":
            print(f"  Método: Picard")
            print(f"  Convección: {convection_form}")
            print(f"  Tolerancia: {tolerance:.2e}")
            print(f"  Max iteraciones: {max_iter}")
        elif solver_type == "2":
            print(f"  Método: Newton")
            print(f"  Tolerancia: {tolerance:.2e}")
            print(f"  Max iteraciones: {max_iter}")
    
    print("="*80)
    
    confirm = get_choice("¿Proceder con esta configuración? (s/n)", ["s", "n", "y"])
    
    if confirm == "n":
        return None
    
    # Build configuration
    config = {
        'is_unsteady': is_unsteady,
        'domain_type': domain_type,
        'h': h,
        'viscosity': viscosity,
        'solver_type': solver_type,
        'convection_form': convection_form,
        'viz_choice': viz_choice
    }
    
    if is_unsteady:
        config['dt'] = dt
        config['T_final'] = T_final
        config['save_frequency'] = save_freq
        if solver_type in ["2", "3"]:
            config['nl_tolerance'] = nl_tol
            config['max_nl_iter'] = max_nl_iter
    else:
        config['tolerance'] = tolerance
        config['max_iter'] = max_iter
    
    return config


def solve_problem(config):
    """Solve the problem with given configuration."""
    print("\n" + "="*80)
    print("RESOLVIENDO...")
    print("="*80)
    
    try:
        from mmfem import (
            rectangle_mesh, L_mesh,
            taylor_hood,
            picard_iteration, newton_iteration,
            stokes_problem,
            time_stepping_semiimplicit, time_stepping_implicit,
            export_time_series
        )
        
        # Create mesh
        print("\n1. Generando malla...")
        if config['domain_type'] == "1":
            mesh = rectangle_mesh(0, 1, 0, 1, h=config['h'])
            dirichlet_boundaries = "left|right|bottom|top"
            bc_label = "top"
        else:
            mesh = L_mesh(h=config['h'])
            dirichlet_boundaries = "side"
            bc_label = "flow_in"
        
        print(f"   ✓ Malla creada: {mesh.ne} elementos, {mesh.nv} vértices")
        
        # Create FEM space
        X = taylor_hood(mesh, dirichlet_boundaries)
        print(f"   ✓ Espacio FEM: {X.ndof} DOFs")
        
        # Boundary condition
        u_bc = CF((1, 0))
        
        # Solve
        if config['is_unsteady']:
            # UNSTEADY PROBLEM
            print("\n2. Resolviendo problema evolutivo de Navier-Stokes...")
            
            # Initial condition: Stokes solution
            print("   Calculando condición inicial (Stokes)...")
            stokes_sol = stokes_problem(mesh, X, bc_label, u_bc, config['viscosity'])
            u0 = stokes_sol.components[0]
            print("   ✓ Condición inicial lista")
            
            # Time stepping
            if config['solver_type'] == "1":
                # Semi-implicit
                print(f"\n   Esquema semi-implícito: {int(config['T_final']/config['dt'])} pasos...")
                solutions, times = time_stepping_semiimplicit(
                    mesh=mesh,
                    fespace=X,
                    initial_velocity=u0,
                    dirichlet_boundaries=bc_label,
                    velocity_bc=u_bc,
                    dt=config['dt'],
                    T_final=config['T_final'],
                    viscosity=config['viscosity'],
                    convection_form=config['convection_form'],
                    save_frequency=config['save_frequency'],
                    verbose=True
                )
                info = {'method': 'semi-implicit', 'n_snapshots': len(solutions)}
                
            else:
                # Fully implicit (Picard or Newton)
                method = 'picard' if config['solver_type'] == "2" else 'newton'
                print(f"\n   Esquema totalmente implícito ({method})...")
                solutions, times = time_stepping_implicit(
                    mesh=mesh,
                    fespace=X,
                    initial_velocity=u0,
                    dirichlet_boundaries=bc_label,
                    velocity_bc=u_bc,
                    dt=config['dt'],
                    T_final=config['T_final'],
                    viscosity=config['viscosity'],
                    method=method,
                    tolerance=config['nl_tolerance'],
                    max_nonlinear_iter=config['max_nl_iter'],
                    save_frequency=config['save_frequency'],
                    verbose=True
                )
                info = {'method': f'implicit-{method}', 'n_snapshots': len(solutions)}
            
            solution = solutions[-1]  # Final solution for visualization
            
        else:
            # STEADY-STATE PROBLEM
            print("\n2. Resolviendo ecuaciones de Navier-Stokes estacionarias...")
            
            if config['solver_type'] == "1":
                # Picard
                solution, info = picard_iteration(
                    mesh, X,
                    bc_label, u_bc,
                    viscosity=config['viscosity'],
                    convection_form=config['convection_form'],
                    tolerance=config['tolerance'],
                    max_iterations=config['max_iter'],
                    verbose=True
                )
                
            else:
                # Newton
                solution, info = newton_iteration(
                    mesh, X,
                    bc_label, u_bc,
                    viscosity=config['viscosity'],
                    tolerance=config['tolerance'],
                    max_iterations=config['max_iter'],
                    verbose=True
                )
            
            solutions = None
            times = None
        
        # Results
        print("\n" + "="*80)
        print("SOLUCIÓN ENCONTRADA")
        print("="*80)
        
        if config['is_unsteady']:
            print(f"  Snapshots guardados: {info['n_snapshots']}")
            print(f"  Tiempo final: {times[-1]:.4f}")
            print(f"  Método: {info['method']}")
        else:
            print(f"  Convergió: {'✓ Sí' if info['converged'] else '✗ No'}")
            print(f"  Iteraciones: {info['iterations']}")
            print(f"  Error final: {info['final_error']:.6e}")
        
        print("="*80)
        
        return solution, info, mesh, solutions, times
        
    except Exception as e:
        print(f"\n❌ Error durante la solución: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None, None, None


def visualize_results(solution, info, mesh, solutions, times, config):
    """Visualize and save results."""
    print("\n" + "="*80)
    print("VISUALIZACIÓN Y EXPORTACIÓN")
    print("="*80)
    
    viz_choice = config['viz_choice']
    
    # Visualize in NGSolve GUI
    if viz_choice in ["1", "3"]:
        print("\nAbriendo ventanas de visualización NGSolve...")
        velocity = solution.components[0]
        pressure = solution.components[1]
        
        Draw(velocity, mesh, "Velocidad")
        Draw(pressure, mesh, "Presión")
        print("  ✓ Ventanas abiertas (ciérralas para continuar)")
    
    # Export to VTK/ParaView
    if viz_choice in ["2", "3"]:
        if config['is_unsteady']:
            # Export time series
            print("\nExportando serie temporal para ParaView...")
            try:
                from mmfem import export_time_series
                
                output_name = input("\nNombre del archivo (default: simulation): ").strip() or "simulation"
                
                pvd_file = export_time_series(
                    mesh=mesh,
                    solutions=solutions,
                    times=times,
                    field_names=["velocity", "pressure"],
                    output_name=output_name,
                    output_dir="vtk_output",
                    subdivision=2,
                    verbose=True
                )
                
                print(f"\n  🎬 Para ver la animación:")
                print(f"     paraview {pvd_file}")
                
            except Exception as e:
                print(f"  ❌ Error exportando: {e}")
        else:
            # Export single snapshot
            print("\nExportando para ParaView...")
            try:
                from ngsolve import VTKOutput
                
                output_name = input("\nNombre del archivo (default: steady_solution): ").strip() or "steady_solution"
                
                velocity = solution.components[0]
                pressure = solution.components[1]
                
                vtk = VTKOutput(
                    ma=mesh,
                    coefs=[velocity, pressure],
                    names=["velocity", "pressure"],
                    filename=output_name,
                    subdivision=2
                )
                vtk.Do()
                
                print(f"  ✓ Archivo guardado: {output_name}.vtu")
                print(f"  📊 Abrir con: paraview {output_name}.vtu")
                
            except Exception as e:
                print(f"  ❌ Error exportando: {e}")


def main():
    """Main interactive loop."""
    while True:
        clear_screen()
        print_header()
        
        print_menu("MENÚ PRINCIPAL", {
            "1": "Nuevo problema",
            "2": "Ayuda / Guía rápida",
            "3": "Acerca de",
            "4": "Salir"
        })
        
        choice = get_choice("Opción", ["1", "2", "3", "4"])
        
        if choice == "1":
            # Solve problem
            config = setup_problem()
            if config is None:
                input("\nPresiona Enter para volver al menú...")
                continue
            
            solution, info, mesh, solutions, times = solve_problem(config)
            
            if solution is not None:
                visualize_results(solution, info, mesh, solutions, times, config)
            
            input("\nPresiona Enter para volver al menú...")
        
        elif choice == "2":
            # Help
            clear_screen()
            print_header()
            print("\n" + "="*80)
            print("GUÍA RÁPIDA")
            print("="*80)
            print("""
MMfem resuelve las ecuaciones de Navier-Stokes usando elementos finitos.

TIPOS DE PROBLEMA:

1. ESTACIONARIO (Steady-state):
   - Busca la solución en equilibrio (∂u/∂t = 0)
   - Más rápido de resolver
   - Bueno para entender el comportamiento final

2. EVOLUTIVO (Time-dependent):
   - Simula la evolución temporal del flujo
   - Muestra cómo el flujo evoluciona desde condición inicial
   - Permite estudiar transitorios y dinámica
   - Genera animaciones visualizables en ParaView

ESQUEMAS TEMPORALES (solo evolutivo):

1. Semi-implícito (IMEX):
   - Rápido (un sistema lineal por paso)
   - Requiere Δt pequeño (condición CFL)
   - Recomendado: Δt ~ h/|u|_max

2. Totalmente implícito (Picard/Newton):
   - Estable (sin condición CFL)
   - Permite Δt más grande (5-10x)
   - Requiere iteraciones no lineales

PARÁMETROS CLAVE:

- Viscosidad (ν): Controla Re = 1/ν
  * 0.01 → Re = 100 (flujo laminar estable)
  * 0.001 → Re = 1000 (transición)

- Paso de tiempo (Δt, solo evolutivo):
  * Semi-implícito: 0.001-0.02
  * Totalmente implícito: 0.01-0.1

- Tiempo final (T, solo evolutivo):
  * 1.0-5.0 para ver desarrollo completo

VISUALIZACIÓN:

- NGSolve: Ventanas interactivas inmediatas
- ParaView: Exporta .pvd para animaciones profesionales
  → Abre con: paraview archivo.pvd
  → Click "Apply", selecciona campo, Play ▶️

RECOMENDACIONES:

Problema estacionario:
  ✓ Comenzar con Picard
  ✓ Re < 200 para garantizar convergencia
  ✓ Usar Newton si Picard es lento

Problema evolutivo:
  ✓ Comenzar con semi-implícito
  ✓ Δt = 0.01, T = 1.0 como prueba
  ✓ Usar implícito para Re alto o Δt grande
  ✓ save_frequency = 5-10 para no saturar disco

Para más información: 
  - UNSTEADY_IMPLEMENTATION.md
  - PARAVIEW_GUIDE.md
  - VTK_SOLUTION.md
            """)
            input("\nPresiona Enter para continuar...")
        
        elif choice == "3":
            # About
            clear_screen()
            print_header()
            print("\n" + "="*80)
            print("ACERCA DE MMfem")
            print("="*80)
            print("""
MMfem v0.2.0
Librería Python para resolver ecuaciones de Navier-Stokes
mediante el Método de Elementos Finitos.

CARACTERÍSTICAS:

Problemas:
  ✓ Estacionarios (steady-state)
  ✓ Evolutivos (time-dependent) ← NUEVO!

Espacios de elementos finitos:
  ✓ Taylor-Hood (P2-P1)
  ✓ MINI elements

Solucionadores estacionarios:
  ✓ Picard (robusto)
  ✓ Newton (convergencia cuadrática)

Solucionadores evolutivos:
  ✓ Semi-implícito (IMEX)
  ✓ Totalmente implícito (Picard/Newton)

Exportación:
  ✓ VTK para ParaView
  ✓ Animaciones temporales (.pvd)
  ✓ Gráficos de energía cinética

Documentación:
  ✓ Ejemplos completos
  ✓ Tests unitarios
  ✓ Guías de usuario

AUTOR: Mauricio Mendiluce
LICENCIA: MIT

Para documentación completa, ver:
  - README.md
  - UNSTEADY_IMPLEMENTATION.md
  - examples/
            """)
            input("\nPresiona Enter para continuar...")
        
        elif choice == "4":
            # Exit
            print("\n¡Hasta luego! 👋")
            sys.exit(0)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n❌ Programa interrumpido por el usuario")
        sys.exit(0)
    except Exception as e:
        print(f"\n\n❌ Error inesperado: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)