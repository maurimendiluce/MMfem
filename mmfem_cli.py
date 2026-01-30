"""
MMfem Interactive Command Line Interface

Interactive tool for solving Navier-Stokes problems without writing code.
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
    
    # 1. Domain selection
    print_menu("1. Selecciona el dominio:", {
        "1": "Rectángulo (cavity flow)",
        "2": "Dominio en L (con singularidad)"
    })
    domain_type = get_choice("Opción", ["1", "2"])
    
    # 2. Mesh size
    print("\n2. Tamaño de malla inicial")
    if domain_type == "1":
        h = get_float("Tamaño de malla h", default=0.05, min_val=0.01)
    else:
        h = get_float("Tamaño de malla h", default=0.3, min_val=0.05)
    
    # 3. Viscosity
    print("\n3. Viscosidad")
    viscosity = get_float("Viscosidad (ν)", default=0.01, min_val=0.0001)
    reynolds = 1.0 / viscosity
    print(f"   → Número de Reynolds: Re ≈ {reynolds:.0f}")
    
    # 4. Solver selection
    print_menu("4. Selecciona el método de solución:", {
        "1": "Picard (robusto, convergencia lineal)",
        "2": "Newton (rápido, convergencia cuadrática)",
        "3": "Adaptativo (refinamiento automático de malla)"
    })
    solver_type = get_choice("Opción", ["1", "2", "3"])
    
    # 5. Solver parameters
    if solver_type in ["1", "2"]:
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
    else:
        # Adaptive parameters
        n_refinements = get_int("Número de ciclos de refinamiento", default=5, min_val=1)
        theta = get_float("Parámetro de marcado θ (0-1)", default=0.7, min_val=0.01)
        print_menu("Estrategia de marcado:", {
            "1": "Maximum (simple, efectiva)",
            "2": "Dörfler (teóricamente óptima)"
        })
        mark_choice = get_choice("Opción", ["1", "2"])
        marking_strategy = "maximum" if mark_choice == "1" else "dorfler"
        tolerance = 1e-8
        max_iter = None
        convection_form = None
    
    # Summary
    print("\n" + "="*80)
    print("RESUMEN DE CONFIGURACIÓN")
    print("="*80)
    print(f"  Dominio: {'Rectángulo' if domain_type == '1' else 'L-shaped'}")
    print(f"  Tamaño de malla: h = {h}")
    print(f"  Viscosidad: ν = {viscosity} (Re ≈ {reynolds:.0f})")
    
    if solver_type == "1":
        print(f"  Método: Picard")
        print(f"  Convección: {convection_form}")
        print(f"  Tolerancia: {tolerance:.2e}")
        print(f"  Max iteraciones: {max_iter}")
    elif solver_type == "2":
        print(f"  Método: Newton")
        print(f"  Tolerancia: {tolerance:.2e}")
        print(f"  Max iteraciones: {max_iter}")
    else:
        print(f"  Método: Adaptativo")
        print(f"  Ciclos de refinamiento: {n_refinements}")
        print(f"  Parámetro θ: {theta}")
        print(f"  Estrategia: {marking_strategy}")
    
    print("="*80)
    
    confirm = get_choice("¿Proceder con esta configuración? (s/n)", ["s", "n", "y"])
    
    if confirm == "n":
        return None
    
    # Build configuration
    config = {
        'domain_type': domain_type,
        'h': h,
        'viscosity': viscosity,
        'solver_type': solver_type,
        'tolerance': tolerance,
        'max_iter': max_iter,
        'convection_form': convection_form
    }
    
    if solver_type == "3":
        config['n_refinements'] = n_refinements
        config['theta'] = theta
        config['marking_strategy'] = marking_strategy
    
    return config


####

###


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
            
            #solution, info, mesh = solve_problem(config)
            
            #if solution is not None:
            #    visualize_results(solution, info, mesh, config)
            
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

CONCEPTOS BÁSICOS:

1. DOMINIO:
   - Rectángulo: Problema clásico de lid-driven cavity
   - L-shaped: Dominio con singularidad (para adaptatividad)

2. VISCOSIDAD (ν):
   - Controla el número de Reynolds: Re = 1/ν
   - Valores típicos: 0.01 (Re=100) a 0.001 (Re=1000)
   - Mayor viscosidad → flujo más lento y estable

3. MÉTODOS:
   - Picard: Robusto, bueno para Re bajo/medio
   - Newton: Rápido si hay buen punto inicial
   - Adaptativo: Refina malla automáticamente

4. PARÁMETROS ADAPTATIVOS:
   - θ (theta): 0.5-0.8, controla cuántos elementos refinar
   - Estrategia: Maximum (simple) o Dörfler (óptima)

RECOMENDACIONES:
- Comenzar con viscosidad alta (0.01-0.1)
- Usar malla gruesa inicialmente (h=0.1-0.2)
- Probar Picard primero, luego Newton si necesario
- Usar adaptativo para problemas con singularidades

Para más información: docs/BEST_PRACTICES.md
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
MMfem v0.1.0
Librería Python para resolver ecuaciones de Navier-Stokes
mediante el Método de Elementos Finitos.

CARACTERÍSTICAS:
✓ Espacios de elementos finitos: Taylor-Hood, MINI
✓ Solucionadores: Picard, Newton
✓ Refinamiento adaptativo de mallas
✓ Documentación completa
✓ Ejemplos y tests

AUTOR: Mauricio Mendiluce
LICENCIA: MIT
REPOSITORIO: https://github.com/maurimendiluce/MMfem

Para documentación completa, ver README.md
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