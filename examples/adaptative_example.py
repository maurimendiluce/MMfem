"""
Ejemplo de Refinamiento Adaptativo de Mallas

Este ejemplo demuestra el uso de refinamiento adaptativo para el dominio en L,
que tiene una singularidad en la esquina reentrante.
"""
from netgen import gui
from ngsolve import CF, Draw
from mmfem import (
    rectangle_mesh,
    adaptive_method,
    boundary_string,
)


def main():
    """Ejemplo de refinamiento adaptativo."""
    
    print("="*80)
    print(" "*15 + "REFINAMIENTO ADAPTATIVO")
    print("="*80)
    
    # Crear malla inicial
    print("Creando malla inicial...")
    h_initial = 0.3
    mesh = rectangle_mesh(0,1,0,1,h=h_initial)
    dirichlet_labels = boundary_string(mesh)
    print(f"   Elementos iniciales: {mesh.ne}")
    print(f"   Vértices iniciales: {mesh.nv}")
    
    # Condiciones de borde
    u_bc = CF((1, 0))  # Velocidad horizontal en la entrada
    
    viscosity_default = 1

    viscosity_input = input("Viscosity (default nu=1) = ")

    if viscosity_input.strip() == "":
        viscosity = viscosity_default
    else:
        viscosity = float(viscosity_input)

    
    print(f"   Viscosidad: {viscosity}")
    print(f"   Reynolds: {1/viscosity:.0f}")

    n_ref_default = 5
    n_ref_input = input("ref cycles (default: 5)) = ")
    if n_ref_input == "":
        n_ref = n_ref_default
    else:
        n_ref = float(n_ref_input)

    theta_default = 0.7
    theta_input = input("theta (default: 0.7)) = ")
    if theta_input == "":
        theta = theta_default
    else:
        theta = float(theta_input)

    domains_types = {'convex': "convex",'non-convex': "non-convex"}
    
    print("Elegir TIPO DE DOMINIO")
    print("─" * 70)
    print(f"{domains_types['convex']}: 0 | {domains_types['non-convex']}: 1 ")
    domain_option = int(input("Option: "))
    if domain_option == 0:
        domain = domains_types["convex"]
    elif domain_option == 1:
        domain = domains_types["non-convex"]

    strategy_options = {'maximum': "maximum",'average': "average"}
    
    print("Elegir ESTRATEGIA DE MARCADO")
    print("─" * 70)
    print(f"{strategy_options['maximum']}: 0 | {strategy_options['average']}: 1 ")
    strategy_option = int(input("Option: "))
    if strategy_option == 0:
        strategy = strategy_options["maximum"]
    elif strategy_option == 1:
        strategy = strategy_options["average"]
    
    # Resolver con refinamiento adaptativo
    print("\n3. Resolviendo con refinamiento adaptativo...")
    print("-"*80)
    
    history = adaptive_method(mesh,dirichlet_labels,"top",u_bc,n_ref,theta,viscosity,domain,strategy)
    #solution, history = adaptive_solve(
    #    mesh=mesh,
    #    dirichlet_boundaries="top",
    #    dirichlet_labels = dirichlet_labels,
    #    velocity_bc=u_bc,
    #    viscosity=viscosity,
    #    n_refinements=5,
    #    theta=0.7,
    #    marking_strategy="maximum",
    #    domain_type="convex",
    #    verbose=True
    #)
    
    # Resultados
    print("\n" + "="*80)
    print("RESULTADOS")
    print("="*80)
    
    #print(f"\nElementos:")
    #print(f"  Inicial: {history['n_elements'][0]}")
    #print(f"  Final:   {history['n_elements'][-1]}")
    #print(f"  Factor:  {history['n_elements'][-1]/history['n_elements'][0]:.1f}x")
    
    #print(f"\nError:")
    #print(f"  Inicial: {history['eta_global'][0]:.6e}")
    #print(f"  Final:   {history['eta_global'][-1]:.6e}")
    #print(f"  Reducción: {history['eta_global'][0]/history['eta_global'][-1]:.1f}x")
    
    #if history['convergence_rate_eta']:
    #    avg_rate = sum(history['convergence_rate_eta']) / len(history['convergence_rate_eta'])
    #    print(f"\nTasa de convergencia promedio: {avg_rate:.2f}")
    #    print(f"  (Teórico óptimo adaptativo: ~1.0)")
    #    print(f"  (Teórico uniforme: ~0.5)")
    
    # Generar gráficos
    #print("\n" + "="*80)
    #print("GENERANDO VISUALIZACIÓN")
    #print("="*80)
    
    #try:
    #    generate_convergence_plot(
    #        history,
    #        filename='adaptive_convergence.png',
    #        show=False
    #    )
    #    print("✓ Gráfico guardado: adaptive_convergence.png")
    #except Exception as e:
    #    print(f"⚠ Error generando gráfico: {e}")
    
    #try:
    #    export_results_latex(history, filename='adaptive_results.tex')
    #    print("✓ Tabla LaTeX guardada: adaptive_results.tex")
    #except Exception as e:
    #    print(f"⚠ Error generando tabla: {e}")
    
    # Visualizar solución
    #print("\n✓ Abriendo ventanas de visualización...")
    #velocity = solution.components[0]
    #pressure = solution.components[1]
    
    #Draw(velocity, mesh, "Velocidad")
    #Draw(pressure, mesh, "Presión")
    
    #print("\n" + "="*80)
    #print("COMPLETO")
    #print("="*80)
    #print("""
#Archivos generados:
#  - adaptive_convergence.png (gráfico de convergencia)
#  - adaptive_results.tex (tabla LaTeX)

#Ventanas de visualización abiertas.
#    """)
#    
    input("Presiona Enter para cerrar...")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n❌ Programa interrumpido")
    except Exception as e:
        print(f"\n\n❌ Error: {e}")
        import traceback
        traceback.print_exc()