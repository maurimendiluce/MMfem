"""
Docstring for examples.non_regular_data
"""
#from netgen import gui
from ngsolve import CF, Draw, x, y, sqrt, Integrate,InnerProduct,IfPos,VTKOutput
from mmfem import rectangle_mesh, taylor_hood, mini_elements, picard_iteration,newton_iteration,stabilization_p1p1,newton_iteration_p1p1
import matplotlib.pyplot as plt
import numpy as np


def save_convergence_table(h, err, filename, title=""):
    """
    Guarda los datos de convergencia en un archivo txt
    
    Parameters:
    -----------
    h : array-like
        Tamaños de malla
    err : array-like
        Errores correspondientes
    filename : str
        Nombre del archivo (ej: 'taylor_hood.txt')
    title : str, optional
        Título descriptivo
    """
    h = np.array(h)
    err = np.array(err)
    
    with open(filename, 'w') as f:
        # Escribir título como comentario
        if title:
            f.write(f"# {title}\n")
        
        # Escribir encabezado
        f.write("# h_max           error_L4\n")
        
        # Escribir datos
        for i in range(len(h)):
            f.write(f"{h[i]:.6e}  {err[i]:.6e}\n")


def print_convergence_table(h, err, title):
    h = np.array(h)
    err = np.array(err)

    print(title)
    print("-" * 60)
    print(f"{'h_max':>12} | {'error L4':>12} | {'order':>8}")
    print("-" * 60)

    for i in range(len(err)):
        if i == 0:
            print(f"{h[i]:12.4e} | {err[i]:12.4e} | {'-':>8}")
        else:
            order = np.log(err[i] / err[i-1]) / np.log(h[i] / h[i-1])
            print(f"{h[i]:12.4e} | {err[i]:12.4e} | {order:8.5f}")

    print("-" * 60)

def run(h_max):

    tolerance = 1e-8
    viscosity = 1
    meshes = []
    solutions_Taylor2 = []
    solutions_Mini = []
    solutions_p1p1 = []

    for j in range(len(h_max)):
        mesh = rectangle_mesh(-1, 1, -1, 1, h=h_max[j])
        meshes.append(mesh)
        X_taylor2 = taylor_hood(mesh, "left|right|bottom|top")
        X_mini = mini_elements(mesh, "left|right|bottom|top")
        X_p1p1 = stabilization_p1p1(mesh, "left|right|bottom|top")
        #u_bc = CF((1,0))
        f = IfPos(x - (-1),IfPos(x - (-1+h_max[j]),IfPos(x - (1-h_max[j]), (1-x)/h_max[j], 1),(x+1)/h_max[j]), 0)
        u_bc = CF((f,0))
        print("="*70)
        print(f"Solving example with h = {h_max[j]}")
        print("="*70)
        print("Solving with Taylor Hood P2P1..")
        print(f"  Velocity DOFs: {X_taylor2.components[0].ndof}")
        solution_taylor2, info = newton_iteration(
            mesh=mesh,
            fespace=X_taylor2,
            dirichlet_boundaries="top",
            velocity_bc=u_bc,
            viscosity=viscosity,
            tolerance=tolerance,
            max_iterations=100,
            verbose=False
            )
        print(f"Converged: {info['converged']}")
        print(f"Iterations: {info['iterations']}")
        print(f"error: {info['final_error']}")
        solutions_Taylor2.append(solution_taylor2.components[0])

        print("Solving with Mini-Elements..")
        print(f"  Velocity DOFs: {X_mini.components[0].ndof}")
        solution_mini, info = newton_iteration(
            mesh=mesh,
            fespace=X_mini,
            dirichlet_boundaries="top",
            velocity_bc=u_bc,
            viscosity=viscosity,
            tolerance=tolerance,
            max_iterations=100,
            verbose=False
            )
        print(f"Converged: {info['converged']}")
        print(f"Iterations: {info['iterations']}")
        print(f"error: {info['final_error']}")
        solutions_Mini.append(solution_mini.components[0])

        print("Solving with p1p1..")
        solution_p1p1, info = newton_iteration_p1p1(
            mesh=mesh,
            fespace=X_p1p1,
            dirichlet_boundaries="top",
            velocity_bc=u_bc,
            viscosity=viscosity,
            tolerance=tolerance,
            max_iterations=100,
            verbose=False
            )
        print(f"Converged: {info['converged']}")
        print(f"Iterations: {info['iterations']}")
        print(f"error: {info['final_error']}")
        solutions_p1p1.append(solution_p1p1.components[0])

    return meshes,solutions_Mini,solutions_Taylor2,solutions_p1p1


def main():
    """Run simulation."""
    
    print("="*70)
    print("LID-DRIVEN CAVITY EXAMPLE REGULARIZATION")
    print("="*70)
    
    h_max = [1/4,1/8,1/16,1/32,1/64,1/128,1/256]
    errorL4_Taylor2 = []
    errorL4_Mini = []
    errorL4_p1p1 = []

    meshes,solutions_Mini,solutions_Taylor2,solutions_p1p1 = run(h_max= h_max)
    #meshes,solutions_Mini,solutions_Taylor2 = run(h_max= h_max)
    for j in range(len(meshes)-1):
        errorL4_Taylor2.append(sqrt(sqrt(Integrate(InnerProduct(solutions_Taylor2[j+1]-solutions_Taylor2[j],solutions_Taylor2[j+1]-solutions_Taylor2[j])**2,meshes[j+1]))))
        errorL4_Mini.append(sqrt(sqrt(Integrate(InnerProduct(solutions_Mini[j+1]-solutions_Mini[j],solutions_Mini[j+1]-solutions_Mini[j])**2,meshes[j+1]))))
        errorL4_p1p1.append(sqrt(sqrt(Integrate(InnerProduct(solutions_p1p1[j+1]-solutions_p1p1[j],solutions_p1p1[j+1]-solutions_p1p1[j])**2,meshes[j+1]))))
    

    h = np.array(h_max[:-1])
    #err_taylor2 = np.array(errorL4_Taylor2)
    #err_mini = np.array(errorL4_Mini)
    #err_p1p1 = np.array(errorL4_p1p1)
    #err_p1p1 = np.array([0.244399,0.171503,0.12124,0.0857048,0.060589])
    
    #print table
    print_convergence_table(h, errorL4_Taylor2, "Taylor-Hood")
    print_convergence_table(h, errorL4_Mini, "Mini-Elements")
    print_convergence_table(h, errorL4_p1p1,"P1P1 Stabilization")

    save_convergence_table(h, errorL4_Taylor2, "taylor-hood.txt", title="")
    save_convergence_table(h, errorL4_Mini, "mini-elements.txt", title="")
    save_convergence_table(h, errorL4_p1p1, "p1p1.txt", title="")
    
    #grafico de las g_h
    #plt.figure(2)
    #plt.grid(True, alpha=0.1)
    #plt.title(r'Regularization functions', fontsize=14)
    #X = np.linspace(-1,1,num=1000)
    #Y = np.ones_like(X) * 0.2
    #f1 = IfPos(x - (-1),IfPos(x - (-1+h_max[1]),IfPos(x - (1-h_max[1]), (1-x)/h_max[1], 1),(x+1)/h_max[1]), 0)
    #f4 = IfPos(x - (-1),IfPos(x - (-1+h_max[4]),IfPos(x - (1-h_max[4]), (1-x)/h_max[4], 1),(x+1)/h_max[4]), 0)
    #plt.plot(X, f1(meshes[1](X, Y)),label = rf"$g_h$ regularization with $h=1/8$")
    #plt.plot(X, f4(meshes[4](X, Y)),label = rf"$g_h$ regularization with $h=1/64$")
    #plt.legend()
    #plt.savefig('gh_plot.png', dpi=150)
    #print("  Functions regularization plots saved to 'gh_plot.png'")

    #vtk = VTKOutput(meshes[3],
    #            coefs=[solutions_Taylor2[3],solutions_Mini[3],solutions_p1p1[3]],
    #            names = ["TaylorHood","MiniElementes","P1P1"],
    #            filename="velocity_reg_example",
    #            subdivision=2)
    # Exporting the results:
    #vtk.Do()
    
    # Final summary
    #print("\n" + "="*70)
    #print("SIMULATION COMPLETE")
    #print("="*70)
    #print(f"✓ Converged in {info['iterations']} iterations")
    #print(f"✓ Final error: {info['final_error']:.2e}")
    #print(f"✓ Solution ready for visualization")
    #print("="*70)
    
    # Keep GUI open
    input("\nPress Enter to close visualization and exit...")


if __name__ == "__main__":
    main()