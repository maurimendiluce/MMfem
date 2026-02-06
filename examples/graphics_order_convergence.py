import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, NullFormatter


def read_convergence_data(filename):
    """
    Lee los datos de convergencia desde un archivo txt
    
    Parameters:
    -----------
    filename : str
        Nombre del archivo a leer
        
    Returns:
    --------
    h : numpy array
        Tamaños de malla
    err : numpy array
        Errores
    """
    h = []
    err = []
    
    with open(filename, 'r') as f:
        for line in f:
            # Saltar líneas de comentarios
            if line.startswith('#'):
                continue
            
            # Leer datos
            values = line.split()
            if len(values) == 2:
                h.append(float(values[0]))
                err.append(float(values[1]))
    
    return np.array(h), np.array(err)

def make_graphics(h,err1,err2,err3):

    #grafico de convergencia
    # recta de referencia O(h^{1/2}) separada
    C = 0.6 * err1.max() / (h.min()**0.5)   # constante elegida a mano
    ref = C * h**0.5

    plt.figure(1)
    plt.loglog(h, err1, marker='o', linestyle='-.',linewidth=0.5, markersize=6, label=r'Taylor-Hood $P_2P_1$')
    plt.loglog(h, err2, marker='*', linestyle='-.',linewidth=0.5, markersize=6, label=r'Mini-Elements')
    plt.loglog(h, err3, marker='^', linestyle='-.',linewidth=0.5, markersize=6, label=r'$P_1P_1$ Stabilization')
    plt.loglog(h, ref, '--', linewidth=2, label=r'$O(h^{1/2})$')

        
    plt.xlabel(r'$log(h)$', fontsize=12)
    plt.ylabel(r'$log(e_{L^4}$)', fontsize=12)
    #plt.title('Order visualization', fontsize=14)
    plt.grid(True,which='minor', linestyle=':', alpha=0.4)
    plt.legend(fontsize=9, framealpha=0.95, loc='best')
    #plt.show()

    plt.savefig('convergence_order_01.png', dpi=300, bbox_inches='tight')
    print("  Convergence plot saved to 'convergence_order.png'")

def main():

    h,err_taylor2 = read_convergence_data("taylor-hood_01.txt")
    h,err_mini = read_convergence_data("mini-elements_01.txt")
    h,err_p1p1 = read_convergence_data("p1p1_01.txt")

    make_graphics(h,err_taylor2,err_mini,err_p1p1)

if __name__ == "__main__":
    main()
