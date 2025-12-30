from triangulation import RectangleMesh
from spaces import FEMSpaces
from problems import Stokes
from ngsolve.webgui import Draw


def main():
    rect = RectangleMesh(0, 1, 0, 1, 0.1)
    mesh = rect.generate()
    etiquetas = rect.etiquetas()

    spaces = FEMSpaces(mesh, dirichlet=rect.etiquetas())
    X = spaces.build(interactive=True)

    problem = Stokes(mesh,X)
    gf = problem.solve_stokes(etiq_dirich=etiquetas,uD=(0,0))
    Draw(gf.components[0],mesh,"Velocity")
    
    input("presione Enter para terminar...")


if __name__ == "__main__":
    main()
