from ngsolve import *

class FEMSpaces:
    """
    Construcción de espacios FEM para problemas incomprensibles.
    """

    def __init__(self, mesh, dirichlet):
        """
        Parameters
        ----------
        mesh : ngsolve.Mesh
        dirichlet : str
            Etiquetas de borde para condición de Dirichlet.
        """
        self.mesh = mesh
        self.dirichlet = dirichlet

    def build(self, method="taylor-hood", interactive=False):
        if interactive:
            method = self._choose_method()

        if method == "taylor-hood":
            V = VectorH1(self.mesh, order=2,dirichlet=self.dirichlet)
            Q = H1(self.mesh, order=1)

        elif method == "mini":
            V = VectorH1(self.mesh, order=1,dirichlet=self.dirichlet)
            V.SetOrder(TRIG, 3)
            V.Update()
            Q = H1(self.mesh, order=1)

        else:
            raise ValueError("Método FEM desconocido.")

        N = NumberSpace(self.mesh)
        return FESpace([V, Q, N])

    def _choose_method(self):
        opciones = {
            "1": ("taylor-hood", "Taylor–Hood P2–P1"),
            "2": ("mini", "Mini-elements"),
        }

        print("Opciones disponibles:")
        for k, (_, name) in opciones.items():
            print(f"{k}) {name}")

        opcion = input("Elegí el método a usar: ")

        if opcion not in opciones:
            raise ValueError("Opción inválida.")

        method, name = opciones[opcion]
        print(f"Elegiste: {name}")
        return method