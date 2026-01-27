from ngsolve import *
from FEMSpaces import *
from solvers import *
import matplotlib.pyplot as plt
import numpy as np

def SolveEstimeMark(mesh,etiquetas_dirichlet,etiqueta_uD,uD,theta=0.7, convex = 0):
    h = specialcf.mesh_size
    n = specialcf.normal(mesh.dim)
    X = taylor_hood(mesh,etiquetas_dirichlet)
    gf = picard_iter_adapt(mesh,X,etiqueta_uD,uD)

    uh = gf.components[0]
    ph = gf.components[1]

    if convex == 0:
        #estimamos
        ph_x = grad(ph)[0]
        ph_y = grad(ph)[1]
        V1 = H1(mesh,order=2,autoupdate=True)
        uh_0 = GridFunction(V1)
        uh_1 = GridFunction(V1)
        uh_0.Set(uh[0])
        uh_1.Set(uh[1])
        laplaciano_0 = Trace(uh_0.Operator("hesse"))
        laplaciano_1 = Trace(uh_1.Operator("hesse"))
        convect_0 = uh_0*grad(uh_0)[0]+uh_1*grad(uh_0)[1]
        convect_1 = uh_0*grad(uh_1)[0]+uh_1*grad(uh_1)[1]
        estima1 = (h**8)*((laplaciano_0-ph_x-convect_0)**4+(laplaciano_1-ph_y-convect_1)**4)*dx

        uh_0_x = grad(uh_0)[0]
        uh_1_y = grad(uh_1)[1]
        estima2 = (h**4)*(uh_0_x+uh_1_y)**4*dx
        
        jump = (h**5)*( (grad(uh_0)-grad(uh_0).Other())*n )**4 *dx(element_vb=BND)+(h**5)*( (grad(uh_1)-grad(uh_1).Other())*n )**4 *dx(element_vb=BND)
        
        eta = Integrate(estima1+estima2+jump,mesh,element_wise=True)
    
    if convex == 1:
        #estimamos
        ph_x = grad(ph)[0]
        ph_y = grad(ph)[1]
        V1 = H1(mesh,order=2,autoupdate=True)
        uh_0 = GridFunction(V1)
        uh_1 = GridFunction(V1)
        uh_0.Set(uh[0])
        uh_1.Set(uh[1])
        laplaciano_0 = Trace(uh_0.Operator("hesse"))
        laplaciano_1 = Trace(uh_1.Operator("hesse"))
        convect_0 = uh_0*grad(uh_0)[0]+uh_1*grad(uh_0)[1]
        convect_1 = uh_0*grad(uh_1)[0]+uh_1*grad(uh_1)[1]
        estima1 = (h**(20/3))*((laplaciano_0-ph_x-convect_0)**4+(laplaciano_1-ph_y-convect_1)**4)*dx

        uh_0_x = grad(uh_0)[0]
        uh_1_y = grad(uh_1)[1]
        estima2 = (uh_0_x+uh_1_y)**4*dx
        
        jump = (h**(11/3))*( (grad(uh_0)-grad(uh_0).Other())*n )**4 *dx(element_vb=BND)+(h**2)*( (grad(uh_1)-grad(uh_1).Other())*n )**4 *dx(element_vb=BND)
        
        eta = Integrate(estima1+estima2+jump,mesh,element_wise=True)

    eta_global = sqrt(sqrt(sum(eta)))
    #eta_global = sum(sqrt(sqrt(eta)))
    nv = mesh.nv #num de vertices
    ne = mesh.ne #num de triangulos

    #marcamos
    eta_max = max(eta)
    for el in mesh.Elements():
        mesh.SetRefinementFlag(el, eta[el.nr] > theta*eta_max)

    return gf,nv,ne,eta,eta_global

def adaptative_method(cant_ref,mesh,etiquetas_dirichlet,etiqueta_uD,uD,theta=0.7):

    datos = {"nv":[],"ne":[],"error":[],"eta":[]}
    ref = 0
    while ref<cant_ref:
        gf,nv,ne,eta,eta_global = SolveEstimeMark(mesh,etiquetas_dirichlet,etiqueta_uD,uD)
        uh_old = gf.components[0]
        mesh.Refine()
        X = taylor_hood(mesh,etiquetas_dirichlet)
        solution = picard_iter_adapt(mesh,X,etiqueta_uD,uD)
        uh,ph = solution.components[0],solution.components[1]
        eL4 = sqrt(sqrt(Integrate(InnerProduct(uh-uh_old,uh-uh_old)**2,mesh)))
        datos["nv"].append(nv)
        datos["ne"].append(ne)
        datos["error"].append(eL4)
        datos["eta"].append(eta_global)
        ref += 1
    
    Draw(uh,mesh,"velocidad")
    Draw(ph,mesh,"presion")
    return datos

def order(datos):
    x = datos["ne"]
    eta = datos["eta"]
    error = datos["error"]
    plt.xlabel(r"$log(elements)$")
    plt.ylabel(r"$log(\eta)$")
    plt.plot(log(x),log(eta), "+")
    
    m,b=np.polyfit(log(x),log(eta),1)
    m1,b1=np.polyfit(log(x),log(error),1)
    x_val = np.linspace(np.min(log(x))-0.1,np.max(log(x))+0.1,100)
    y_val = m*x_val+b
    y_error = m1*x_val+b1
    plt.plot(x_val,y_val,linestyle = 'dashed',label=f"orden en eta: {np.round(m,2)}")
    plt.plot(x_val,y_error,label=f"orden del error L4: {np.round(m1,2)}")
    plt.legend()
    plt.show()

def results(datos):

    x = np.array(datos["ne"])
    eta = np.array(datos["eta"])
    error = datos["error"]

    orden_eta = []
    orden_error = []
    for i in range(1, len(x)):
        m, b = np.polyfit(np.log(x[:i+1]), np.log(eta[:i+1]), 1)
        m1, b1 = np.polyfit(np.log(x[:i+1]), np.log(error[:i+1]), 1)
        orden_eta.append(m)
        orden_error.append(m1)
    
    data = np.column_stack((x[1:],error[1:],orden_error, eta[1:], orden_eta))
    
    # --------- archivo LaTeX ----------
    with open("results.txt", "w") as f:
        f.write(r"\begin{tabular}{rcccc}" + "\n")
        f.write(r"\hline" + "\n")
        f.write(r"$n_e$ & error & ord(error) & $\eta$ & ord($\eta$) \\" + "\n")
        f.write(r"\hline" + "\n")

        for xi, ei, oei, etai, oi in data:
            f.write(
                f"{int(xi)} & "
                f"{ei:.3e} & "
                f"{oei:.2f} & "
                f"{etai:.3e} & "
                f"{oi:.2f} \\\\\n"
            )

        f.write(r"\hline" + "\n")
        f.write(r"\end{tabular}" + "\n")

    #np.savetxt(
    #    "results.txt",
    #    data,
    #    header="x    error    orden_error    eta    orden_eta",
    #    fmt="%.6f",
    #)

    print("\nResultados:")
    print(f"{'ne':>10} {'error':>15} {'orden_error':>15} {'eta':>15} {'orden_eta':>15}")
    print("-" * 60)

    for xi, ei, oei, etai, oi in data:
        print(f"{xi:10.0f} {ei:15.6e} {oei:15.4f} {etai:15.6e} {oi:15.4f}")