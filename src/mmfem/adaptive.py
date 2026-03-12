"""
Adaptive mesh refinement for Navier-Stokes problems.

This module provides a posteriori error estimation and adaptive mesh refinement
strategies for improving the accuracy of finite element solutions.

Functions:
    ACTUALIZAR
    estimate_error: Compute a posteriori error estimates
    mark_elements: Mark elements for refinement using Dörfler strategy
    adaptive_solve: Solve with adaptive mesh refinement
    compute_convergence_rates: Analyze convergence rates
    generate_convergence_plot: Visualize convergence behavior
    export_results_latex: Export results as LaTeX table
"""
import numpy as np
import matplotlib.pyplot as plt
from ngsolve import *

from .spaces import taylor_hood
from .solvers import picard_iteration

#FUNCIONES

def solve_estime_mark(mesh,dirichlet_labels, dirichlet_boundary, uD, viscosity = 1, theta = 0.7, domain = "convex", strategy = "maximum"):
    
    h = specialcf.mesh_size
    n = specialcf.normal(mesh.dim)
    X = taylor_hood(mesh,dirichlet_labels)
    gf,info = picard_iteration(mesh,X,dirichlet_boundary,uD,viscosity,verbose=False)

    uh = gf.components[0]
    ph = gf.components[1]

    V1 = H1(mesh,order=2,autoupdate=True)
    uh_0 = GridFunction(V1)
    uh_1 = GridFunction(V1)
    uh_0.Set(uh[0])
    uh_1.Set(uh[1])
    ph_x = grad(ph)[0]
    ph_y = grad(ph)[1]
    laplacian_0 = Trace(uh_0.Operator("hesse"))
    laplacian_1 = Trace(uh_1.Operator("hesse"))
    convect_0 = uh_0*grad(uh_0)[0]+uh_1*grad(uh_0)[1]
    convect_1 = uh_0*grad(uh_1)[0]+uh_1*grad(uh_1)[1]
    residual = ((laplacian_0-ph_x-convect_0)**4+(laplacian_1-ph_y-convect_1)**4)*dx

    uh_0_x = grad(uh_0)[0]
    uh_1_y = grad(uh_1)[1]

    div = (uh_0_x+uh_1_y)**4*dx
    jump = ( (grad(uh_0)-grad(uh_0).Other())*n )**4 *dx(element_vb=BND)+(h**5)*( (grad(uh_1)-grad(uh_1).Other())*n )**4 *dx(element_vb=BND)

    if domain == "convex":
        alpha, beta, gamma = 8, 4, 5
    elif domain == "non-convex":
        alpha, beta, gamma = 20/3, 0, 11/3
    else:
        raise ValueError(f"Invalid domain type")

    eta_residual = (h**alpha) * residual

    if beta > 0:
        eta_div = (h**beta) * div
    else:
        eta_div = div

    eta_jump = (h**gamma) * jump

    eta_local = Integrate(eta_residual+eta_div+eta_jump,mesh,element_wise=True)
    eta_global = sqrt(sqrt(sum(eta_local)))
    nv = mesh.nv #num de vertices
    ne = mesh.ne #num de triangulos

    if strategy == "maximum":
        eta_max = max(eta_local)
        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, eta_local[el.nr] > theta*eta_max)

    elif strategy == "average":
        eta_avg = (1/ne)*sum(eta_local)
        for el in mesh.Elements(): 
            mesh.SetRefinementFlag(el, eta_local[el.nr] > theta*eta_avg)

    return gf,nv,ne,eta_local,eta_global
    
def adaptive_method(mesh,dirichlet_labels,dirichlet_boundary,uD,n_ref=5,theta=0.7,viscosity=1,domain = "convex",strategy="maximum",verbose = True):
    ngsglobals.msg_level = 0
    if verbose:
        print("="*80)
        print("ADAPTIVE MESH REFINEMENT FOR NAVIER-STOKES")
        print("="*80)
        print(f"Initial mesh: {mesh.ne} elements, {mesh.nv} vertices")
        print(f"Viscosity: {viscosity}")
        print(f"Max refinement cycles: {n_ref}")
        print("-"*80)

    history = {
        'n_vertices': [],
        'n_elements': [],
        'n_dofs': [],
        'eta_global': [],
        'error_L4': [],
    }
    
    ref = 0
    while ref<n_ref:
        if verbose:
            print(f"\nRef {ref+1}/{n_ref}")
            print("-"*80)
        
        gf,nv,ne,eta_local,eta_global = solve_estime_mark(mesh,dirichlet_labels,dirichlet_boundary,uD,viscosity,theta,domain,strategy)
        uh_old = gf.components[0]
        mesh.Refine()
        X = taylor_hood(mesh,dirichlet_labels)
        solution = picard_iteration(mesh,X,dirichlet_boundary,uD,viscosity,verbose=False)
        uh,ph = solution.components[0],solution.components[1]
        error_L4 = sqrt(sqrt(Integrate(InnerProduct(uh-uh_old,uh-uh_old)**2,mesh)))
        history["nv"].append(nv)
        history["ne"].append(ne)
        history["error_L4"].append(error_L4)
        history["eta_global"].append(eta_global)
        ref += 1
    
    Draw(uh,mesh,"velocidad")
    Draw(ph,mesh,"presion")
    
    return history