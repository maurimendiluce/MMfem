module MMfem

using DelaunayTriangulation, CairoMakie, LinearAlgebra, SparseArrays, LinearSolve

#Funciones para generar la malla
export rectangle_mesh
"""
    rectangle_mesh(a,b,c,d,N)

This function generates a triangulation in the rectangle [a,b]x[c,d] using N subdivisions in each interval.
"""
function rectangle_mesh(a,b,c,d,N)

    tri = triangulate_rectangle(a, b, c, d, N, N,delete_ghosts = true)
    return tri
end

export nodes_tri
"""
    nodes_tri(tri::Triangulation)

This function returns the array of nodes of the triangulation.
"""
function nodes_tri(tri)
    # función que devuelve la matriz de nodos de la triangulación

    points = get_points(tri)

    nodes = zeros(size(points)[1],2)
    for i=1:size(nodes)[1]
        nodes[i,1] = points[i][1]
        nodes[i,2] = points[i][2]
    end
    return nodes
end

export elements_tri
"""
    elements_tri(tri::Triangulation)

This function returns the numbering of each element of the triangulation.
"""    
function elements_tri(tri)
    #función que devuelve los triangulos numerados de la triangulación

    triangles = get_triangles(tri)
    elements = zeros(length(triangles),3)
    j=1
    for k in triangles
        elements[j,1],elements[j,2],elements[j,3] = k[1], k[2], k[3]
        j=j+1
    end
    
    elements = Int.(elements)

    return elements
end


export boundary_nodes
"""
    boundary_nodes(tri::Triangulation)

This function returns the numbering of the edge nodes.
"""
function boundary_nodes(tri)
    #función que devuelve la numeración de los nodos de borde

    get_boundary = get_boundary_nodes(tri)
    boundary = [get_boundary[1];get_boundary[2];get_boundary[3];get_boundary[4]]    
    
    return unique(boundary)
end

export mesh_data
"""
    mesh_data(tri::Triangulation)

"""
function mesh_data(tri)

    p = nodes_tri(tri)
    t = elements_tri(tri)
    b = boundary_nodes(tri)
    return p,t,b   
end

#Funciones para graficar
export plot_mesh
"""
    plot_mesh(tri::DelaunayTriangulation)

Plot the mesh
"""
function plot_mesh(tri)
    display(triplot(tri))
end

export scalar_plot
"""
    scalar_plot

This function graphs a scalar function.
"""
function scalar_plot(u,tri;title = "title")
    f, ax, tr  = tricontourf(tri, u,levels = 100,axis = (; aspect = 1, title))
    Colorbar(f[1, 2], tr)
    f
end

export vectorial_plot
"""
    vectorial_plot

This function graphs a vectorial function.
"""
function vectorial_plot(F,points)
    display(Plots.quiver(points[:,1],points[:, 2],quiver=(F[:, 1]/norm(F[:, 1]),F[:, 2]/norm(F[:, 2]))))
end


#funciones bases

export β
"""
    β(j,v,points,triangle)

Functions corresponding to the Lagrange basis
""" 
function β(j,v,points,triangle)
    x,y=v[1],v[2]
    vert=points[triangle,:]
    x1,x2,x3=vert[:,1]
    y1,y2,y3=vert[:,2]
    if j==1 
        area=(x3-x2)*(y1-y2)+(y2-y3)*(x1-x2)
        return ((x3-x2)*(y-y2)+(y2-y3)*(x-x2))/area
    elseif j==2 
        area=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1)
        return ((x3-x1)*(y-y1)+(y1-y3)*(x-x1))/area
    else 
        area=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1)
        return ((x2-x1)*(y-y1)+(y1-y2)*(x-x1))/area 
    end
end

export ∇β
"""
    ∇β(j,points,triangle)

"""
function ∇β(j,points,triangle)
    vert=points[triangle,:]
    x1,x2,x3=vert[:,1]
    y1,y2,y3=vert[:,2]
    if j==1 
        area=(x3-x2)*(y1-y2)+(y2-y3)*(x1-x2)
        return [(y2-y3)/area,(x3-x2)/area]
    elseif j==2
        area=(x3-x1)*(y2-y1)+(y1-y3)*(x2-x1)
        return [(y1-y3)/area,(x3-x1)/area]
    else 
        area=(x2-x1)*(y3-y1)+(y1-y2)*(x3-x1)
        return [(y1-y2)/area,(x2-x1)/area]
    end
end

#cuatraturas
export quad
"""
    quad(f,points,triangle;type="midpoints")

"""
function quad(f,points,triangle;type="midpoints")
    v1,v2,v3=points[triangle[1],:],points[triangle[2],:],points[triangle[3],:]
    area_T = abs(det([v2-v1 v3-v1])/2)
    if type == "midpoints"
        return (area_T/3)*(f((v1+v2)/2)+f((v1+v3)/2)+f((v2+v3)/2))
    end
    if type == "barycenter"
        return (area_T)*(f((v1+v2+v3)/3))
    end
end

#Matrices de rigidez
export Stiffnes_Δ
"""
    Stiffnes_Δ(tri::Triangulation)


"""
function Stiffnes_Δ(tri)
    p,t,b = mesh_data(tri)
    nn = size(p)[1]
    nt = size(t)[1]
    A = spzeros(nn,nn)
    for k=1:nt
        elem = t[k,:]
        v1,v2,v3=p[elem[1],:],p[elem[2],:],p[elem[3],:]
        area_k = abs(det([v2-v1 v3-v1])/2)
        M = zeros(2,3)
        M[:,1] = ∇β(1,p,elem)
        M[:,2] = ∇β(2,p,elem)
        M[:,3] = ∇β(3,p,elem)
        M = area_k*M'*M 
        A[elem,elem] = A[elem,elem] + M
    end
    return A
end

#Lado derecho
export fβ
"""
   fβ(f,j,v,points,triangle)

""" 
function fβ(f,j,v,points,triangle)
    return f(v)*β(j,v,points,triangle)
end

export LoadVector
"""
    LoadVector(f,tri)

"""
function LoadVector(f,tri)
    p,t,b = mesh_data(tri)
    nn = size(p)[1]
    nt = size(t)[1]
    F = spzeros(nn,1) 
    for k=1:nt 
        elem = t[k,:]
        local_vector = [quad(v->fβ(f,1,v,p,elem),p,elem),quad(v->fβ(f,2,v,p,elem),p,elem),quad(v->fβ(f,3,v,p,elem),p,elem)]
        F[elem] = F[elem]+local_vector
    end
    return F
end

#Boundary conditions

export dirichlet
"""
    dirichlet(tri,dirichlet_nodes,uD)

"""
function dirichlet(tri,dirichlet_nodes,uD)
    p,t,b = mesh_data(tri)
    u = spzeros(size(p)[1],1)
    u[dirichlet_nodes] = uD.(p[dirichlet_nodes],:)
    return u
end


#Solver 
export solver_problem
"""
    solver_problem(A,F)

"""
function solver_problem(A,F)
    A = Matrix(A)
    F = Matrix(F)
    return A\F
end


end # module MMfem
