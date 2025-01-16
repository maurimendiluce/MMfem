using DelaunayTriangulation, CairoMakie, LinearAlgebra, SparseArrays, Plots

#Funciones para generar la malla
"""
    rectangle_mesh(a,b,c,d,N)

This function generates a triangulation in the rectangle [a,b]x[c,d] using N subdivisions in each interval.
"""
function rectangle_mesh(a,b,c,d,N)

    tri = triangulate_rectangle(a, b, c, d, N, N,delete_ghosts = true)
    return tri
end

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

"""
    plot_mesh(tri::DelaunayTriangulation)

Plot the mesh
"""
function plot_mesh(tri)
    display(triplot(tri))
end


"""
    scalar_plot

This function graphs a scalar function.
"""
function scalar_plot(u,tri;title = "title")
    f, ax, tr  = tricontourf(tri, u,levels = 100,axis = (; aspect = 1, title))
    Colorbar(f[1, 2], tr)
    f
end


"""
    vectorial_plot

This function graphs a vectorial function.
"""
function vectorial_plot(F,points)
    display(Plots.quiver(points[:,1],points[:, 2],quiver=(F[:, 1]/norm(F[:, 1]),F[:, 2]/norm(F[:, 2]))))
end


#funciones bases


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

function int_div_x(i,j,p,elem) #beta*∇ₓbeta
    v₁ = p[elem[1],:]
    v₂ = p[elem[2],:]
    v₃ = p[elem[3],:]
    B = [v₂-v₁ v₃-v₁]
    βₓ = ∇β(j,p,elem)[1]

    a₁ = β(i,0.5*(v₁+v₂),p,elem)*βₓ
    a₂ = β(i,0.5*(v₁+v₃),p,elem)*βₓ
    a₃ = β(i,0.5*(v₂+v₃),p,elem)*βₓ

    return (a₁+a₂+a₃)*abs(det(B))/6
end

function int_div_y(i,j,p,elem) #beta*∇_y beta
    v₁ = p[elem[1],:]
    v₂ = p[elem[2],:]
    v₃ = p[elem[3],:]
    B = [v₂-v₁ v₃-v₁]
    β_y = ∇β(j,p,elem)[2]

    a₁ = β(i,0.5*(v₁+v₂),p,elem)*β_y
    a₂ = β(i,0.5*(v₁+v₃),p,elem)*β_y
    a₃ = β(i,0.5*(v₂+v₃),p,elem)*β_y

    return (a₁+a₂+a₃)*abs(det(B))/6
end

function Stiffnes_div(tri)
    p,t,b = mesh_data(tri)
    n_nodes = size(p,1)
    n_elem = size(t,1) 

    Bx = spzeros(n_nodes,n_nodes)
    By = spzeros(n_nodes,n_nodes)

    for k=1:n_elem
        elem = t[k,:]
        Mx = zeros(3,3)
        My = zeros(3,3)
        for i=1:3
            for j=1:3
                Mx[i,j] = int_div_x(i,j,p,elem)
                My[i,j] = int_div_y(i,j,p,elem)
            end
        end
        Bx[elem,elem] = Bx[elem,elem] + Mx
        By[elem,elem] = Bx[elem,elem] + My 
    end
    
    return Bx,By
end

function int_G(i,j,p,elem)
    v₁ = p[elem[1],:]
    v₂ = p[elem[2],:]
    v₃ = p[elem[3],:]
    B = [v₂-v₁ v₃-v₁]
    Π = 1/3
    a₁ = (β(i,0.5*(v₁+v₂),p,elem)-Π)*(β(j,0.5*(v₁+v₂),p,elem)-Π)
    a₂ = (β(i,0.5*(v₁+v₃),p,elem)-Π)*(β(j,0.5*(v₁+v₃),p,elem)-Π)
    a₃ = (β(i,0.5*(v₂+v₃),p,elem)-Π)*(β(j,0.5*(v₂+v₃),p,elem)-Π)

    return (a₁+a₂+a₃)*abs(det(B))/6
end

function Stiffnes_G(tri) #matriz de estabilizacion
    p,t,b = mesh_data(tri)
    n_nodes = size(p,1)
    n_elem = size(t,1) 
    G = spzeros(n_nodes,n_nodes)

    for k=1:n_elem
        elem = t[k,:]
        M = zeros(3,3)
        for i=1:3
            for j=1:3
                M[i,j] = int_G(i,j,p,elem)
            end
        end
        G[elem,elem] = G[elem,elem] + M
    end
    return G
end

#Lado derecho
"""
   fβ(f,j,v,points,triangle)

""" 
function fβ(f,j,v,points,triangle)
    return f(v)*β(j,v,points,triangle)
end

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
"""
    dirichlet(tri,dirichlet_nodes,uD)

"""
function dirichlet(tri,dirichlet_nodes,uD)
    p,t,b = mesh_data(tri)
    u = spzeros(size(p)[1],1)
    u[dirichlet_nodes] = uD.(p[dirichlet_nodes],:)
    return u
end


#Calculo de Gradiente 
function search_node(nodo,elements)
    n_elem = size(elements)[1]
    for k = 1:n_elem
        elem = elements[k,:]
        if nodo == elem[1]
            return [k,1]
        elseif nodo == elem[2]
            return [k,2]
        elseif nodo == elem[3]
            return [k,3]
        end
    end
end

function grad_calculate(uₕ,tri)
    
    p,t,b = mesh_data(tri)
    n_elem = size(t)[1]
    n_nodes = size(p)[1]

    grad = zeros(n_nodes,2)

    for i=1:n_nodes
        el,lugar = search_node(i,t)
        a₁ = i
        if lugar == 1
            a₂ = t[el,2]
            a₃ = t[el,3]
            r₁ = p[a₂,:]-p[a₁,:]
            r₂ = p[a₃,:]-p[a₁,:]
        elseif lugar == 2
            a₂ = t[el,1]
            a₃ = t[el,3]
            r₁ = p[a₂,:]-p[a₁,:]
            r₂ = p[a₃,:]-p[a₁,:]
        elseif lugar == 3
            a₂ = t[el,1]
            a₃ = t[el,2]
            r₁ = p[a₂,:]-p[a₁,:]
            r₂ = p[a₃,:]-p[a₁,:]
        end
        A = [r₁./norm(r₁) r₂./norm(r₂)]
        A = A'
        b = [(uₕ[a₂]-uₕ[a₁])/norm(r₁);(uₕ[a₃]-uₕ[a₁])/norm(r₂)]
        parcial = A\b
        grad[i,:] = parcial
    end
    return grad
end

