#Packege

using DelaunayTriangulation, CairoMakie, LinearAlgebra, SparseArrays, Plots

#### Mesh functions ####

function rectangle_mesh(a,b,c,d,N)

    tri = triangulate_rectangle(a, b, c, d, N, N,delete_ghosts = true)
    return tri
end

function nodes_tri(tri)
    
    points = get_points(tri)

    nodes = zeros(size(points)[1],2)
    for i=1:size(nodes)[1]
        nodes[i,1] = points[i][1]
        nodes[i,2] = points[i][2]
    end
    return nodes
end

function elements_tri(tri)

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

function boundary_nodes(tri)

    get_boundary = get_boundary_nodes(tri)
    boundary = [get_boundary[1];get_boundary[2];get_boundary[3];get_boundary[4]]    
    
    return unique(boundary)
end


struct Mesh 
    nodes::Matrix
    elements::Matrix
    boundary::Vector
    nv::Int
    nt::Int
end

function create_mesh(a,b,c,d,N;graph=true)

    tri = rectangle_mesh(a,b,c,d,N)
    nodes = nodes_tri(tri)
    elements = elements_tri(tri)
    boundary = boundary_nodes(tri)
    nv = size(nodes_tri(tri))[1]
    nt = size(elements_tri(tri))[1]
    mesh = Mesh(nodes,elements,boundary,nv,nt)

    if graph == true
        display(triplot(tri))
    end

    display(statistics(tri))

    return mesh,tri
end

### Base functions ####

function ϕ(j,v,elem,mesh)
    x = v[1]
    y = v[2]
    vert = mesh.nodes[elem,:]
    x1,x2,x3 = vert[:,1]
    y1,y2,y3 = vert[:,2]
    if j==1 
        δ = (x3-x2)*(y1-y2)+(y2-y3)*(x1-x2)
        return ((x3-x2)*(y-y2)+(y2-y3)*(x-x2))/δ
    elseif j==2 
        δ = (x3-x1)*(y2-y1)+(y1-y3)*(x2-x1)
        return ((x3-x1)*(y-y1)+(y1-y3)*(x-x1))/δ
    else 
        δ = (x2-x1)*(y3-y1)+(y1-y2)*(x3-x1)
        return ((x2-x1)*(y-y1)+(y1-y2)*(x-x1))/δ
    end
end

function ∇ϕ(j,elem,mesh)
    vert = mesh.nodes[elem,:]
    x1,x2,x3 = vert[:,1]
    y1,y2,y3 = vert[:,2]
    if j==1 
        δ = (x3-x2)*(y1-y2)+(y2-y3)*(x1-x2)
        return [(y2-y3)/δ , (x3-x2)/δ]
    elseif j==2 
        δ = (x3-x1)*(y2-y1)+(y1-y3)*(x2-x1)
        return [(y1-y3)/δ , (x3-x1)/δ]
    else 
        δ = (x2-x1)*(y3-y1)+(y1-y2)*(x3-x1)
        return [(y1-y2)/δ , (x2-x1)/δ]
    end
end

function midpoints(elem,mesh)
    vert = mesh.nodes[elem,:]
    x1,x2,x3 = vert[:,1]
    y1,y2,y3 = vert[:,2]
    return 0.5*[[x1+x2,y1+y2],[x1+x3,y1+y3],[x2+x3,y2+y3]]
end

function quad(elem,mesh)
    vert = mesh.nodes[elem,:]
    x1,x2,x3 = vert[:,1]
    y1,y2,y3 = vert[:,2]
    return (x2-x1)*(y3-y1)+(y1-y2)*(x3-x1)
end

function ∫ϕ(j,elem,mesh)
    p_med=midpoints(elem,mesh)
    a1=ϕ(j, p_med[1],elem,mesh)
    a2=ϕ(j, p_med[2],elem,mesh)
    a3=ϕ(j, p_med[3],elem,mesh)
    
    return abs(quad(elem,mesh)/6)*(a1+a2+a3)
end
    
### Matrixs ###

function stiffness_Δ(mesh)
    
    A = zeros(mesh.nv,mesh.nv)
    
    for t=1:mesh.nt
        elem = mesh.elements[t,:]
    
        S = zeros(2,3)
        S[:,1]=∇ϕ(1,elem,mesh)
        S[:,2]=∇ϕ(2,elem,mesh)
        S[:,2]=∇ϕ(3,elem,mesh)
        M=(abs(quad(elem,mesh))/2)*S'*S
        A[elem,elem] = A[elem,elem] + M 
    end
    return A
end

function ∫Bx(i,j,elem,mesh)
    p_med=midpoints(elem,mesh)

    ϕx=∇ϕ(j,elem,mesh)[1]
    a1=ϕ(i, p_med[1],elem,mesh)*ϕx
    a2=ϕ(i, p_med[2],elem,mesh)*ϕx
    a3=ϕ(i, p_med[3],elem,mesh)*ϕx
    
    return (a1+a2+a3)*abs(quad(elem,mesh))/6
end

function ∫By(i,j,elem,mesh)
    p_med=midpoints(elem,mesh)

    ϕy=∇ϕ(j,elem,mesh)[2]
    a1=ϕ(i, p_med[1],elem,mesh)*ϕy
    a2=ϕ(i, p_med[2],elem,mesh)*ϕy
    a3=ϕ(i, p_med[3],elem,mesh)*ϕy
    
    return (a1+a2+a3)*abs(quad(elem,mesh))/6
end

function stiffness_div(mesh)

    Bx = zeros(mesh.nv,mesh.nv)
    By = zeros(mesh.nv,mesh.nv)

    for t=1:mesh.nt
        elem=mesh.elements[t,:]
        Sx = zeros(3,3)
        Sy = zeros(3,3)
        for i=1:3
            for j=1:3
                Sx[i,j]=∫Bx(i,j,elem,mesh)
                Sy[i,j]=∫By(i,j,elem,mesh)
            end
        end
        Bx[elem,elem] = Bx[elem,elem] + Sx
        By[elem,elem] = By[elem,elem] + Sy
    end
                
    return Bx,By
end

function ∫convect(w₁,w₂,i,j,elem,mesh)

    p_med = midpoints(elem,mesh)

    w₁_12 = 0.5*(w₁[1]+w₁[2])
    w₁_13 = 0.5*(w₁[1]+w₁[3])
    w₁_23 = 0.5*(w₁[2]+w₁[3])

    ϕx=∇ϕ(i,elem,mesh)[1]
    a1 = w₁_12*ϕ(j, p_med[1],elem,mesh)*ϕx
    a2 = w₁_13*ϕ(j, p_med[2],elem,mesh)*ϕx
    a3 = w₁_23*ϕ(j, p_med[3],elem,mesh)*ϕx

    int_1 = (a1+a2+a3)*abs(quad(elem,mesh))/6

    w₂_12 = 0.5*(w₂[1]+w₂[2])
    w₂_13 = 0.5*(w₂[1]+w₂[3])
    w₂_23 = 0.5*(w₂[2]+w₂[3])

    ϕy=∇ϕ(i,elem,mesh)[2]
    b1 = w₂_12*ϕ(j, p_med[1],elem,mesh)*ϕy
    b2 = w₂_13*ϕ(j, p_med[2],elem,mesh)*ϕy
    b3 = w₂_23*ϕ(j, p_med[3],elem,mesh)*ϕy

    int_2 = (b1+b2+b3)*abs(quad(elem,mesh))/6

    return int_1+int_2

end



function stiffness_convec(mesh,w₁,w₂)
    
    C = zeros(mesh.nv,mesh.nv)

    for t=1:mesh.nt
        elem = mesh.elements[t,:]
        M = zeros(3,3)
        for i=1:3
            for j=1:3
                M[i,j] = ∫convect(w₁,w₂,i,j,elem,mesh)
            end
        end
        C[elem,elem] = C[elem,elem] +M 
    end
    return C 
end

function ∫G(i,j,elem,mesh)
    p_med=midpoints(elem,mesh)
    Π=1/3
    a1=(ϕ(i, p_med[1],elem,mesh)-Π)*(ϕ(j, p_med[1],elem,mesh)-Π)
    a2=(ϕ(i, p_med[2],elem,mesh)-Π)*(ϕ(j, p_med[2],elem,mesh)-Π)
    a3=(ϕ(i, p_med[3],elem,mesh)-Π)*(ϕ(j, p_med[3],elem,mesh)-Π)
    
    return (abs(quad(elem,mesh))/6)*(a1+a2+a3)
end

function stiffness_G(mesh)

    G = zeros(mesh.nv,mesh.nv)
    
    for t=1:mesh.nt
        elem=mesh.elements[t,:]
        S = zeros(3,3)
        for i=1:3
            for j=1:3
                S[i,j]=∫G(i, j,elem,mesh)
            end
        end
        
        G[elem,elem] = G[elem,elem] + S
    end

    return G
end

### Rigth hand ###

function ∫f(f,i,elem,mesh)
    
    p_med=midpoints(elem,mesh)
    a1=f(p_med[1])*ϕ(i,p_med[1],elem,mesh)
    a2=f(p_med[2])*ϕ(i,p_med[2],elem,mesh)
    a3=f(p_med[3])*ϕ(i,p_med[3],elem,mesh)
    
    return abs(quad(elem,mesh)/6)*(a1+a2+a3)
end

function LoadVector(f,mesh)

    F = zeros(2*mesh.nv)
    
    for t=1:mesh.nt
        elem=mesh.elements[t,:]
        
        for i=1:3
            F[elem[i]]=F[elem[i]]+∫f(f,i,elem,mesh)[1]
            F[mesh.nv+elem[i]]=F[mesh.nv+elem[i]]+∫f(f,i,elem,mesh)[2]
        end
    end
    
    return [F;zeros(mesh.nv)]
end

function ∫p(ph,mesh)

    integral=0
    for t=1:mesh.nt
        elem=mesh.elements[t,:]
        integral=integral+(ph[elem[1]]*∫ϕ(1,elem,mesh)+ph[elem[2]]*∫ϕ(2,elem,mesh)+ph[elem[3]]*∫ϕ(3,elem,mesh))
    end

    return integral
end

### Plots functions ###

function contour_plot(triangulation,u)
    f, ax, tr  = tricontourf(triangulation, u,levels = 100,axis = (; aspect = 1))
    Colorbar(f[1, 2], tr)
    f
end

function quiver_plot(mesh,u₁,u₂)
    
    Plots.quiver(mesh.nodes[:,1],mesh.nodes[:, 2],quiver=(u₁/norm(u₁),u₂/norm(u₂)))

end

function suface_plot(mesh,u)

    Plots.surface(mesh.nodes[:,1],mesh.nodes[:,2],u)
end


### solvers ###

function solve_Stokes(f,u₀,mesh)
    

    #print("MMfem Library: Stokes Problem","\n")
    #print("\n")
    #print("Number of vertices: ", mesh.nv, "\n")
    #print("Number of triangles: ", mesh.nt,"\n")


    uₕ = zeros(2*mesh.nv)
    pₕ = zeros(mesh.nv)
    
    
    X_boundary=mesh.nodes[mesh.boundary,:]
    uₕ[mesh.boundary]=u₀(X_boundary,1) 
    uₕ[mesh.boundary.+mesh.nv]=u₀(X_boundary,2) 
    
    
    A = stiffness_Δ(mesh)
    Bx , By = stiffness_div(mesh)
    G = stiffness_G(mesh)
    
    M₀ = zeros(mesh.nv,mesh.nv)
    
    M = [A M₀ -Bx';M₀ A -By';-Bx -By -G]
    
    F₀ = LoadVector(f,mesh)
    
    NFree_u = setdiff(1:mesh.nv,mesh.boundary)
    nv_free_u = length(NFree_u)
    
    node_p = NFree_u[1]
    NFree_p = setdiff(1:mesh.nv,node_p)
    
    pₕ[node_p]=0

    F₀ = F₀-M*[uₕ;pₕ]
    
    K = [A[NFree_u,NFree_u] M₀[NFree_u,NFree_u] -Bx[NFree_p,NFree_u]';
        M₀[NFree_u,NFree_u] A[NFree_u,NFree_u] -By[NFree_p,NFree_u]';
        -Bx[NFree_p,NFree_u] -By[NFree_p,NFree_u] -G[NFree_p,NFree_p]]
    
 
    F1 = F₀[NFree_u]
    F2 = F₀[NFree_u.+mesh.nv]
    F3 = F₀[NFree_p.+2*mesh.nv]
    F = [F1;F2;F3]

    
    ũ = K\F
    
    uₕ[NFree_u]=ũ[1:nv_free_u]
    uₕ[NFree_u.+mesh.nv]=ũ[nv_free_u+1:2*nv_free_u]

    uₕ_1 = uₕ[1:mesh.nv]
    uₕ_2 = uₕ[mesh.nv+1:2*mesh.nv]

    pₕ[NFree_p]=ũ[2*nv_free_u+1:2*nv_free_u+mesh.nv-1]
    
    integral = ∫p(pₕ,mesh)
    pₕ=pₕ.-integral 
    
    return uₕ_1,uₕ_2,pₕ
end

function solve_Ossen(f,u₀,w₁,w₂,mesh)

    uₕ = zeros(2*mesh.nv)
    pₕ = zeros(mesh.nv)
    
    X_boundary=mesh.nodes[mesh.boundary,:]
    uₕ[mesh.boundary]=u₀(X_boundary,1) 
    uₕ[mesh.boundary.+mesh.nv]=u₀(X_boundary,2) 
    
    
    A₀ = stiffness_Δ(mesh)
    C = stiffness_convec(mesh,w₁,w₂)
    A = A₀+C
    Bx , By = stiffness_div(mesh)
    G = stiffness_G(mesh)
    
    M₀ = zeros(mesh.nv,mesh.nv)
    
    M = [A M₀ -Bx';M₀ A -By';-Bx -By -G]
    
    F₀ = LoadVector(f,mesh)
    
    NFree_u = setdiff(1:mesh.nv,mesh.boundary)
    nv_free_u = length(NFree_u)
    
    node_p = NFree_u[1]
    NFree_p = setdiff(1:mesh.nv,node_p)
    
    pₕ[node_p]=0

    F₀ = F₀-M*[uₕ;pₕ]
    
    K = [A[NFree_u,NFree_u] M₀[NFree_u,NFree_u] -Bx[NFree_p,NFree_u]';
        M₀[NFree_u,NFree_u] A[NFree_u,NFree_u] -By[NFree_p,NFree_u]';
        -Bx[NFree_p,NFree_u] -By[NFree_p,NFree_u] -G[NFree_p,NFree_p]]
    
 
    F1 = F₀[NFree_u]
    F2 = F₀[NFree_u.+mesh.nv]
    F3 = F₀[NFree_p.+2*mesh.nv]
    F = [F1;F2;F3]

    
    ũ = K\F
    
    uₕ[NFree_u]=ũ[1:nv_free_u]
    uₕ[NFree_u.+mesh.nv]=ũ[nv_free_u+1:2*nv_free_u]

    uₕ_1 = uₕ[1:mesh.nv]
    uₕ_2 = uₕ[mesh.nv+1:2*mesh.nv]

    pₕ[NFree_p]=ũ[2*nv_free_u+1:2*nv_free_u+mesh.nv-1]
    
    integral = ∫p(pₕ,mesh)
    pₕ=pₕ.-integral 
    
    return uₕ_1,uₕ_2,pₕ
end

function solve_NavierStokes(f,u₀,mesh;max_iter = 10,tol =1e-8)
    
    w₁,w₂,pₕ = solve_Stokes(f,u₀,mesh)

    for i=1:max_iter

        w₁,w₂,pₕ = solve_Ossen(f,u₀,w₁,w₂,mesh)

    end
    uₕ_1 = w₁
    uₕ_2 = w₂

    return uₕ_1,uₕ_2,pₕ

end


### Errors functions ###

function error_L2(mesh,u,uₕ)
    error = 0
    
    for t=1:mesh.nt
        elem = mesh.elements[t,:]
        
        p_med = midpoints(elem,mesh)
        
        uₕ_12=0.5*(uₕ[elem[1]]+uₕ[elem[2]])
        uₕ_23=0.5*(uₕ[elem[2]]+uₕ[elem[3]])
        uₕ_13=0.5*(uₕ[elem[1]]+uₕ[elem[3]])

        B=[(v₂-v₁) (v₃-v₁)]
        M=(abs(quad(elem,mesh))/6)*((u(p_med[1])-uₕ12)^2+(u(p_med[2])-uₕ13)^2+(u(p_med[3])-uₕ23)^2)
        
        error =error + M
    end
    return sqrt(error)
end

function search_node(node,mesh)
    
    for t = 1:mesh.nt
        elem = mesh.elements[t,:]
        if node == elem[1]
            return [t,1]
        elseif node == elem[2]
            return [t,2]
        elseif node == elem[3]
            return [t,3]
        end
    end
end

function ∇u(u,mesh)

    grad = zeros(mesh.nv,2)

    for i=1:mesh.nv
        elem , k = search_node(i,mesh)
        a₁ = i
        if k == 1
            a₂ = mesh.elements[elem,2]
            a₃ = mesh.elements[elem,3]
            r₁ = mesh.nodes[a₂,:]-mesh.nodes[a₁,:]
            r₂ = mesh.nodes[a₃,:]-mesh.nodes[a₁,:]
        elseif k == 2
            a₂ = mesh.elements[elem,1]
            a₃ = mesh.elements[elem,3]
            r₁ = mesh.nodes[a₂,:]-mesh.nodes[a₁,:]
            r₂ = mesh.nodes[a₃,:]-mesh.nodes[a₁,:]
        elseif k == 3
            a₂ = mesh.elements[elem,1]
            a₃ = mesh.elements[elem,2]
            r₁ = mesh.nodes[a₂,:]-mesh.nodes[a₁,:]
            r₂ = mesh.nodes[a₃,:]-mesh.nodes[a₁,:]
        end
        A = [r₁./norm(r₁) r₂./norm(r₂)]
        A = A'
        b = [(u[a₂]-u[a₁])/norm(r₁);(u[a₃]-u[a₁])/norm(r₂)]
        parcial = A\b
        grad[i,:] = parcial
    end
    
    return grad
end