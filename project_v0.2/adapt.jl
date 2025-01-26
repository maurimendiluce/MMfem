#creamos la estructura de la malla
using CairoMakie, DelaunayTriangulation

mutable struct data_mesh
    elem_vertices::Matrix{Int64}
    elem_neighbours::Matrix{Int64}
    elem_boundaries::Matrix{Int64}
    vertex_coordinates::Matrix{Float64}
    mark::Vector{Int64}
end

#mutable struct size_mesh
#    n_elem::Int64
#    n_vertices::Int64
#end

mesh=data_mesh([3 1 2;1 3 4],[0 0 2;0 0 1],[-1 1 0;1 -1 0],[0.0 0.0;1.0 0.0;1.0 1.0;0.0 1.0],[1,1])
tri_old = triangulate([0.0 0.0;1.0 0.0;1.0 1.0;0.0 1.0]')

#n_elem=size(mesh.elem_vertices,1)
#n_vertices=size(mesh.vertex_coordinates,1)

#uh=zeros(n_elem,1)
#fh=zeros(n_elem,1)

function refine_element(elem_to_refine,mesh)
    
    n_refined = 0
    if maximum(mesh.mark)==0
        #si no hay elementos marcados no hace nada.
        return
    end

    vecino = mesh.elem_neighbours[elem_to_refine, 3]

    if vecino > 0 #no es lado del borde
        if mesh.elem_neighbours[vecino, 3] != elem_to_refine  #no compatible
            n_refined = n_refined + refine_element(vecino)
        end
    end
    
    n_elem=size(mesh.elem_vertices,1)
    n_vertices=size(mesh.vertex_coordinates,1)

    #ahora el vecino es compatible y viene el refinamiento
    #la informacion del vecino puede haber cambiado
    vecino = mesh.elem_neighbours[elem_to_refine, 3]

    v_elem = mesh.elem_vertices[elem_to_refine,:]
    neighs_elem = mesh.elem_neighbours[elem_to_refine,:] 
    bdries_elem = mesh.elem_boundaries[elem_to_refine,:] 

    n_vertices = n_vertices + 1
    aux=0.5*(mesh.vertex_coordinates[v_elem[1],:] + mesh.vertex_coordinates[v_elem[2],:])
    aux=hcat(aux...)
    mesh.vertex_coordinates =vcat(mesh.vertex_coordinates,aux)
    #uh[n_vertices, :] = 0.5*(uh[v_elem[1], :] + uh[v_elem[2], :])
    #fh[n_vertices, :] = 0.5*(fh[v_elem[1], :] + fh[v_elem[2], :])
    
    if vecino == 0 #si el marcado es un lado de borde
        n_elem = n_elem + 1
        mesh.elem_vertices[elem_to_refine, :] = [v_elem[3] v_elem[1] n_vertices]
        mesh.elem_vertices =vcat(mesh.elem_vertices,[v_elem[2] v_elem[3] n_vertices])

        #actualiza la informacion del vecino
        mesh.elem_neighbours[elem_to_refine, :] = [neighs_elem[3] n_elem neighs_elem[2]]
        mesh.elem_neighbours = vcat(mesh.elem_neighbours,[elem_to_refine neighs_elem[3] neighs_elem[1]])
        #actualiza la informacion de borde
        mesh.elem_boundaries[elem_to_refine, :] = [bdries_elem[3] 0 bdries_elem[2]]
        mesh.elem_boundaries = vcat(mesh.elem_boundaries,[0 bdries_elem[3] bdries_elem[1]])

        #actualiza la informacion del vecino en los datos sobre los
        #vecinos
        if neighs_elem[1] > 0
            f=findfirst(x->x==elem_to_refine,mesh.elem_neighbours[neighs_elem[1],:])
            mesh.elem_neighbours[neighs_elem[1],f] = n_elem
        end

        mesh.mark=vcat(mesh.mark,maximum(mesh.mark[elem_to_refine] - 1, 0))
        mesh.mark[elem_to_refine] = maximum(mesh.mark[elem_to_refine] - 1, 0)
        n_refined = n_refined + 1;
    else
        v_neigh = mesh.elem_vertices[vecino,:]
        neighs_neigh = mesh.elem_neighbours[vecino,:] 
        bdries_neigh = mesh.elem_boundaries[vecino,:] 

        n_elem = n_elem + 2
        mesh.elem_vertices[elem_to_refine, :] = [v_elem[3]  v_elem[1] n_vertices]
        mesh.elem_vertices=vcat(mesh.elem_vertices,[v_elem[2]  v_elem[3] n_vertices])
        #mesh.elem_vertices(n_elem-1, :)       = [v_elem(2)  v_elem(3) n_vertices];
        mesh.elem_vertices[vecino, :]      = [v_neigh[3] v_neigh[1] n_vertices]
        #mesh.elem_vertices(n_elem, :)         = [v_neigh(2) v_neigh(3) n_vertices];
        mesh.elem_vertices=vcat(mesh.elem_vertices,[v_neigh[2] v_neigh[3] n_vertices])

        #actualiza la informacion de los vecinos
        mesh.elem_neighbours[elem_to_refine, :] = [n_elem n_elem-1 neighs_elem[2]]
        mesh.elem_neighbours=vcat(mesh.elem_neighbours,[elem_to_refine vecino neighs_elem[1]])
        #mesh.elem_neighbours(n_elem-1, :) = [elem_to_refine vecino neighs_elem(1)];
        mesh.elem_neighbours[vecino, :] = [n_elem-1 n_elem neighs_neigh[2]]
        mesh.elem_neighbours=vcat(mesh.elem_neighbours,[vecino elem_to_refine neighs_neigh[1]])
        #mesh.elem_neighbours(n_elem, :) = [vecino elem_to_refine neighs_neigh(1)];
        
        #actualiza la informacion del borde
        mesh.elem_boundaries[elem_to_refine, :] = [0 0 bdries_elem[2]]
        mesh.elem_boundaries=vcat(mesh.elem_boundaries,[0 0 bdries_elem[1]])
        #mesh.elem_boundaries(n_elem-1, :)       = [0 0 bdries_elem(1)];
        mesh.elem_boundaries[vecino, :]      = [0 0 bdries_neigh[2]]
        #mesh.elem_boundaries(n_elem, :)         = [0 0 bdries_neigh(1)];
        mesh.elem_boundaries=vcat(mesh.elem_boundaries,[0 0 bdries_neigh[1]])

        #update neighbour information of the neighbours
        if neighs_elem[1] > 0
            #f = find(mesh.elem_neighbours(neighs_elem(1),:) == elem_to_refine);
            f=findfirst(x->x==elem_to_refine,mesh.elem_neighbours[neighs_elem[1],:])
            mesh.elem_neighbours[neighs_elem[1],f] = n_elem-1
        end
        if neighs_neigh[1] > 0
            #f = find(mesh.elem_neighbours(neighs_neigh(1),:) == vecino);
            f=findfirst(x->x==vecino,mesh.elem_neighbours[neighs_neigh[1],:])
            mesh.elem_neighbours[neighs_neigh[1],f] = n_elem
        end

        #mesh.mark(n_elem-1) = max(mesh.mark(elem_to_refine) - 1, 0);
        mesh.mark=vcat(mesh.mark,maximum([mesh.mark[elem_to_refine] - 1, 0]))
        mesh.mark[elem_to_refine] = maximum([mesh.mark[elem_to_refine] - 1, 0])

        #mesh.mark(n_elem) = max(mesh.mark(vecino) - 1, 0);
        mesh.mark=vcat(mesh.mark,maximum([mesh.mark[vecino] - 1, 0]))
        mesh.mark[vecino] = maximum([mesh.mark[vecino] - 1, 0])
        n_refined = n_refined + 2
    end

    return n_refined,mesh
    #mesh.n_elem = n_elem
    #mesh.n_vertices = n_vertices;
    
end

function refine_mesh(mesh)
    
    n_refined = 0

    if maximum(mesh.mark)==0
      #no hay elementos marcados, termina.
      return
    end
    
    while maximum(mesh.mark) > 0
        first_marked = minimum(findfirst(x->x>0,mesh.mark))
        N,mesh=refine_element(first_marked,mesh)
        n_refined = n_refined + N
    end
    return n_refined,mesh
end



#len=size_mesh(size(mesh.elem_vertices,1),size(mesh.elem_vertices,1))
N,mesh=refine_mesh(mesh)
tri_new = triangulate(mesh.vertex_coordinates') 
triplot(tri_new)