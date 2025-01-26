#Packege

using DelaunayTriangulation, CairoMakie, LinearAlgebra, Plots, GMT

# initial meshes #


function inital_mesh(type)

    if type == "square_all_dirichlet"
    
        nodes = [0.0  0.0;1.0  0.0;1.0  1.0;0.0  1.0]
        elem_vertices =  [3  1  2;1  3  4]
        elem_neighbours = [0  0  2;0  0  1]
        elem_boundaries =  [1  1  0;1  1  0]
        nv = size(nodes)[1]
        nt = size(elem_vertices)[1]
    end

    return nodes,elem_vertices,elem_neighbours,elem_boundaries,nv,nt

end


function refine_element(elem_to_refine,mesh)
    
    n_refined = 0

    neighbour = mesh.neighbours[elem_to_refine, 3]

    if neighbour> 0 #no es lado del borde
        if mesh.neighbours[neighbour, 3] != elem_to_refine  #no compatible
            n_refined = n_refined + refine_element(neighbour,mesh)
        end
    end
    
    n_elem=mesh.nt
    n_vertices=mesh.nv

    #ahora el vecino es compatible y viene el refinamiento
    #la informacion del vecino puede haber cambiado
    neighbour = mesh.neighbours[elem_to_refine, 3]

    v_elem = mesh.elements[elem_to_refine,:]
    neighs_elem = mesh.neighbours[elem_to_refine,:] 
    bdries_elem = mesh.boundaries[elem_to_refine,:] 

    n_vertices = n_vertices + 1
    mesh.nodes[n_vertices,:] = 0.5*(mesh.nodes[v_elem[1],:] + mesh.nodes[v_elem[2],:])
    #uh[n_vertices, :] = 0.5*(uh[v_elem[1], :] + uh[v_elem[2], :])
    #fh[n_vertices, :] = 0.5*(fh[v_elem[1], :] + fh[v_elem[2], :])
    
    if neighbour == 0 #si el marcado es un lado de borde
        n_elem = n_elem + 1
        mesh.elements[elem_to_refine, :] = [v_elem[3] v_elem[1] n_vertices]
        mesh.elements[n_elem,:] = [v_elem[2] v_elem[3] n_vertices]

        #actualiza la informacion del vecino
        mesh.neighbours[elem_to_refine, :] = [neighs_elem[3] n_elem neighs_elem[2]]
        mesh.neighbours[n_elem,:] = [elem_to_refine neighs_elem[3] neighs_elem[1]]
        #actualiza la informacion de borde
        mesh.boundaries[elem_to_refine, :] = [bdries_elem[3] 0 bdries_elem[2]]
        mesh.boundaries[n_elem,:] = [0 bdries_elem[3] bdries_elem[1]]

        #actualiza la informacion del vecino en los datos sobre los
        #vecinos
        if neighs_elem[1] > 0
            f=findfirst(x->x==elem_to_refine,mesh.neighbours[neighs_elem[1],:])
            mesh.neighbours[neighs_elem[1],f] = n_elem
        end

        mesh.mark[n_elem] = maximum([mesh.mark[elem_to_refine] - 1, 0])
        mesh.mark[elem_to_refine] = maximum([mesh.mark[elem_to_refine] - 1, 0])
        n_refined = n_refined + 1;
    else
        v_neigh = mesh.elements[neighbour,:]
        neighs_neigh = mesh.neighbours[neighbour,:] 
        bdries_neigh = mesh.boundaries[neighbour,:] 

        n_elem = n_elem + 2
        mesh.elements[elem_to_refine, :] = [v_elem[3]  v_elem[1] n_vertices]
        mesh.elements[n_elem-1,:] = [v_elem[2]  v_elem[3] n_vertices]
        mesh.elements[neighbour, :] = [v_neigh[3] v_neigh[1] n_vertices]
        mesh.elements[n_elem,:] = [v_neigh[2] v_neigh[3] n_vertices]

        #actualiza la informacion de los vecinos
        mesh.neighbours[elem_to_refine, :] = [n_elem n_elem-1 neighs_elem[2]]
        mesh.neighbours[n_elem-1,:] = [elem_to_refine neighbour neighs_elem[1]]
        mesh.neighbours[neighbour, :] = [n_elem-1 n_elem neighs_neigh[2]]
        mesh.neighbours[n_elem,:] = [neighbour elem_to_refine neighs_neigh[1]]
        
        #actualiza la informacion del borde
        mesh.boundaries[elem_to_refine, :] = [0 0 bdries_elem[2]]
        mesh.boundaries[n_elem-1,:] = [0 0 bdries_elem[1]]
        mesh.boundaries[neighbour, :] = [0 0 bdries_neigh[2]]
        mesh.boundaries[n_elem,:] = [0 0 bdries_neigh[1]]

        #update neighbour information of the neighbours
        if neighs_elem[1] > 0
            f=findfirst(x->x==elem_to_refine,mesh.neighbours[neighs_elem[1],:])
            mesh.neighbours[neighs_elem[1],f] = n_elem-1
        end
        if neighs_neigh[1] > 0
            f=findfirst(x->x==neighbour,mesh.neighbours[neighs_neigh[1],:])
            mesh.neighbours[neighs_neigh[1],f] = n_elem
        end

        mesh.mark[n_elem-1] = maximum([mesh.mark[elem_to_refine] - 1, 0])
        mesh.mark[elem_to_refine] = maximum([mesh.mark[elem_to_refine] - 1, 0])

        mesh.mark[n_elem] = maximum([mesh.mark[neighbour] - 1, 0])
        mesh.mark[neighbour] = maximum([mesh.mark[neighbour] - 1, 0])
        n_refined = n_refined + 2
    end

    mesh.nt = n_elem
    mesh.nv = n_vertices;
    return n_refined
    
end

function refine_mesh(mesh)
    
    if (maximum(mesh.mark)==0)
        #no elements marked, doing nothing
        return
    end

    if true
        #we first extend the matrices and vectors to save 
        # time in memory allocation
        # estimation of final number of elements
        nelem_new = mesh.nt*2*maximum(mesh.mark)
        mesh.elements = vcat(mesh.elements,zeros(nelem_new-size(mesh.elements)[1],3))
        mesh.neighbours = vcat(mesh.neighbours,zeros(nelem_new-size(mesh.neighbours)[1],3))
        mesh.boundaries = vcat(mesh.boundaries,zeros(nelem_new-size(mesh.boundaries)[1],3))
        mesh.mark = vcat(mesh.mark,zeros(nelem_new-length(mesh.mark)))
        
        #estimation of final number of vertices
        # improve this computation
        nvertices_new = 4*mesh.nv
        mesh.nodes = vcat(mesh.nodes,zeros(nvertices_new-size(mesh.nodes)[1],2))
        #uh(nvertices_new,1) = 0;
        #fh(nvertices_new,1) = 0;

    end

    n_refined = 0

    while maximum(mesh.mark) > 0
        first_marked = minimum(findall(x -> x>0,mesh.mark))
        n_refined = n_refined + refine_element(first_marked,mesh);
    end

    if true
        #we now clean the unused space
        
        mesh.elements = mesh.elements[1:mesh.nt,:]
        mesh.neighbours = mesh.neighbours[1:mesh.nt,:]
        mesh.boundaries = mesh.boundaries[1:mesh.nt,:]
        mesh.mark = mesh.mark[1:mesh.nt]

        mesh.nodes = mesh.nodes[1:mesh.nv,:]
        
        #uh(mesh.n_vertices+1:nvertices_new,:) = [];
        #fh(mesh.n_vertices+1:nvertices_new,:) = [];
    end

    #return mesh
end

