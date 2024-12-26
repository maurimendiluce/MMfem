module MMfem

using DelaunayTriangulation, CairoMakie

#Funciones para generar la malla
export rectangle_mesh
"""
    rectangle_mesh

This function generates a mesh in the rectangle [a,b]x[c,d] using N subdivisions in each interval.
"""
function rectangle_mesh(a,b,c,d,N;graficar=true)

    tri = triangulate_rectangle(a, b, c, d, N, N,delete_ghosts = true)

    if graficar == true
        display(triplot(tri))
    end
    return tri
end

export nodes_tri
"""
    nodes_tri

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
    elements_tri

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
    boundary_nodes

This function returns the numbering of the edge nodes.
"""
function boundary_nodes(tri)
    #función que devuelve la numeración de los nodos de borde

    get_boundary = get_boundary_nodes(tri)
    boundary = [get_boundary[1];get_boundary[2];get_boundary[3];get_boundary[4]]    
    
    return unique(boundary)
end

#Funciones para graficar
export scalar_plot
"""
    scalar_plot

This function graphs a scalar function.
"""
function scalar_plot(u,tri;title = "titulo")
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

end # module MMfem
