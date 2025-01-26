############################################################
# AFEM in 2d for the elliptic problem
#  -div(a∇u) + b⋅∇u + c u = f   in Ω
#       u = g_D  in Γ_D
#   ∂u/∂n = g_N  in Γ_N
############################################################

include("MM_adapt.jl")

#initial mesh
domain = "square_all_dirichlet"


mutable struct data_mesh
    nodes::Matrix
    elements::Matrix{Int64}
    neighbours::Matrix{Int64}
    boundaries::Matrix{Int64}
    nv::Int
    nt::Int
    mark::Vector{Int64}
end

struct Adapt

    tolerance::Float64 #tolerance for the adaptive strategy
    max_iterations::Int64 #maximum number of iterations of adaptive strategy
    strategy::String #marking_strategy
    n_refine::Int #number of refinements of each marked element
    MS_gamma::Float64 #parameters of the different marking strategies

end

vertex_coordinates,elem_vertices,elem_neighbours,elem_boundaries,n_vert,n_elem = inital_mesh(domain)
mesh = data_mesh(vertex_coordinates,elem_vertices,elem_neighbours,elem_boundaries,n_vert,n_elem,[1])
#init_tri = triangulate(mesh.nodes')
#display(triplot(init_tri))

#data
global_refinements = true
a = 1 
b = [0 0]
c = 0
#f(v) = 
#u(v) = 
#grad_u(u) = 




# resolve

#primer refinamiento
if global_refinements
    for i=1:10
        mesh.mark = Int.(ones(mesh.nt))
        refine_mesh(mesh)
        GMT.triplot(mesh.nodes,lc=:blue,show=true)
    end
end

#tri = triangulate(mesh.nodes')
#display(triplot(tri))
