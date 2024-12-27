# MMfem: Finite Element Library

Welcome to the documentation page. 

!!! note "Note"
    Page under construction

# Mesh generate
```@docs
rectangle_mesh(a,b,c,d,N)
nodes_tri(tri)
elements_tri(tri)
boundary_nodes(tri)
mesh_data(tri)
```

# Plot functions
```@docs
plot_mesh(tri)
scalar_plot(u,tri;title = "titulo")
vectorial_plot(F,points)
```

# Base Functions
```@docs
β(j,v,points,triangle)
∇β(j,points,triangle)
```

# Cuadratures
```@docs
quad(f,points,triangle;type="midpoints")
```

# Operators
## Stiffnes_Matrix
```@docs
Stiffnes_Δ(tri)
```
## Load Vectors
```@docs
fβ(f,j,v,points,triangle)
LoadVector(f,tri)
```

# Boundary conditions
```@docs
dirichlet(tri,dirichlet_nodes,uD)
```

# Solver 
```@docs
solver_problem(A,F)
```


