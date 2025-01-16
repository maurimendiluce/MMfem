include("MMfem.jl")

mesh = rectangle_mesh(-1,1,-1,1,16)
#plot_mesh(mesh)

#Funcion exacta
U(x,y)=(x-0.5*x^2)*(y-0.5*y^2)
U(v)=U(v[1],v[2])

#dato de borde
function g(x,y)
    if x==-1
        return -(3/2)*(y-0.5*y^2)
    elseif y==-1
        return -(3/2)*(x-0.5*x^2)
    else
        return 0
    end
end
g(v) = g(v[1],v[2])

#fuente
f(x,y) = y-0.5*y^2+x-0.5*x^2
f(v) = f(v[1],v[2])


p,t,b = mesh_data(mesh)
uₕ = zeros(size(p,1),1) #vector solucion

#Assamble
A = Stiffnes_Δ(mesh)
F = LoadVector(f,mesh)

NodesFree = setdiff(1:size(p,1),b)

#Incorporación dato Dirichlet
for i=1:size(b,1)
    v=p[b[i],:]
    uₕ[b[i]]=g(v)
end 

F = F - A*uₕ

#Solve
A = Matrix(A)
F = Matrix(F)
uₕ[NodesFree] = A[NodesFree,NodesFree]\F[NodesFree]
uₕ = vec(uₕ)
∇uₕ = grad_calculate(uₕ,mesh)

#Graficos
p1 = scalar_plot(uₕ,mesh;title = "title")
p2 = vectorial_plot(-∇uₕ,p)
display(p1)
display(p2)
