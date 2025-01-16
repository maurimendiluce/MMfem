include("MMfem.jl")

#genero la malla
mesh = rectangle_mesh(0,1,0,1,16)
p,t,b = mesh_data(mesh) 


#datos
μ = 1
f1(v) = 0
f2(v) = 0

function g(x,y)
    if y==1 && 0<x<1
        return 1
    else
        return 0
    end
end
g(v) = g(v[1],v[2])

#reservo el espacio de solucion
uₕ = zeros(2*size(p,1),1)
pₕ = zeros(size(p,1),1)

NodesFree = setdiff(1:size(p,1),b)

#incorporo dato de borde
for i=1:size(b,1)
    v=p[b[i],:]
    uₕ[b[i]]=g(v)
end 

#Matrix Assamble

A = μ*Stiffnes_Δ(mesh)
Bx,By = Stiffnes_div(mesh)
G = Stiffnes_G(mesh)

K = [A spzeros(size(p,1),size(p,1)) -Bx';spzeros(size(p,1),size(p,1)) A -By';-Bx -By -G]
#uₕ_1 = uₕ[1:size(p,1)]
#uₕ_2 = uₕ[size(p,1)+1:2*size(p,1)]

F1 = LoadVector(f1,mesh)
F2 = LoadVector(f2,mesh)
F = [F1;F2;spzeros(size(p,1))]

#Fijamos nodo de la presión para incorporar promedio 0
node_p = NodesFree[1]
pₕ[node_p] = 0

F = F - K*[uₕ;pₕ]