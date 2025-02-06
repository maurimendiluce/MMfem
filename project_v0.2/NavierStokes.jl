include("MM_mixedfem.jl")

print("MMfem Library: Navier-Stokes Problem","\n")
#print("\n")
print("------------------------------------", "\n")


#mesh

mesh,triangulation = create_mesh(0,1,0,1,32,graph=false)

#data
f(v)=[0,0]
μ = 1

## boundary date ###
function u₀(X,component)
    u = zeros(length(X[:,1]))
    if component == 1
        for j=1:length(X[:,1])
            if X[:,2][j] == 1
                u[j] = 1
            elseif X[:,2][j] == 0
                u[j] = 0
            elseif X[:,1][j] == 1
                u[j] = 0
            elseif X[:,1][j] == 0
                u[j] = 0
            end
        end
        return u
    else
        for j=1:length(X[:,1])
            if X[:,2][j] == 1
                u[j]=0
            elseif X[:,2][j] == 0
                u[j]=0
            elseif X[:,1][j] == 1
                u[j] = 0
            elseif X[:,1][j] == 0
                u[j]=0
            end
        end
        return u
    end
end

uₕ_1,uₕ_2,pₕ = solve_NavierStokes(f,u₀,mesh)


#Plots

trisurf([mesh.nodes[:,1] mesh.nodes[:,2] pₕ], show=true)
#quiver(mesh.nodes[:,1], mesh.nodes[:,2], uₕ_1, uₕ_2, fill=:black, lc=:red, show=true)

#contour_plot(triangulation,pₕ)
#quiver_plot(mesh,uₕ_1,uₕ_2)
#surface_plot(mesh,pₕ)