export plot_function
function plot_function(domain::Domain_euclidean, mesh::Mesh, values::VecOrMat{Float64}, dim3::Bool=false, verbose::Bool=false)
    verb(verbose, "Plotting on $(domain.id)...")

    if (domain.dim == 1)
        p = plot([p.p[1] for p in mesh.points], values[:], legend=false, aspect_ratio=:equal)
        display(p)

    elseif (domain.dim == 2)
        
        if dim3
            p = scatter3d([p.p[1] for p in mesh.points], [p.p[2] for p in mesh.points], values[:], marker_z=values[:], 
                xaxis="x", yaxis="y", aspect_ratio=:equal, ms=0.5, legend=false, camera=(-40,0))
            display(p)
        else

            p = scatter([p.p[1] for p in mesh.points], [p.p[2] for p in mesh.points], marker_z=values[:], 
                xaxis="x", yaxis="y", legend=false, colorbar=true, aspect_ratio=:equal, ms=2, markershape=:square)
            display(p)
        end
    else
        throw(ErrorException("dimension $(domain.dim) not implemented in $(fname)"))
    end

    verb(verbose, "... done.")
    return p
end