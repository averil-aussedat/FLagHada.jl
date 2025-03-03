export get_mesh
"""
$(SIGNATURES)

Mesh the domain ``D`` with points ``(x_k)_{k \\in K}`` such that 
```math
    \\sup_{x \\in D} \\inf_{k \\in K} d(x, x_k) \\leqslant x_{\\text{step}}.
```

"""
function get_mesh(domain::Domain_euclidean, step::Float64, verbose::Bool=false)
    verb(verbose, "Meshing $(domain.id) with step $step...")
    
    # expected : bounds contains two points a and b with a .< b
    layers = round.(Int, (domain.bounds[2].p .- domain.bounds[1].p) ./ step) .+ 2
    step = maximum((domain.bounds[2].p .- domain.bounds[1].p) ./ (layers .- 1))

    # first dimension varies slower 
    points = Vector{RdPoint}([RdPoint(zeros(domain.dim)) for _ = 1:prod(layers)])

    for dd=1:domain.dim
        pattern = collect(LinRange(domain.bounds[1].p[dd], domain.bounds[2].p[dd], layers[dd]))
        pointctr = 1
        for outeri ∈ 1:prod(layers[1:dd-1])
            for pp ∈ pattern
                for inneri ∈ 1:prod(layers[dd+1:end])
                    points[pointctr].p[dd] = pp
                    pointctr += 1
                end
            end
        end
    end

    verb(verbose, "... done with $(prod(layers)) points.")
    return EuMesh(prod(layers), step, points, layers)
end

export search_in_mesh 
function search_in_mesh(domain::Domain_euclidean, layers::Vector{Int}, p::Vector{Float64})
    xlayers = 1 .+ round.(Int, (layers .- 1) .* (p .- domain.bounds[1].p) ./ (domain.bounds[2].p .- domain.bounds[1].p))
    xlayers = max.(1, min.(layers, xlayers))
    return 1 + sum((xlayers[i]-1) * prod(layers[i+1:end]) for i in 1:domain.dim)
end

function search_in_mesh(domain::Domain_euclidean, mesh::EuMesh, x::RdPoint)
    return search_in_mesh(domain, mesh.layers, x.p)
end

export get_mesh_convex
function get_mesh_convex(domain::Domain_euclidean, meshlayers::Vector{Int}, meshpoints, points::Vector{<:EuPoint}, verbose::Bool=false)::Vector{Int}

    if length(points) == 1
        return [search_in_mesh(domain, meshlayers, points[1].p)]

    elseif domain.dim == 1
        indmin = 1 + round(Int, (meshlayers[1] - 1) * (minimum([pt.p[1] for pt in points]) - domain.bounds[1].p[1]) / (domain.bounds[2].p[1] - domain.bounds[1].p[1]))
        indmax = 1 + round(Int, (meshlayers[1] - 1) * (maximum([pt.p[1] for pt in points]) - domain.bounds[1].p[1]) / (domain.bounds[2].p[1] - domain.bounds[1].p[1]))
        # println("indmin : ", indmin, ", indmax : ", indmax, ", interval : ", max(1,indmin):min(length(meshpoints),indmax))
        return collect(max(1,indmin):min(length(meshpoints),indmax))

    elseif domain.dim == 2
        meshproj = union([search_in_mesh(domain, meshlayers, pt.p) for pt in points])
        polygon = VPolygon(convex_hull([meshpoints[mp].p for mp in meshproj]))
        lowerinds = round.(Int, 1 .+ (meshlayers .- 1) .* (low(polygon)  .- domain.bounds[1].p) ./ (domain.bounds[2].p .- domain.bounds[1].p))
        upperinds = round.(Int, 1 .+ (meshlayers .- 1) .* (high(polygon) .- domain.bounds[1].p) ./ (domain.bounds[2].p .- domain.bounds[1].p))
        boxindexes = [1 + (ix-1) * meshlayers[2] + (iy-1) for ix in lowerinds[1]:upperinds[1] for iy in lowerinds[2]:upperinds[2]]
        return meshproj ∪ [ipt for ipt in boxindexes if meshpoints[ipt].p ∈ polygon]

    else 
        meshproj = union([search_in_mesh(domain, meshlayers, pt.p) for pt in points])
        polytope = VPolytope(convex_hull([meshpoints[mp].p for mp in meshproj]))
        lowerinds = round.(Int, 1 .+ (meshlayers .- 1) .* (low(polytope)  .- domain.bounds[1].p) ./ (domain.bounds[2].p .- domain.bounds[1].p))
        upperinds = round.(Int, 1 .+ (meshlayers .- 1) .* (high(polytope) .- domain.bounds[1].p) ./ (domain.bounds[2].p .- domain.bounds[1].p))
        # heavy, I know, but always less than checking the whole mesh... I assume
        boxindexes = [1 + sum((boxpt[i] - 1) * prod(meshlayers[i+1:end]) for i in 1:domain.dim) for boxpt in 
            zip([repeat(lowerinds[i]:upperinds[i], inner=prod(meshlayers[i+1:end]), outer=prod(meshlayers[1:i-1])) for i in 1:domain.dim]...)
        ]

        # if maximum(boxindexes) > length(meshpoints)
        #     println("meshlayers : ", meshlayers)
        #     println("meshproj : ", meshproj)
        #     println("lowerinds : ", lowerinds)
        #     println("upperinds : ", upperinds)
        #     # println("boxiterator : ", boxiterator, ", of type ", typeof(boxiterator))
        #     println("boxindexes : ", boxindexes)
        # end

        return meshproj ∪ [ipt for ipt in boxindexes if meshpoints[ipt].p ∈ polytope]
    end
end

function get_mesh_convex(domain::Domain_euclidean, mesh::EuMesh, points::Vector{RdPoint}, verbose::Bool=false)::Vector{Int}
    return get_mesh_convex(domain, mesh.layers, mesh.points, points, verbose)
end