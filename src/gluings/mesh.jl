export get_mesh
function get_mesh(domain::Gluing, h::Float64, verbose::Bool=false)
    verb(verbose, "Meshing domain $(domain.id) with step ", h, "...")

    # simple solution: mesh every component independantly of the junctions.
    submeshes = [get_mesh(comp, h) for comp in domain.comps]
    compstops = cumsum([m.npoints for m in submeshes])

    verb(verbose, "... done with $(compstops[end]) points.")
    return GMesh(compstops[end], maximum([m.step for m in submeshes]), [GPoint(icomp, pt.p) for (icomp,m) in enumerate(submeshes) for pt in m.points ], 
                ([0] ∪ compstops[1:end-1]) .+ 1, compstops, [m.layers for m in submeshes])
end

export search_in_mesh 
function search_in_mesh(domain::Gluing, mesh::GMesh, x::GPoint)
    return mesh.compinits[x.comp] + search_in_mesh(domain.comps[x.comp], mesh.complayers[x.comp], x.p) - 1
end

export get_mesh_convex
function get_mesh_convex(domain::Gluing, mesh::GMesh, points::Vector{GPoint}, verbose::Bool=false)

    if length(union([p.comp for p in points])) == 1 # only one component
        return @views mesh.compinits[points[1].comp] .- 1 .+ 
                    get_mesh_convex(domain.comps[points[1].comp], mesh.complayers[points[1].comp], 
                    mesh.points[mesh.compinits[points[1].comp]:mesh.compstops[points[1].comp]], points, verbose)

    else
        allvertices =  [cross for p1 in points for p2 in points for cross in domain.get_crossing(p1, p2)]
        comps = union([p.comp for p in allvertices], [p.comp for p in points])
        return union([get_mesh_convex(domain, mesh, [pt for pt in allvertices ∪ points if pt.comp == c], verbose) for c in comps]...)
    end
end