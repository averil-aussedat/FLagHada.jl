export pt_to_str
function pt_to_str(pt::GPoint)
    return "C$(pt.comp)=>[" * join(numtostr.(pt.p, 5), ",") * "]"
end

export plot_gluing
"""
    For obvious reasons, we only consider 
    plots for gluings of dimensions <= 3.
    For plotting reasons, get_plotting_coords
    must return coordinates in dim ∈ [2,3].
"""
function plot_gluing(domain::Gluing, plotter::Gluing_plotter; lw=:auto, ms=:auto)

    vertices = Vector{Vector{Float64}}() # list of points 
    edges = Vector{Vector{Vector{Float64}}}() # list of pairs of points

    for (icomp, comp) in enumerate(domain.comps) 

        bd = [plotter.get_plotting_coords(domain, icomp, p) for p in comp.bounds]

        if comp.dim == 1
            append!(vertices, bd)
            append!(edges, [bd])

        elseif comp.dim == 2
            upperleft = plotter.get_plotting_coords(domain, icomp, RdPoint([comp.bounds[1].p[1],comp.bounds[2].p[2]]))
            lowerrigh = plotter.get_plotting_coords(domain, icomp, RdPoint([comp.bounds[2].p[1],comp.bounds[1].p[2]]))
            [append!(vertices, [bb]) for bb in [bd[1], upperleft, lowerrigh, bd[2]]]
            for (i,j) in [(1,2),(1,3),(2,4),(3,4)]
                append!(edges, [[vertices[end-4+i],vertices[end-4+j]]])
            end
            
        elseif comp.dim == 3 # so plotter.dim = 3
            [append!(vertices, [[bd[i][1],bd[j][2],bd[k][3]]]) for (i,j,k) in [(1,1,1),(1,1,2),(1,2,1),(1,2,2),(2,1,1),(2,1,2),(2,2,1),(2,2,2)]]
            [append!(edges, [[a,b]]) for a in vertices[end-7:end] for b in vertices[end-7:end] if sum(a .≈ b) == 2]

        else
            throw(ErrorException("can't plot component $(comp.id) of dimension $(comp.dim) > 3"))
        end

    end

    p = plot(legend=false, aspect_ratio=1.0)

    if plotter.dim == 2
        minx = minimum([v[1] for v in vertices])
        miny = minimum([v[2] for v in vertices])
        maxx = maximum([v[1] for v in vertices])
        maxy = maximum([v[2] for v in vertices])
        if ms == :auto 
            ms = 3
        end
        scatter!(p, [v[1] for v ∈ vertices], [v[2] for v ∈ vertices], xlims=[minx-0.5,maxx+0.5], ylims=[miny-0.5,maxy+0.5], color=:black, ms=ms)
        for edge in edges
            plot!(p, [edge[1][1],edge[2][1]], [edge[1][2],edge[2][2]], color=:black, lw=lw)
        end
    elseif plotter.dim == 3
        mins = minimum(minimum(vertices) .- 0.5)
        maxs = maximum(maximum(vertices) .+ 0.5)
        if ms == :auto
            ms = 1.5
        end
        scatter3d!(p, [v[1] for v ∈ vertices], [v[2] for v ∈ vertices], [v[3] for v ∈ vertices], xlims=[mins,maxs], ylims=[mins,maxs], zlims=[mins,maxs], color=:black, ms=ms)
        for edge in edges
            plot!(p, [edge[1][1],edge[2][1]], [edge[1][2],edge[2][2]], [edge[1][3],edge[2][3]], color=:black, lw=lw)
        end
        plot!(camera=(-30,40), axis=false)
    else
        throw(ErrorException("can't plot in dimension $(plotter.dim)"))
    end

    return p
end

export plot_points!
function plot_points!(p, domain::Gluing, plotter::Gluing_plotter, points::Vector{GPoint}; color=:blue, ms=1.5)

    toplot = [plotter.get_plotting_coords(domain, p.comp, p) for p in points]

    if plotter.dim == 2
        scatter!(p, [tp[1] for tp in toplot], [tp[2] for tp in toplot], color=color, ms=ms)
    else
        scatter3d!(p, [tp[1] for tp in toplot], [tp[2] for tp in toplot], [tp[3] for tp in toplot], color=color, ms=ms)
    end
end

export plot_values!
function plot_values!(p, domain::Gluing, plotter::Gluing_plotter, mesh::GMesh, values::Vector{Float64}, ms::Float64=2.0, verbose::Bool=false)
    minv = minimum(values)
    maxv = maximum(values)
    if maxv <= minv + 1e-8
        colors = :black 
        verb(verbose, "Constant values, falling back to default color")
    else 
        colors = RGB[cgrad(:inferno)[0.8 * (v - minv) / (maxv - minv)] for v ∈ values]
    end

    toplot = [plotter.get_plotting_coords(domain, p.comp, p) for p in mesh.points]
    if plotter.dim == 2
        scatter!(p, [tp[1] for tp in toplot], [tp[2] for tp in toplot], mc=colors, msc=colors, ms=ms, shape=:rect)
    else
        if domain.id == "ISS"
            toplot .= vcat(toplot[mesh.compstops[end]:-1:mesh.compinits[end]], toplot[begin:mesh.compinits[end]-1])
            colors .= vcat(colors[mesh.compstops[end]:-1:mesh.compinits[end]], colors[begin:mesh.compinits[end]-1])
        end

        scatter3d!(p, [tp[1] for tp in toplot], [tp[2] for tp in toplot], [tp[3] for tp in toplot], mc=colors, msc=colors, ms=ms)
        plot!(p, title = "Value between $(numtostr(minv,6)) and $(numtostr(maxv,6))")
        verb(verbose, "value between $minv and $maxv")
    end
end

export plot_feedback_map!
function plot_feedback_map!(p, domain::Gluing, plotter::Gluing_plotter, dynamic::MathFunction{Dynamic}, value, h::Float64, T::Float64, verbose::Bool=false)

    if domain.id == "planeline"
        meshstep = 0.4
        larr = 0.2
        arrlw = 2.0
    else # if domain.id == "ISS"
        meshstep = 0.6
        larr = 0.25
        arrlw = 1.0
    end

    # independent submesh 
    hmesh = get_mesh(domain, meshstep, true)
    directions = [get_feedback_direction(domain, value, dynamic.call(pt), pt, h, T) for pt in hmesh.points] 
    nextcrossd = Vector{GPoint}(undef, length(directions))
    for (ipt, pt,dir) in zip(1:hmesh.npoints, hmesh.points, directions)
        if dir.comp == pt.comp 
            nextcrossd[ipt] = dir
        else
            nextcrossd[ipt] = domain.get_crossing(pt, dir)[1]
        end
    end
    coords = [plotter.get_plotting_coords(domain, pt.comp, pt) for pt in hmesh.points]
    arrows = [plotter.get_plotting_coords(domain, ncd.comp, ncd) .- co for (co, ncd) in zip(coords, nextcrossd)]
    for arr in arrows 
        if norm(arr) > 1e-6
            arr .*= larr / norm(arr)
        end
    end

    col = :orangered
    for (co,arr) in zip(coords, arrows)
        arrow0!(p, co..., arr...; lc=col, ec=2.0, as=0.3, lw=arrlw)
        if norm(arr) <= 1e-5
            if domain.id == "planeline"
                scatter!(p, [[co[i]] for i in 1:2]..., color=col, msc=col, msw=0.3, ms=0.8, lw=arrlw)
            else # if domain.id == "ISS"
                scatter3d!(p, [[co[i]] for i in 1:3]..., color=col, msc=col, msw=0.3, ms=0.8, lw=arrlw)
            end
        end
    end

end