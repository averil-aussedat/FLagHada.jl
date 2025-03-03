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
        # println("bd : ", bd)

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
    # println("vertices : ", vertices)
    # println("edges : ", edges)

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
            # println("plotting edge ", edge)
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
        # if domain.id == "planeline"
        #     subnum = 40
        #     # subnum = sqrt(mesh.compstops[1]-mesh.compinits[1]) / 50
        #     nxy = mesh.complayers[1][1]
        #     step = ceil(Int,nxy/subnum)
        #     ind1 = [1 + (i-1) * nxy + j-1 for i in collect(1:step:nxy) ∪ [nxy] for j in collect(1:step:nxy) ∪ [nxy]]
        #     # ind1 = collect(mesh.compinits[1]:ceil(Int,sqrt(mesh.compstops[1]-mesh.compinits[1])/subnum):mesh.compstops[1]) ∪ [mesh.compstops[1]]
        #     # println("ind1 : ", ind1)
        #     ind2 = collect(mesh.compinits[2]:ceil(Int,(mesh.compstops[2]-mesh.compinits[2])/subnum):mesh.compstops[2]) ∪ [mesh.compstops[2]]
        #     surface!(p, [tp[1] for tp in toplot[ind1]], [tp[2] for tp in toplot[ind1]], values[ind1])
        #     # surface!(p, [tp[1] for tp in toplot[ind2]], [tp[2] for tp in toplot[ind2]], values[ind2])
        #     scatter3d!(p, [tp[1] for tp in toplot[ind2]], [tp[2] for tp in toplot[ind2]], values[ind2], ms=ms, mc=colors[ind2], msc=colors[ind2])
        #     # plot!(p, [tp[1] for tp in toplot[ind2]], [tp[2] for tp in toplot[ind2]], values[ind2], ms=3.0)
        # else
            scatter!(p, [tp[1] for tp in toplot], [tp[2] for tp in toplot], mc=colors, msc=colors, ms=ms, shape=:rect)
        # end
    else
        if domain.id == "ISS"
            toplot .= vcat(toplot[mesh.compstops[end]:-1:mesh.compinits[end]], toplot[begin:mesh.compinits[end]-1])
            colors .= vcat(colors[mesh.compstops[end]:-1:mesh.compinits[end]], colors[begin:mesh.compinits[end]-1])
        end

        scatter3d!(p, [tp[1] for tp in toplot], [tp[2] for tp in toplot], [tp[3] for tp in toplot], mc=colors, msc=colors, ms=ms)
        plot!(p, title = "Value between $(numtostr(minv,6)) and $(numtostr(maxv,6))")
        verb(verbose, "value between $minv and $maxv")
    end
    # plot!(p, colorbar=true)
end

export plot_feedback_map!
function plot_feedback_map!(p, domain::Gluing, plotter::Gluing_plotter, dynamic::MathFunction{Dynamic}, value, h::Float64, T::Float64, verbose::Bool=false)

    # if domain.id == "planeline"
    #     # subgrid of Ngrid^2 
    #     Ngrid = 15

    #     # getting the optimal feedback 
    #     lrange = LinRange(0.05,0.95,Ngrid)
    #     points1 = [GPoint(1,[xx,yy]) for xx in lrange for yy in lrange]
    #     points2 = [GPoint(2,[rr]) for rr in lrange]

    #     directions1 = [get_feedback_direction(domain, value, dynamic.call(pt), pt, h, T) for pt in points1]
    #     directions2 = [get_feedback_direction(domain, value, dynamic.call(pt), pt, h, T) for pt in points2]
    #     # truedirect = [get_exponential(domain, dynamic.call(pt)[1 + (pt.alpha < pt.gamma)], pt, h) for pt in points]

    #     zz1 = [dir.p .- pt.p for (dir,pt) in zip(directions1, points1)]
    #     zz2 = [sign(dir.p[1] .- pt.p[1]) for (dir,pt) in zip(directions2, points2)]
    #     # tz = [[dir.alpha - pt.alpha, dir.beta  - pt.beta] for (dir,pt) in zip(truedirect, points)]
    #     # tzb = [dir.beta  - pt.beta  for (dir,pt) in zip(truedirect, points)]

    #     for z in zz1 
    #         if norm(z) > 1e-7
    #             z ./= norm(z)
    #         end
    #     end
    #     # for z in tz
    #     #     if norm(z) > 1e-7
    #     #         z ./= norm(z)
    #     #     end
    #     # end
    #     zz1 .*= 0.3 * 2.0/Ngrid
    #     zz2 .*= 0.3 * 2.0/Ngrid
    #     # tz .*= 0.45 * 2.0/Ngrid

    #     # plotting 
    #     plot!(p, xlabel="x", ylabel = "y", dpi=300, legend=false, xlims=[-0.1,2.1], ylims=[-0.1,1.1])
    #     # function get_col(inn)
    #     #     if inn 
    #     #         return :black 
    #     #     else 
    #     #         return :grey
    #     #     end
    #     # end
    #     i1 = 1; i2 = 1
    #     for xx in lrange
    #         for yy in lrange
    #             arrow0!(p, xx, yy, zz1[i1][1], zz1[i1][2], lw=1.5, as=0.3, lc=:black, ec=2.0)
    #             i1 += 1
    #         end
    #         arrow0!(p, 1.0 + xx, 0.5, zz2[i2], 0.0, lw=2, as=0.3, lc=:black, ec=2.0)
    #         i2 += 1
    #     end

    # else

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
    # nextcrossd = [domain.get_crossing(pt, dir)[1] for (pt,dir) in zip(mesh.points, directions)]
    coords = [plotter.get_plotting_coords(domain, pt.comp, pt) for pt in hmesh.points]
    arrows = [plotter.get_plotting_coords(domain, ncd.comp, ncd) .- co for (co, ncd) in zip(coords, nextcrossd)]
    for arr in arrows 
        if norm(arr) > 1e-6
            arr .*= larr / norm(arr)
        end
    end
    # println("coords : ", coords)
    # println("arrows : ", arrows)

    col = :orangered
    for (co,arr) in zip(coords, arrows)
        # scatter3d!([[co[i]] for i in 1:3]..., color=:brown, ms=0.5)
        arrow0!(p, co..., arr...; lc=col, ec=2.0, as=0.3, lw=arrlw)
        if norm(arr) <= 1e-5
            scatter3d!(p, [[co[i]] for i in 1:3]..., color=col, msc=col, msw=0.3, ms=0.8, lw=arrlw)
        end
    end

end