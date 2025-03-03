export pt_to_str
function pt_to_str(pt::NetPoint)
    if typeof(pt) == JuncP 
        return "junction " * numtostr(pt.id,1+floor(Int64,pt.id/10))
    else
        dmin = pt.dist_left * (pt.left < pt.right) + pt.dist_right * (pt.left > pt.right)
        dmax = pt.dist_right * (pt.left < pt.right) + pt.dist_left * (pt.left > pt.right)
        jmin = min(pt.left,pt.right)
        jmax = max(pt.left,pt.right)
        return "edge [" * numtostr(jmin,1+floor(Int64,jmin/10)) * " " * numtostr(dmin,5) * " x " * numtostr(dmax,5) * " " * numtostr(jmax,1+floor(Int64,jmax/10)) * "]" 
    end
end

function pt_to_coords(domain::Network, plotter::Network_plotter, pt::Point)
    if typeof(pt) == JuncP
        coords = plotter.junccoords[pt.id]
    else
        α = 1.0 - pt.dist_left / domain.distances[pt.left, pt.right]
        coords = α .* plotter.junccoords[pt.left] .+ (1.0 - α) .* plotter.junccoords[pt.right]
    end
    return coords
end

function plot_edges!(p, bounds::Matrix{Float64}, domain::Network, plotter::Network_plotter, dim3::Bool=true, printid::Bool=true)
    for istart = 1:domain.njunctions
        for istop = domain.juncpoints[istart].neighbours
            xcoords = [plotter.junccoords[istart][1], plotter.junccoords[istop][1]]
            ycoords = [plotter.junccoords[istart][2], plotter.junccoords[istop][2]]
            if dim3
                zcoords = [0.0 for _ in 1:length(xcoords)]
                Plots.plot3d!(p, xcoords, ycoords, zcoords, color=:gray)
            else
                plot!(p, xcoords, ycoords, color=:gray)
            end
            bounds[:,1] .= min.(bounds[:,1], [minimum(xcoords), minimum(ycoords)])
            bounds[:,2] .= max.(bounds[:,2], [maximum(xcoords), maximum(ycoords)])
        end
        if (printid && !dim3)
            annotate!(plotter.junccoords[istart][1]-0.1, plotter.junccoords[istart][2]+0.1, "$istart")
        end
    end
end

function plot_junction_points!(p, bounds::Matrix{Float64}, domain::Network, plotter::Network_plotter, idtoplot::Vector{<:Any}, vals::VecOrMat{Float64}, dim3::Bool=true, ms::Float64=5.0)
    if (length(idtoplot) > 0)
        juncpoints_x = [plotter.junccoords[ii][1] for ii in idtoplot]
        juncpoints_y = [plotter.junccoords[ii][2] for ii in idtoplot]
        if dim3 
            scatter3d!(p, juncpoints_x, juncpoints_y, 0.0 .* vals, color=:gray, ms=2, shape=:diamond)
            scatter3d!(p, juncpoints_x, juncpoints_y, vals, marker_z=vals, ms=2)
        else
            scatter!(p, juncpoints_x, juncpoints_y, marker_z=vals, ms=ms)
        end
        bounds[:,1] .= min.(bounds[:,1], [minimum(juncpoints_x), minimum(juncpoints_y)])
        bounds[:,2] .= max.(bounds[:,2], [maximum(juncpoints_x), maximum(juncpoints_y)])
    end
end

function plot_edge_points!(p, domain::Network, plotter::Network_plotter, mesh::Mesh, subvals::VecOrMat{Float64}, submesh::Vector{Int64}, dim3::Bool=true)
    for (ipt, pt) in enumerate(mesh.points[submesh])
        if (typeof(pt) == EdgeP)
            α = 1.0 - pt.dist_left / domain.distances[pt.left, pt.right]
            coords = α .* plotter.junccoords[pt.left] .+ (1.0 - α) .* plotter.junccoords[pt.right]
            if dim3
                scatter3d!(p, [coords[1]], [coords[2]], [subvals[ipt]], marker_z=[subvals[ipt]], ms=1)
            else
                scatter!(p, [coords[1]], [coords[2]], marker_z=[subvals[ipt]], ms=2, markershape=:square)
            end
        end
    end
end

export plot_function 
function plot_function(domain::Network, plotter::Network_plotter, vals::VecOrMat{Float64}, dim3::Bool=true, verbose::Bool=false, ms::Float64=5.0)
    verb(verbose, "Plotting on $(domain.id) ...")

    printid = false

    p = plot(legend=false, aspect_ratio=:equal, camera=(60,10), axis=false); # colorbar=true, 
    bounds = [-0.5 0.5; -0.5 0.5]
    plot_edges!(p, bounds, domain, plotter, dim3, printid)
    plot_junction_points!(p, bounds, domain, plotter, collect(1:domain.njunctions), vals, dim3, ms)
    bounds .+= 0.2 .* (bounds[:,2] .- bounds[:,1]) .* [-1.0 1.0; -1.0 1.0]
    plot!(p, xlims=bounds[1,:], ylims=bounds[2,:])

    # display(p)
    verb(verbose, "... done.")
    return p
end

function plot_function(domain::Network, plotter::Network_plotter, mesh::Mesh, vals::VecOrMat{Float64}, submesh::Vector{Int64}=[0], dim3::Bool=true, verbose::Bool=false)
    verb(verbose, "Plotting on $(domain.id) ...")

    if submesh == [0]
        submesh = collect(1:mesh.npoints)
        subvals = vals
    else 
        subvals = vals[submesh]
    end
    subvals[subvals .== -Inf] .= -1e-6

    bounds = [-0.5 0.5; -0.5 0.5]
    printid = true

    p = plot(legend=false, aspect_ratio=:equal, camera=(25,30)); # colorbar=true, 
    plot_edges!(p, bounds, domain, plotter, dim3, printid)
    plot_edge_points!(p, domain, plotter, mesh, subvals, submesh, dim3)

    idtoplot   = [pt.id for pt in mesh.points[submesh] if typeof(pt) == JuncP]
    juncvals = [subvals[ipt] for (ipt, pt) in enumerate(mesh.points[submesh]) if typeof(pt) == JuncP]

    plot_junction_points!(p, bounds, domain, plotter, idtoplot, juncvals, dim3)
    bounds .+= 0.2 .* (bounds[:,2] .- bounds[:,1]) .* [-1.0 1.0; -1.0 1.0]
    plot!(p, xlims=bounds[1,:], ylims=bounds[2,:])

    # display(p)

    verb(verbose, "... done.")
    return p
end

export plot_trajectory!
function plot_trajectory!(p, domain::Network, plotter::Network_plotter, haty::Vector{Point}, col, verbose::Bool=false)
    verb(verbose, "Plotting trajectory issued from ", pt_to_str(haty[1]), "...")

    # plotting the network itself
    # p = plot_function(domain, plotter, zeros(domain.njunctions), false, false)

    coords = zeros(2,length(haty))
    for (ipt,pt) in enumerate(haty)
        coords[:,ipt] = pt_to_coords(domain, plotter, pt)
    end

    plot!(p, coords[1,:], coords[2,:], color=col, lw=2)

    verb(verbose, "... done.")
end

export plot_feedback_map
function plot_feedback_map(domain::Network, plotter::Network_plotter, dynamic::MathFunction{Dynamic}, value::Function, h::Float64, T::Float64, verbose::Bool=false)
    verb(verbose, "Plotting feedback map...")

    # space step of the plotting mesh
    if domain.id == "intricate"
        dxplot = 0.65
        coordfact = 0.4
        ec = 2.0
        as = 0.6
        lw = 2.0
    else
        dxplot = 0.2
        coordfact = 0.12
        ec = 1.0
        as = 0.2
        lw = 3.0
    end
    meshplot = get_mesh(domain, dxplot)

    # getting the optimal feedback 
    # feedback = [argmin([value(0.0, get_exponential(domain, vv, pt, h)) for vv in dynamic.call(pt)]) for pt in meshplot.points]
    directions = [domain.juncpoints[get_direction(domain, pt, get_feedback_direction(domain, value, dynamic.call(pt), pt, h, T))] for pt in meshplot.points]
    coordspts = [pt_to_coords(domain, plotter, pt) for pt in meshplot.points]
    coordsdir = [pt_to_coords(domain, plotter, dir) .- pt_to_coords(domain, plotter, pt) for (dir,pt) in zip(directions, meshplot.points)]

    for coorddir in coordsdir
        if norm(coorddir) >= 1e-7
            coorddir .= coordfact .* coorddir ./ norm(coorddir)
        end
    end

    # plotting 
    p = plot_function(domain, plotter, zeros(domain.njunctions), false, false, 3.0)
    plot!(p, dpi=300)

    for (cpt, arr) in zip(coordspts, coordsdir)
        arrow0!(p, cpt[1], cpt[2], arr[1], arr[2], lw=lw, as=as, lc=:blue, ec=ec)
    end

    verb(verbose, "... done.")
    return p
end