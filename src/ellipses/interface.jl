export pt_to_str
function pt_to_str(pt::EPoint)
    return "alpha $(numtostr(pt.alpha,5)), beta $(numtostr(pt.beta,5)), gamma $(numtostr(pt.gamma,5))"
end

export plot_EPoint
function plot_EPoint(domain::Ellipses, p::EPoint, verbose::Bool=false)
    verb(verbose, "Plotting EPoint...")

    NN = 200
    eigs = eigen([p.c11 p.c12; p.c12 p.c22])
    xx = [eigs.values[1] * cos(2*pi*th) * eigs.vectors[1,1] + eigs.values[2] * sin(2*pi*th) * eigs.vectors[1,2] for th in LinRange(0.0, 1.0, NN)]
    yy = [eigs.values[1] * cos(2*pi*th) * eigs.vectors[2,1] + eigs.values[2] * sin(2*pi*th) * eigs.vectors[2,2] for th in LinRange(0.0, 1.0, NN)]
    p = plot(xx, yy, legend=false, aspect_ratio=:equal)

    verb(verbose, "... done.")
    return p
end

export plot_EPoints
function plot_EPoints(domain::Ellipses, pp::Vector{EPoint}, verbose::Bool=false)
    verb(verbose, "Plotting EPoints...")

    NN = 200
    q = plot(legend=false, aspect_ratio=:equal)
    cols = range(colorant"red", stop=colorant"blue", length=length(pp))
    for (col,p) in zip(cols,pp)
        eigs = eigen([p.c11 p.c12; p.c12 p.c22])
        xx = [eigs.values[1] * cos(2*pi*th) * eigs.vectors[1,1] + eigs.values[2] * sin(2*pi*th) * eigs.vectors[1,2] for th in LinRange(0.0, 1.0, NN)]
        yy = [eigs.values[1] * cos(2*pi*th) * eigs.vectors[2,1] + eigs.values[2] * sin(2*pi*th) * eigs.vectors[2,2] for th in LinRange(0.0, 1.0, NN)]
        plot!(q, xx, yy, color=col)
    end 

    verb(verbose, "... done.")
    return q
end

thearg(a,b,g) = cosh(0.5*(a-g))^2 * cosh(b/sqrt(2)) - sinh(0.5*(a-g))^2
function distabg(a,b,g) 
    arg = thearg(a,b,g)
    return sqrt(log(arg + sqrt(abs(arg^2 - 1.0)))^2 + log(arg - sqrt(abs(arg^2 - 1.0)))^2)
end

export plot_mesh
function plot_mesh(domain::Ellipses, mesh::EllMesh, verbose::Bool=false)
    if mesh.npoints <= 1001

        cols = range(colorant"blue", stop=colorant"red", length=mesh.nb)
        algae = LinRange(-domain.alga, domain.alga, mesh.nag)

        if domain.ctedet 

            p = plot(legend=false, xlabel="α + γ = 0") #, camera=(-45,0))

            for (beta,col) in zip(mesh.betas,cols)
                plane_b = [sign(beta) * distabg(aa,beta,-aa) for aa in algae]
                scatter!(p, algae, plane_b, msc=col, color=col)
            end

        else
            
            plane_a = [aa for aa in algae for gg in algae]
            plane_g = [gg for aa in algae for gg in algae]

            p = plot(legend=false, xlabel="a", ylabel="g") #, camera=(-45,0))

            for (beta,col) in zip(mesh.betas,cols)
                plane_b = [sign(beta) * distabg(aa,beta,gg) for aa in algae for gg in algae]
                scatter3d!(p, plane_a, plane_g, plane_b, msc=col, color=col)
                surface!(p, plane_a, plane_g, plane_b, color=col, alpha=0.2)
            end

        end

    else
        throw(ErrorException("will not plot $(mesh.npoints) points"))
    end

    return p
end

export plot_convex_in_mesh
function plot_convex_in_mesh(domain::Ellipses, mesh::EllMesh, convex::Vector{Int64}, verbose::Bool=false)

    p = plot_mesh(domain, mesh, verbose)
    aas = [mesh.points[c].alpha for c in convex]
    ggs = [mesh.points[c].gamma for c in convex]
    bbs = [sign(mesh.points[c].beta) * distabg(mesh.points[c].alpha, mesh.points[c].beta, mesh.points[c].gamma) for c in convex]

    if domain.ctedet 
        scatter!(p, aas, bbs, msc=:black, color=:black)
    else
        scatter3d!(p, aas, ggs, bbs, msc=:black, color=:black)
    end

    return p
end

export plot_function_plane
"""
    Plotting the values `values` over the 2-dim grid of constant beta = mesh.betas[ibeta]
"""
function plot_function_plane(domain::Ellipses, mesh::EllMesh, values::Vector{Float64}, ibeta::Int64, verbose::Bool=false)

    p = plot(xlims=[-domain.alga,domain.alga], ylims=[-domain.alga,domain.alga])
    
    algae = LinRange(-domain.alga, domain.alga, mesh.nag)
    aas = [aa for aa in algae for gg in algae]
    ggs = [gg for aa in algae for gg in algae]
    # scatter3d!(p, aas, ggs, values)
    surface!(p, aas, ggs, values)

    return p
end

export plot_function_ctedet
"""
    Plotting the values `values` over the 2-dim grid of constant determinant
"""
function plot_function_ctedet(domain::Ellipses, mesh::EllMesh, values::Vector{Float64}, verbose::Bool=false)

    if mesh.npoints <= 1000

        xag = [aa for bb in 1:mesh.nb for aa in LinRange(-domain.alga, domain.alga, mesh.nag)]
        ybb = [sign(bb)*get_distance(domain, EPoint(aa,bb,-aa), EPoint(aa,0.0,-aa)) for bb in LinRange(-domain.beta, domain.beta, mesh.nb) for aa in LinRange(-domain.alga, domain.alga, mesh.nag)]
        p = plot(xlims=[-domain.alga,domain.alga], ylims=[minimum(ybb),maximum(ybb)])
        surface!(p, xag, ybb, values)
        # scatter3d!(p, xag, ybb, values)

    else

        Nnewmesh = 100

        indices_ag = 1:round(Int,mesh.nag/Nnewmesh):mesh.nag
        indices_bb = 1:round(Int,mesh.nb/Nnewmesh):mesh.nb
        verb(verbose, "Selecting a submesh of ", length(indices_ag) * length(indices_bb), " points")

        xag = [aa for bb in indices_bb for aa in LinRange(-domain.alga, domain.alga, mesh.nag)[indices_ag]]
        ybb = [sign(bb)*get_distance(domain, EPoint(aa,bb,-aa), EPoint(aa,0.0,-aa)) for bb in LinRange(-domain.beta, domain.beta, mesh.nb)[indices_bb] 
                for aa in LinRange(-domain.alga, domain.alga, mesh.nag)[indices_ag]]
        p = plot(xlims=[-domain.alga,domain.alga], ylims=[minimum(ybb),maximum(ybb)])
        surface!(p, xag, ybb, values[[mesh.nag * (ibb - 1) + iag for ibb in indices_bb for iag in indices_ag]])
        # scatter3d!(p, xag, ybb, values[[mesh.nag * (ibb - 1) + iag for ibb in indices_bb for iag in indices_ag]],ms=1)

    end

    return p
end

export plot_feedback_map_ctedet
"""
    Plotting an approximation of the feedback over the 2-dim grid of constant determinant.

    The feedback issued from a point is computed as 
    ```math
        \\argmin_{u \\in U} V(0, \\vartheta_{u}(h, G)).
    ```
"""
function plot_feedback_map_ctedet(domain::Ellipses, dynamic::MathFunction{Dynamic}, value::Function, h::Float64, T::Float64, verbose::Bool=false)

    # subgrid of Ngrid^2 in the variables (α, β)
    Ngrid = 20

    # subsampling 
    aas = LinRange(-domain.alga, domain.alga, Ngrid)
    bbs = LinRange(-domain.beta, domain.beta, Ngrid)
    xag = [aa for bb in bbs for aa in aas]
    ybb = [sign(bb)*get_distance(domain, EPoint(aa,bb,-aa), EPoint(aa,0.0,-aa)) for bb in bbs for aa in aas]

    # getting the optimal feedback 
    points = [EPoint(aa,bb,-aa) for bb in bbs for aa in aas]
    inners = [(x.c11 + x.c22)/2.0 + sqrt((x.c11 + x.c22)^2/4 - x.det) < 2.0 for x in points]

    # feedback = [argmin([value(0.0, get_exponential(domain, vv, pt, h)) for vv in dynamic.call(pt)]) for pt in points]
    # truefeed = [1 + (pt.alpha < pt.gamma) for pt in points]
    # directions = [get_exponential(domain, dynamic.call(pt)[ipt], pt, h) for (ipt,pt) in zip(feedback,points)]
    # truedirect = [get_exponential(domain, dynamic.call(pt)[ipt], pt, h) for (ipt,pt) in zip(truefeed,points)]

    directions = [get_feedback_direction(domain, value, dynamic.call(pt), pt, h, T) for pt in points]
    truedirect = [get_exponential(domain, dynamic.call(pt)[1 + (pt.alpha < pt.gamma)], pt, h) for pt in points]

    zz = [[dir.alpha - pt.alpha, dir.beta  - pt.beta] for (dir,pt) in zip(directions, points)]
    # zbb = [  for (dir,pt) in zip(directions, points)]
    tz = [[dir.alpha - pt.alpha, dir.beta  - pt.beta] for (dir,pt) in zip(truedirect, points)]
    # tzb = [dir.beta  - pt.beta  for (dir,pt) in zip(truedirect, points)]

    for z in zz 
        if norm(z) > 1e-7
            z ./= norm(z)
        end
    end
    for z in tz
        if norm(z) > 1e-7
            z ./= norm(z)
        end
    end
    zz .*= 0.45 * (2.0*domain.alga)/Ngrid
    tz .*= 0.45 * (2.0*domain.alga)/Ngrid

    # fact = min(1.0, 0.8 * (2.0*domain.alga) / (Ngrid * maximum(sqrt.(zag.^2 .+ zbb.^2))))
    # verb(verbose, "Multiplication by a factor of ", fact)
    # zag .*= fact
    # zbb .*= fact
    # tza .*= fact
    # tzb .*= fact

    # plotting 
    p = plot(xlabel="α = - γ", ylabel = "β", dpi=300, legend=false)
    function get_col(inn)
        if inn 
            return :black 
        else 
            return :grey
        end
    end
    scatter!(p, xag, ybb, ms=2, color=get_col.(inners), msc=get_col.(inners))
    for (inner, x,y,v0,w0) in zip(inners,xag, ybb,[z[1] for z in tz],[z[2] for z in tz])
        if inner 
            col = :red
            arrow0!(p, x, y, v0, w0, lw=2, as=0.3, lc=col, ec=2.0)
        # else 
        #     col = :salmon
        end
    end
    for (inner, x,y,v1,w1) in zip(inners,xag,ybb,[z[1] for z in zz],[z[2] for z in zz])
        if inner 
            col = :black
        else 
            col = :grey 
        end
        arrow0!(p, x, y, v1, w1, lw=2, as=0.3, lc=col, ec=2.0)
    end
    plot!(p, xlims=[-domain.alga*(1.0+3/Ngrid),domain.alga*(1.0+3/Ngrid)], ylims=[-maximum(ybb)*(1.0+2/Ngrid),maximum(ybb)*(1.0+2/Ngrid)])
    # quiver!(p, xag, ybb, quiver=(tza, tzb), lw=2, color=:red)
    # quiver!(p, xag, ybb, quiver=(zag, zbb), lw=2, color=:black)

    return p
end