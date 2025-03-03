export semiLag
"""
$(SIGNATURES)

Compute the semi-Lagrangian with N steps in time over the given space mesh.

"""
function semiLag(domain::Domain, dynamic::MathFunction, terminalcost::MathFunction, T::Float64, N::Int, mesh::Mesh, verbose::Bool=false)
    
    dt = T / N
    verb(verbose, "Beginning Step 1 with dt = $(numtostr(dt,5)): computation of the reachable sets...")
    reachable_sets = [get_mesh_convex(domain, mesh, union([get_exponential(domain, dyn, xx, dt) for dyn in dynamic.call(xx)]), false) for xx in ProgressBar(mesh.points)]
    # reachable_sets = [[search_in_mesh(domain, mesh, get_exponential(domain, dyn, xx, dt)) for dyn in dynamic.call(xx)] for xx in ProgressBar(mesh.points)]

    # println("reachable sets : ")
    # for (ipt,pt) in enumerate(mesh.points)
    #     if pt.comp == 2
    #         println("pt ", pt_to_str(pt), " : ")
    #         println("\texponentials ", [pt_to_str(qt) for qt in union([get_exponential(domain, dyn, mesh.points[ipt], dt) for dyn in dynamic.call(mesh.points[ipt])])])
    #         println("\treachableset ", [pt_to_str(mesh.points[iqt]) for iqt in reachable_sets[ipt]])
    #     end
    # end

    # for (irr, rr) in enumerate(reachable_sets)
    #     if isempty(rr)
    #         pt = mesh.points[irr]
    #         exps = union([get_exponential(domain, dyn, pt, dt) for dyn in dynamic.call(pt)])
    #         println("Point $irr : ", pt_to_str(pt), " yields exp ", exps, " and convex ", get_mesh_convex(domain, mesh, exps, false))
    #     end
    # end

    # test with optimal dynamic 
    # reachable_sets = [[search_in_mesh(domain, mesh, EPoint([exp(-2*dt*u) * xx.c11, xx.c12, exp(2*dt*u)*xx.c22]))] for xx in mesh.points for u in [sign(xx.alpha - xx.gamma)]]

    mincard = minimum([length(rs) for rs in reachable_sets])
    moycard = mean([length(rs) for rs in reachable_sets])
    maxcard = maximum([length(rs) for rs in reachable_sets])
    verb(verbose, "end of Step 1.")
    verb(verbose, "Extremal cardinals in reachable sets : ", mincard, ", and ", maxcard, ", mean ", moycard, ", for ", mesh.npoints, " mesh points")

    # randpoint = 400
    # counters = [0,0,0]

    jplast = mesh.npoints

    verb(verbose, "Beginning Step 2: propagation...")
    UU = zeros(mesh.npoints, N+1)
    @views UU[:,N+1] .= terminalcost.call.(mesh.points)
    for i in ProgressBar(N:-1:1)
        for jp in 1:mesh.npoints
            # println("reachable set : ", reachable_sets[jp])
            UU[jp,i] = minimum(@views UU[reachable_sets[jp],i+1])

            # ctr = argmin(UU[reachable_sets[jp],i+1])
            # counters[ctr] += 1

            # if jp == randpoint
            #     println("Point $randpoint, i=$i/$N, value ", UU[jp,i], " minimizing ", UU[reachable_sets[jp],i+1])
            # end
            # if jp == jplast
            #     # jplast = reachable_sets[jp][argmin(UU[reachable_sets[jp],i+1])] # new move 
            #     println("choice of the last at step $i / $N : ", pt_to_str(mesh.points[reachable_sets[jp][argmin(UU[reachable_sets[jp],i+1])]]), ", value ", UU[jp,i])
            #     println("reachable set ", [pt_to_str(mesh.points[k]) for k in reachable_sets[jp]])
            # end
        end
    end

    # for i in 1:N
    #     println("Current point ", pt_to_str(mesh.points[jplast]) ,", 
    #         exponentials ", [pt_to_str(pt) for pt in [get_exponential(domain, dyn, mesh.points[jplast], dt) for dyn in dynamic.call(mesh.points[jplast])]], ", 
    #         reachable set ", [pt_to_str(mesh.points[k]) for k in reachable_sets[jplast]], ", 
    #         value ", UU[jplast, i])
    #     jplast = reachable_sets[jplast][argmin(UU[reachable_sets[jplast],i+1])] # new move 
    # end

    # println("Counters : ", counters)

    verb(verbose, "...done.")
    return UU
end


export semiLagInterp
"""
$(SIGNATURES)

Compute the semi-Lagrangian with N steps in time over the given space mesh.

"""
function semiLagInterp(domain::Domain, dynamic::MathFunction, terminalcost::MathFunction, T::Float64, N::Int, mesh::Mesh, verbose::Bool=false)

    dt = T / N
    verb(verbose, "Beginning Step 1: computation of the interpolants...")
    interpolants = [[get_mesh_interpolant(domain, mesh, get_exponential(domain, dyn, xx, dt), false) for dyn in dynamic.call(xx)] for xx in ProgressBar(mesh.points)]
    verb(verbose, "end of Step 1.")

    # println("first interpolant : point ", get_exponential(domain, dynamic.call(mesh.points[1])[1], mesh.points[1], dt))
    # println("interpolant : ", interpolants[1][1])

    verb(verbose, "Beginning Step 2: propagation...")
    UU = zeros(mesh.npoints, N+1)
    @views UU[:,N+1] .= terminalcost.call.(mesh.points)
    for i in ProgressBar(N:-1:1)
        for (jp,interps) in enumerate(interpolants)
            # println("interp : ", interp)
            @views UU[jp,i] = minimum([sum([weight * UU[kp,i+1] for (kp,weight) in interp]) for interp in interps])
        end
    end

    verb(verbose, "...done.")
    return UU
end
    
export get_control
"""
$(SIGNATURES)

Approximate the optimal control. 

"""
function get_control(domain::Domain, dynamic::MathFunction, hatV::Function, T::Float64, N::Int, x::Point, verbose::Bool=false)
    verb(verbose, "Getting optimal control and trajectory from point $(pt_to_str(x))...")
    
    dt = T / N
    haty = Vector{Point}(undef, N+1)
    hatu = Vector{Int64}(undef, N)

    haty[1] = x

    for n = 1:N
        exponentials = [get_exponential(domain, dyn, haty[n], dt) for dyn in dynamic.call(haty[n])]
        
        if n <= N
            hatu[n] = argmin([hatV(n*dt, gf) for gf in exponentials])
        end
        haty[n+1] = exponentials[hatu[n]]
        verb(verbose, "Exponentials for n = $n : ", pt_to_str.(exponentials), ", values ", [hatV(n*dt, gf) for gf in exponentials], ", the control $(hatu[n]) was chosen")
        # verb(verbose, "dt = $dt, dist ", get_distance(domain, haty[n], haty[n+1]))
    end

    verb(verbose, "... done.")
    return haty, hatu
end

function get_feedback_direction(domain::Domain, value::Function, dynamics::Vector{Velocity}, pt::Point, h::Float64, T::Float64)
    accessible_values = [value(T-h, get_exponential(domain, dyn, pt, h)) for dyn in dynamics]
    res = get_exponential(domain, dynamics[argmin(accessible_values)], pt, h)
    return res
end