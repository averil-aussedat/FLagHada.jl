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

    mincard = minimum([length(rs) for rs in reachable_sets])
    maxcard = maximum([length(rs) for rs in reachable_sets])
    verb(verbose, "end of Step 1.")
    verb(verbose, "Extremal cardinals in reachable sets : ", mincard, ", and ", maxcard, ", for ", mesh.npoints, " mesh points")

    verb(verbose, "Beginning Step 2: propagation...")
    UU = zeros(mesh.npoints, N+1)
    @views UU[:,N+1] .= terminalcost.call.(mesh.points)
    for i in ProgressBar(N:-1:1)
        for jp in 1:mesh.npoints
            UU[jp,i] = minimum(@views UU[reachable_sets[jp],i+1])
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