export get_distance 
function get_distance(domain::Gluing, p::GPoint, q::GPoint)
    if p.comp == q.comp
        return norm(p.p .- q.p)
    else
        crossings = domain.get_crossing(p,q)
        return norm(p.p .- crossings[1].p) + sum([norm(crossings[i].p .- crossings[i+1].p) for i in 2:2:(length(crossings)-2)]) + norm(crossings[end].p .- q.p)
    end
end

export get_exponential

function get_exponential(domain::Gluing, vv::Velocity, p::GPoint, h::Float64)
    h = get_kappa_time(vv, h, get_distance(domain, vv.target, p))
    if get_distance(domain, vv.target, p) <= 1e-8 || vv.scale <= 1e-8
        return p 
    elseif vv.target.comp == p.comp 
        return GPoint(p.comp, max.(
            domain.comps[p.comp].bounds[1].p, min.(
                domain.comps[p.comp].bounds[2].p,
                p.p .+ min(1.0, h * vv.scale ./ norm(vv.target.p .- p.p)) .* (vv.target.p .- p.p)
            )))
    else 
        # not the same component, so crossing has at least 2 elements
        crossings = domain.get_crossing(p,vv.target)

        if h * vv.scale <= norm(p.p .- crossings[1].p) # we stay within this component 
            return GPoint(p.comp, max.(
                domain.comps[p.comp].bounds[1].p, min.(
                    domain.comps[p.comp].bounds[2].p,
                    p.p .+ h * vv.scale .* (crossings[1].p .- p.p) ./ norm(crossings[1].p .- p.p)
                )))
        else 
            crossed_distance = norm(p.p .- crossings[1].p)
            goon = true; ceinture = 1; bretelles = domain.ncomps + 2
            while goon && (ceinture <= bretelles)
                if 2*ceinture == length(crossings) # last crossing before entering the component of the target 
                    new_distance = norm(crossings[2*ceinture].p .- vv.target.p)
                else
                    new_distance = norm(crossings[2*ceinture].p .- crossings[2*ceinture+1].p)
                end
                if h * vv.scale > crossed_distance + new_distance 
                    crossed_distance += new_distance
                    ceinture += 1
                else
                    goon = false 
                end
            end
            if 2*ceinture == length(crossings) # we got to the last component
                dir = vv.target 
            else # somewhere intermediate
                dir = crossings[2*ceinture+1]
            end
            return GPoint(crossings[2*ceinture].comp, max.(
                domain.comps[crossings[2*ceinture].comp].bounds[1].p, min.(
                    domain.comps[crossings[2*ceinture].comp].bounds[2].p,
                    crossings[2*ceinture].p .+ (h * vv.scale - crossed_distance) .* (dir.p .- crossings[2*ceinture].p) ./ norm(dir.p .- crossings[2*ceinture].p)
                )))
        end
    end
end
