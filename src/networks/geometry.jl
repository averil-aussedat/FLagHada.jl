export get_distance
function get_distance(domain::Network, point1::JuncP, point2::JuncP)
    return domain.distances[point1.id, point2.id]
end

function get_distance(domain::Network, point1::JuncP, point2::EdgeP)
    return min(domain.distances[point1.id, point2.left] + point2.dist_left, domain.distances[point1.id, point2.right] + point2.dist_right)
end

function get_distance(domain::Network, point1::EdgeP, point2::JuncP)
    return min(domain.distances[point2.id, point1.left] + point1.dist_left, domain.distances[point2.id, point1.right] + point1.dist_right)
end

function get_distance(domain::Network, point1::EdgeP, point2::EdgeP)
    if ((point1.left == point2.left) && (point1.right == point2.right))
        return abs(point1.dist_left - point2.dist_left)
    elseif ((point1.left == point2.right) && (point1.right == point2.left))
        return abs(point1.dist_left - point2.dist_right)
    else
        return min(
            domain.distances[point1.left,  point2.left]  + point1.dist_left  + point2.dist_left ,
            domain.distances[point1.right, point2.left]  + point1.dist_right + point2.dist_left ,
            domain.distances[point1.left,  point2.right] + point1.dist_left  + point2.dist_right,
            domain.distances[point1.right, point2.right] + point1.dist_right + point2.dist_right
        )
    end
end

"""
Return the number of the neighbour of p so that 
the geodesic [pq] points towards that neighbour.
"""
function get_direction(domain::Network, p::JuncP, q::JuncP)
    if q.id in p.neighbours
        return q.id 
    else 
        return p.neighbours[argmin([domain.distances[r, q.id] for r in p.neighbours])]
    end
end

function get_direction(domain::Network, p::JuncP, q::EdgeP)
    if q.left == p.id 
        return q.right 
    elseif q.right == p.id 
        return q.left 
    else 
        return p.neighbours[argmin([domain.distances[r, q.left] for r in p.neighbours])]
    end
end

"""
Return p.left if q is towards the left neighbour, p.right otherwise.
"""
function get_direction(domain::Network, p::EdgeP, q::JuncP)
    if domain.distances[p.left, q.id] < domain.distances[p.right, q.id]
        return p.left 
    else 
        return p.right
    end
end

function get_direction(domain::Network, p::EdgeP, q::EdgeP)
    if p.left == q.left 
        if p.right == q.right 
            # same edge 
            if p.dist_left <= q.dist_left 
                #   p.left --- p --------- p.right
                #   q.left -------- q ---- q.right
                return p.right 
            else 
                #   p.left -------- p ---- p.right
                #   q.left --- q --------- q.right
                return p.left 
            end
        else 
            #                   p.left ---- p ---- p.right 
            # q.right --- q --- q.left 
            return p.left 
        end 
    elseif p.right == q.right
        #                  p.right ---- p ---- p.left
        # q.left --- q --- q.right  
        return p.right
    elseif p.left == q.right 
        if p.right == q.left 
            # same edge, inverted order 
            if p.dist_left <= q.dist_right
                #   p.left  --- p --------- p.right
                #   q.right -------- q ---- q.left 
                return p.right 
            else 
                #   p.left  -------- p ---- p.right
                #   q.right --- q --------- q.left
                return p.left 
            end 
        else 
            #                  p.left  ---- p ---- p.right 
            # q.left --- q --- q.right 
            return p.left
        end 
    elseif p.right == q.left 
        #                   p.right ---- p ---- p.left
        # q.right --- q --- q.left 
        return p.right 
    else # no neighbour in common
        if domain.distances[p.right,q.left] < domain.distances[p.left,q.left]
            #   p.left ----- p ----- p.right
            #                           |
            #   q.right ----- q ----- q.left
            return p.right 
        else 
            #   p.left ----- p ----- p.right
            #      | 
            #   q.right ----- q ----- q.left
            return p.left
        end
    end 
end

export get_convex
"""
Return the extremal points of the convex hull of "points". 
It is assumed that points are all disjoints.
"""
function get_convex(domain::Network, points::Vector{<:NetPoint}, trim=true)
    # We only have to trim the points that lie in the interior of the hull. 
    # Which means, in 1D CAT(0) networks, that for any point p in Extremals, 
    # all geodesics linking p to other points in Extremals are going in the 
    # same direction. 

    if trim
        ceinture = 0; bretelles = 10000; ip = 0; last_problem = 0
        goon = true 
        while (goon && (ceinture <= bretelles))
            ip = ip % length(points) + 1 # next point
            p = points[ip]
            # compare all directions to the first one 
            iq = ip % length(points) + 1 # point after p 
            firstdir = get_direction(domain, p, points[iq])
            is_extremal = true; hiq = 1
            while is_extremal && (hiq <= length(points)-2)
                if firstdir == get_direction(domain, p, points[(iq+hiq-1)%length(points)+1])
                    hiq += 1
                else
                    is_extremal = false
                    last_problem = ceinture
                    points = points[union(1:(ip-1),(ip+1):end)]
                end
            end
            # if we still not have made one full turn without removing points
            goon = (ceinture - last_problem < length(points)+1)
            ceinture += 1
        end
    end

    # return Network_convex(points)
    return points
end

# export get_distance_to_convex
# """
# $(SIGNATURES)

# Return the distance to a convex set in 1D CAT(0) networks. 

# """
# function get_dist_to_convex(domain::Network, convex::Network_convex, yj::NetPoint)
#     # Test if yj lies within the convex hull 
#     firstdir = get_direction(domain, yj, convex.points[1])
#     goon = true; iq = 2
#     while goon && (iq <= convex.thelen)
#         if firstdir == get_direction(domain, yj, convex.points[iq])
#             iq += 1
#         else
#             goon = false
#             return 0.0
#         end
#     end
#     return minimum([get_distance(domain, p, yj) for p in convex.points])
# end

export get_exponential
"""
Return the point ``exp(h \\cdot v)`` for ``v`` the given velocity.
Based on the idea that velocity is pointing towards a boundary point. 
If we would have to get after the boundary point, returns an error. 
"""
function get_exponential(domain::Network, velocity::Velocity, p::JuncP, h::Float64; recur::Bool=false)
    if !recur
        h = get_kappa_time(velocity, h, get_distance(domain, velocity.target, p))
    end
    if (h * velocity.scale > get_distance(domain, p, velocity.target))
        # throw(OutOfDomainException())
        return velocity.target
    elseif (velocity.scale <= 1e-8) || ((typeof(velocity.target)==JuncP) && (velocity.target.id == p.id))
        return p
    else
        dir = get_direction(domain, p, velocity.target)
        # if we don't have to change edge  
        if h * velocity.scale <= domain.distances[p.id, dir]
            if abs(h * velocity.scale - domain.distances[p.id, dir]) <= 1e-8
                return domain.juncpoints[dir]
            else 
                return EdgeP(p.id, dir, h * velocity.scale, domain.distances[p.id,dir] - h * velocity.scale)
            end
        # if we have to change edge, recursive
        else
            return get_exponential(domain, velocity, domain.juncpoints[dir], h - domain.distances[p.id,dir]/velocity.scale, recur=true)
        end
    end
end

function get_exponential(domain::Network, velocity::Velocity, p::EdgeP, h::Float64; recur::Bool=false)
    if !recur
        # print("changing h = ", h)
        h = get_kappa_time(velocity, h, get_distance(domain, velocity.target, p))
        # println(" into ", h)
    end
    if (h * velocity.scale > get_distance(domain, p, velocity.target))
        # throw(OutOfDomainException())
        return velocity.target
    elseif (velocity.scale <= 1e-8) || (get_distance(domain,p,velocity.target) <= 1e-8)
        return p
    else
        # if we don't have to change edge  
        dir = get_direction(domain, p, velocity.target)
        dist_dir = p.dist_left * (dir == p.left) + p.dist_right * (dir == p.right)
        if h * velocity.scale <= dist_dir
            if abs(h * velocity.scale - dist_dir) <= 1e-8
                return domain.juncpoints[dir]
            else
                if dir == p.left 
                    return EdgeP(dir, p.right, p.dist_left - h * velocity.scale, p.dist_right + h * velocity.scale)
                else 
                    return EdgeP(dir, p.left, p.dist_right - h * velocity.scale, p.dist_left + h * velocity.scale)
                end
            end
        # if we have to change edge, recursive
        else 
            return get_exponential(domain, velocity, domain.juncpoints[dir], h - dist_dir/velocity.scale, recur=true)
        end
    end
end