export get_distance 
function get_distance(domain::Domain_euclidean, point1::RdPoint, point2::RdPoint)
    return norm(point2.p .- point1.p)
end

# export get_convex
# """
# $(SIGNATURES)

# Compute the convex hull of points. 
# NOT TESTED
# It is assumed that points are all disjoints.

# """
# function get_convex(domain::Domain_euclidean, points::Vector{RdPoint}, trim=true)
#     if (domain.dim == 1)
#         return Euclidean_convex_1D(minimum([p.p[1] for p in points]), maximum([p.p[1] for p in points]))
#     elseif (domain.dim == 2)
#         thelen = length(points)
#         if (thelen == 1)
#             return Euclidean_convex_2D_singlepoint(points[1])
#         elseif (thelen == 2)
#             return Euclidean_convex_2D_segment(points[1], points[2])
#         else 
#             baryp = sum([p.p for p in points]) ./ thelen

#             # order the points monotonically around bary 
#             angles = [atan(p.p[2]-baryp[2],p.p[1]-baryp[1]) for p in points]
#             points .= points[sortperm(angles)]
#             angles = sort(angles)

#             if trim
#                 # remove the points that would have the same angles 
#                 ceinture = 0; bretelles = 10000; ip = 0; last_problem = 0
#                 goon = true 
#                 while (goon && (ceinture <= bretelles))
#                     ip = ip % length(points) + 1 # next point
#                     iq = ip % length(points) + 1 # next next point
#                     if min(abs(angles[ip]-angles[iq]), abs(abs(angles[ip]-angles[iq])-2pi)) <= 1e-6
#                         last_problem = ceinture
#                         if norm(points[ip].p .- baryp) <= norm(points[iq].p .- baryp)
#                             angles = angles[union(1:ip-1,ip+1:end)]
#                             points = points[union(1:ip-1,ip+1:end)]
#                         else
#                             angles = angles[union(1:iq-1,iq+1:end)]
#                             points = points[union(1:iq-1,iq+1:end)]
#                         end 
#                     end
#                     goon = (ceinture - last_problem <= length(points)+1)
#                     ceinture += 1
#                 end

#                 # remove the points in the interior of the convex hull 
#                 # quite naive 
#                 ceinture = 0; bretelles = 10000; ip = 0; last_problem = 0
#                 goon = true 
#                 while (goon && (ceinture <= bretelles))
#                     ip = ip % length(points) + 1 # next point
#                     p = points[ip]
#                     angles = [atan(
#                         # < q - p, R (bary - p) >, where R is the rotation by -pi/2
#                         (q.p[1]-p.p[1])*(baryp[2] - p.p[2]) + (q.p[2]-p.p[2])*(p.p[1] - baryp[1]),
#                         # < q - p, bary - p >
#                         (q.p[1]-p.p[1])*(baryp[1] - p.p[1]) + (q.p[2]-p.p[2])*(baryp[2] - p.p[2])
#                         ) for q in @views circshift(points, -ip)[1:end-1]]
#                     # if the points are all extremal points of the hull, 
#                     # then angles should be decreasing. 
#                     # so we remove the ones before the *last* maximum
#                     agmax = length(angles) - argmax(@views angles[end:-1:1]) + 1
#                     # record the fact that we had to remove points 
#                     if agmax > 1
#                         last_problem = ceinture 
#                     end
#                     points = points[union(1:ip, (ip+agmax):end)]
#                     baryp = sum([p.p for p in points]) ./ length(points)
#                     # if we still not have made one full turn without removing points
#                     goon = (ceinture - last_problem <= length(points)+1)
#                     ceinture += 1
#                 end
#             end # trim

#             return Euclidean_convex_2D_fullset(points)
#         end
#     else 
#         fname = nameof(var"#self#")
#         throw(ErrorException("dimension $(domain.dim) not implemented in $fname"))
#     end
# end

# function dist_to_2d_segment(p1::RdPoint, p2::RdPoint, yj::RdPoint)
#     hh = max(0.0,min(1.0, (yj.p .- p1.p)' * (p2.p .- p1.p) / norm(p2.p .- p1.p)^2))
#     return norm(yj.p .- p1.p .+ hh .* (p1.p .- p2.p))
# end

# export get_dist_to_convex
# function get_dist_to_convex(domain::Domain_euclidean, convex::Euclidean_convex_1D, yj::RdPoint)
#     return max(0.0, convex.min - yj.p[1], yj.p[1] - convex.max)
# end

# function get_dist_to_convex(domain::Domain_euclidean, convex::Euclidean_convex_2D_singlepoint, yj::RdPoint)
#     return norm(convex.point.p .- yj.p)
# end

# function get_dist_to_convex(domain::Domain_euclidean, convex::Euclidean_convex_2D_segment, yj::RdPoint)
#     return dist_to_2d_segment(convex.p1, convex.p2, yj)
# end

# """
# $(SIGNATURES)

# Return the distance to a convex set in 2D. 

# Computed as the min of the distances between yj and the triangles (bary,p_i,p_{i+1}), 
# where (p_i)_i are the extremal points, ordered to that the piecewise linear curve 
# gluing the [p_i, p_{i+1}] turns around the barycenter monotonically. 

# """
# function get_dist_to_convex(domain::Domain_euclidean, convex::Euclidean_convex_2D_fullset, yj::RdPoint)
#     res = Inf 
#     for (ii,p1) in enumerate(convex.points)
#         p2 = convex.points[ii % convex.thelen + 1]
#         # first check if yj is in the triangle (bary, p_i, p_{i+1})
#         ll = convex.matrices[ii] * yj.p .- convex.vectors[ii]
#         if ((ll[1] >= 0.0) && (ll[2] >= 0.0) && (ll[1]+ll[2] <= 1.0))
#             return 0.0
#         # else, the distance is attained on the boundary, hence compute the min 
#         # between the distances to the edges 
#         else
#             res = min(res,
#                 dist_to_2d_segment(p1, p2, yj),
#                 dist_to_2d_segment(p1, convex.bary, yj),
#                 dist_to_2d_segment(p2, convex.bary, yj)
#             )
#         end
#     end
#     return res
# end

"""
Return the point ``exp(h \\cdot v)`` for ``v`` the given velocity.
"""
function get_exponential(domain::Domain_euclidean, velocity::Velocity, p::RdPoint, h::Float64)
    if norm(p.p .- velocities.target) <= 1e-8
        return p
    else
        h = get_kappa_time(velocity, h, norm(velocity.target.p .- p.p))
        return RdPoint(p.p .+ h .* velocity.scale .* (velocity.target.p .- p.p) ./ norm(velocity.target.p .- p.p))
    end
end