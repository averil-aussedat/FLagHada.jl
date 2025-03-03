abstract type EuPoint <: Point end

export RdPoint
@with_kw struct RdPoint <: EuPoint
    "Coordinates in ``\\mathbb{R}^d``"
    p::Vector{Float64}
end
Base.:(==)(p1::RdPoint, p2::RdPoint) = p1.p == p2.p

export Domain_euclidean
@with_kw struct Domain_euclidean <: Domain
    "Unique identifier"
    id::String 
    "Dimension of the surrounding euclidean space"
    dim::Int 
    "Bounds of the cube"
    bounds::Vector{RdPoint}
end
Base.:(==)(d1::Domain_euclidean, d2::Domain_euclidean) = (d1.id == d2.id) && (d1.dim == d2.dim) && (d1.bounds == d2.bounds)

# export Euclidean_convex_2D
# abstract type Euclidean_convex_2D <: Convex end

# export Euclidean_convex_1D
# @with_kw struct Euclidean_convex_1D <: Convex 
#     min::Float64
#     max::Float64
# end
# Base.:(==)(c1::Euclidean_convex_1D, c2::Euclidean_convex_1D) = (c1.min == c2.min) && (c1.max == c2.max)

# export Euclidean_convex_2D_singlepoint
# @with_kw struct Euclidean_convex_2D_singlepoint <: Euclidean_convex_2D
#     point::RdPoint
# end
# Base.:(==)(c1::Euclidean_convex_2D_singlepoint, c2::Euclidean_convex_2D_singlepoint) = c1.point == c2.point

# export Euclidean_convex_2D_segment
# @with_kw struct Euclidean_convex_2D_segment <: Euclidean_convex_2D
#     p1::RdPoint
#     p2::RdPoint
# end
# Base.:(==)(c1::Euclidean_convex_2D_segment, c2::Euclidean_convex_2D_segment) = (c1.p1 == c2.p1) && (c1.p2 == c2.p2)

# export Euclidean_convex_2D_fullset
# @with_kw struct Euclidean_convex_2D_fullset <: Euclidean_convex_2D
#     "Number of extremal points"
#     thelen::Int
#     "Extremal points"
#     points::Vector{RdPoint}
#     "Barycenter"
#     bary::RdPoint
#     "Linear transformations of the triangles (bary,p_i,p_{i+1}) towards ((0.0,0.0), (1.0, 0.0), (0.0, 1.0))"
#     matrices::Vector{Matrix{Float64}}
#     vectors::Vector{Vector{Float64}}

#     function Euclidean_convex_2D_fullset(points::Vector{RdPoint})
#         thelen = length(points)
#         bary = RdPoint(sum([p.p for p in points]) ./ length(points))
#         matrices = [zeros(2,2) for _ in 1:thelen]
#         vectors = [zeros(2) for _ in 1:thelen]
#         for (ii,p1) in enumerate(points)
#             p2 = points[ii % thelen + 1]
#             AA = [norm(bary.p .- p1.p)^2 (bary.p .- p1.p)' * (bary.p .- p2.p); (bary.p .- p1.p)' * (bary.p .- p2.p) norm(bary.p .- p2.p)^2]
#             matrices[ii] .= AA \ [(p1.p .- bary.p)'; (p2.p .- bary.p)']
#             vectors[ii] .= matrices[ii] * bary.p
#         end
#         return new(thelen, points, bary, matrices, vectors)
#     end
# end 
# function Base.:(==)(c1::Euclidean_convex_2D_fullset, c2::Euclidean_convex_2D_fullset) 
#     return (c1.thelen == c2.thelen) && (c1.points == c2.points) && (c1.bary == c2.bary) && (c1.matrices == c2.matrices) && (c1.vectors == c2.vectors)
# end

export EuMesh  
@with_kw struct EuMesh <: Mesh
    # common part
    "Number of points"
    npoints::Int
    "min_{points} max_{other points} d(thepoint, theotherpoint)"
    step::Float64
    "List of points"
    points::Vector{RdPoint}
    # specific part
    "Number of layers along each dimension"
    layers::Vector{Int}
end
Base.:(==)(m1::EuMesh, m2::EuMesh) = (m1.npoints == m2.npoints) && (m1.step == m2.step) && (m1.points == m2.points) && (m1.layers == m2.layers)
