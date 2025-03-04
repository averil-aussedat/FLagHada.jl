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
