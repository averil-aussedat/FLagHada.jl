export NetPoint
abstract type NetPoint <: Point end

export JuncP
@with_kw struct JuncP <: NetPoint 
    "Unique identifier"
    id::Int
    "Identifiers of the connected junction points"
    neighbours::Vector{Int}
end
Base.:(==)(p1::JuncP, p2::JuncP) = (p1.id == p2.id) && (p1.neighbours == p2.neighbours)

export EdgeP
@with_kw struct EdgeP <: NetPoint 
    "Identifiers of the beginning and end of the edge"
    left::Int
    right::Int
    "Distances to the beginning and end of the edge"
    dist_left::Float64
    dist_right::Float64
end
Base.:(==)(p1::EdgeP, p2::EdgeP) = (p1.left == p2.left) && (p1.right == p2.right) && (p1.dist_left == p2.dist_left) && (p1.dist_right == p2.dist_right)

export Network
@with_kw struct Network <: Domain 
    "Unique identifier"
    id::String 
    "Number of junctions"
    njunctions::Int
    "Collection of junction points"
    juncpoints::Vector{JuncP}
    "``\\sup_{(x,y) \\in D^2} d(x,y)``"
    diameter::Float64
    "Half of the smallest loop length"
    CAT0radius::Float64
    "Matrix of distances between junction points"
    distances::Matrix{Float64}
end
Base.:(==)(d1::Network, d2::Network) = (d1.id == d2.id) && (d1.njunctions == d2.njunctions) && (d1.juncpoints == d2.juncpoints) && (d1.CAT0radius == d2.CAT0radius) && (d1.distances == d2.distances)

export Network_plotter
"""
$(TYPEDEF)

A structure to map a network into a plottable (2D) euclidian space.
"""
@with_kw struct Network_plotter 
    "Coordinates of the junction points"
    junccoords::Vector{Vector{Float64}}
end
Base.:(==)(n1::Network_plotter, n2::Network_plotter) = n1.junccoords == n2.junccoords


export NetMesh
@with_kw struct NetMesh <: Mesh
    # common part
    "Number of points"
    npoints::Int64
    "min_{points} max_{other points} d(thepoint, theotherpoint)"
    step::Float64
    "List of points"
    points::Vector{NetPoint}
    # specific part
    "Indices of the junctions in points"
    juncind::Vector{Int64}
    "(Start, stop) indices of the (linear) mesh of an edge"
    edgeind::Matrix{Tuple{Int64,Int64}}
end
Base.:(==)(m1::NetMesh, m2::NetMesh) = (m1.npoints == m2.npoints) && (m1.step == m2.step) && (m1.points == m2.points) && (m1.juncind == m2.juncind) && (m1.edgeind == m2.edgeind)