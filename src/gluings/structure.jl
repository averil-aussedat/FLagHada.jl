export GPoint 
@with_kw struct GPoint <: EuPoint
    "Component in the gluing"
    comp::Int64 
    "Euclidean coordinates within this component"
    p::Vector{Float64}
end
Base.:(==)(p1::GPoint, p2::GPoint) = ((p1.comp == p2.comp) && (p1.p ≈ p2.p))

export Gluing
@with_kw struct Gluing <: Domain 
    "Unique identifier"
    id::String 
    "Number of components"
    ncomps::Int64
    "Components"
    comps::Vector{Domain_euclidean}
    """
    Crossing map:
        To x,y in different components, associates the (possibly empty) sequence (zᵢ)ᵢ 
        of points in the geodesic [xy] jumping from a comp. to another.

        Each point zᵢ is duplicated in its representation in both components.

        If called with (x,x), associates the sequence of points identified with x, i.e. 
        [x] if x belongs to the interior of a component, or [x,y,z]... if d(x,y) = d(x,z) = 0.
    """
    get_crossing::Function # signature: (p1::GPoint, p2::GPoint)::Float64

    Gluing(id::String, comps::Vector{Domain_euclidean}, get_crossing::Function) = new(id, length(comps), comps, get_crossing)
end
Base.:(==)(d1::Gluing, d2::Gluing) = ((d1.id == d2.id) && (d1.comps == d2.comps) && (d1.get_crossing == d2.get_crossing))

export Gluing_plotter
@with_kw struct Gluing_plotter 
    "dimension of the ambiant Euclidean space"
    dim::Int64 
    "translator of local to plotting coordinates"
    get_plotting_coords::Function # signature: (domain::Gluing, p::GPoint)::Vector{Float64}
end

export GMesh 
@with_kw struct GMesh <: Mesh 
    # common part
    "Number of points"
    npoints::Int
    "min_{points} max_{other points} d(thepoint, theotherpoint)"
    step::Float64
    "List of points"
    points::Vector{GPoint}
    # specific part
    "Index ranges of each component in the mesh"
    compinits::Vector{Int64}
    compstops::Vector{Int64}
    "Layers in each components"
    complayers::Vector{Vector{Int}}
end
Base.:(==)(m1::GMesh, m2::GMesh) = (m1.npoints == m2.npoints) && (m1.step == m2.step) && (m1.points == m2.points) && 
                                   (m1.compinits == m2.compinits) && (m1.compstops == m2.compstops) && (m1.complayers == m2.complayers) 

