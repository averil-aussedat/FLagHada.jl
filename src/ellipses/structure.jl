export EPoint
@with_kw struct EPoint <: Point 
    "Coefficients"
    c11::Float64
    c12::Float64
    c22::Float64
    "Determinant"
    det::Float64
    "Coordinates in the exponential base"
    alpha::Float64
    beta::Float64
    gamma::Float64
    expalpha::Float64 
    coshfrac::Float64 # cosh (b / (2 sqrt(2)))
    sinhfrac::Float64 # sinh (b / (2 sqrt(2)))
    expgamma::Float64

    function EPoint(alpha::Float64, beta::Float64, gamma::Float64)
        coshfrac = cosh(beta / (2*sqrt(2)))
        sinhfrac = sinh(beta / (2*sqrt(2)))
        expalpha = exp(alpha)
        expgamma = exp(gamma)
        c11 = expalpha * coshfrac^2 + expgamma * sinhfrac^2
        c22 = expgamma * coshfrac^2 + expalpha * sinhfrac^2
        c12 = (expalpha + expgamma) * coshfrac * sinhfrac
        new(c11, c12, c22, c11*c22-c12^2, alpha, beta, gamma, expalpha, coshfrac, sinhfrac, expgamma)
    end

    function EPoint(A::Vector{Float64}) 
        beta = sqrt(2.0) * atanh(2.0 * A[2] / (A[1] + A[3]))
        coshfrac = cosh(beta / (2*sqrt(2)))
        sinhfrac = sinh(beta / (2*sqrt(2)))
        expalpha = (A[1] * coshfrac^2 - A[3] * sinhfrac^2) / (coshfrac^4 - sinhfrac^4)
        expgamma = (A[3] * coshfrac^2 - A[1] * sinhfrac^2) / (coshfrac^4 - sinhfrac^4)
        new(A[1], A[2], A[3], A[1]*A[3]-A[2]^2, log(expalpha), beta, log(expgamma), expalpha, coshfrac, sinhfrac, expgamma)
    end

    function EPoint(A::Matrix{Float64}) 
        return EPoint([A[1,1], A[1,2], A[2,2]])
    end
end
Base.:(==)(p1::EPoint, p2::EPoint) = (p1.c11 == p2.c11) && (p1.c12 == p2.c12) && (p1.c22 == p2.c22) && (p1.det == p2.det)
Base.:(+)(p1::EPoint, p2::EPoint) = EPoint([p1.c11+p2.c11, p1.c12+p2.c12, p1.c22+p2.c22])

export Ellipses 
"""
    It is assumed that the domain corresponds to 
    - either a centered parallelepiped in the 3D space 
    of coordinates (alpha,beta,gamma), given by 
    [-alga,alga] x [-beta,beta] x [-alga,alga],
    - either the intersection of the above with the 
    surface alpha + gamma = 0.
    In the first case, ctedet = false, and ctedet = true in the second case.
"""
@with_kw struct Ellipses <: Domain 
    "Unique identifier"
    id::String 
    "Maximal coordinates"
    alga::Float64
    beta::Float64
    "True if constant determinant (equal to 1), false otherwise"
    ctedet::Bool
end

export EllMesh
@with_kw struct EllMesh <: Mesh
    # common part
    "Number of points"
    npoints::Int
    "min_{points} max_{other points} d(thepoint, theotherpoint)"
    step::Float64
    "List of points"
    points::Vector{EPoint}
    # specific part
    "Number of divisions in the coordinates alpha and gamma"
    nag::Int64
    "Number of divisions in the coordinate beta"
    nb::Int64
    "Coordinates beta"
    betas::Vector{Float64}
end
Base.:(==)(m1::EllMesh, m2::EllMesh) = (m1.npoints == m2.npoints) && (m1.step == m2.step) && (m1.points == m2.points) && (m1.nag == m2.nag) && (m1.nb == m2.nb) && (m1.betas == m2.betas)
