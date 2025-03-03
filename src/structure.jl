using Parameters # pretty-print composite structures if their declaration is decorated with @with_kw

export Point
abstract type Point end

export Velocity
@with_kw struct Velocity 
    "Point towards which we are moving"
    target::Point 
    "Speed of the motion"
    scale::Float64
    "Treshold of κ"
    thresh::Float64
end
Base.:(==)(v1::Velocity, v2::Velocity) = (v1.target == v2.target) && (v1.scale == v2.scale) && (v1.thresh == v2.thresh)

export Domain
abstract type Domain end

export Convex
abstract type Convex end 

export FuncClass
export Dynamic
export Scalar
abstract type FuncClass end
struct Dynamic <: FuncClass end
struct Scalar <: FuncClass end

export MathFunction
@with_kw struct MathFunction{C<:FuncClass}
    "Unique identifier"
    id::String
    "get the values. Input : Point. Output : list of velocities."
    call::Function
end 

export Mesh
abstract type Mesh end # must all contain the same common part

export OutOfDomainException
struct OutOfDomainException <: Exception 
end