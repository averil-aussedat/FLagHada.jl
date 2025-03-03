export get_Dynamic
"""
$(SIGNATURES)

Return the function corresponding to id. 

The return set V = [v1, ... ,vn] is supposed to be ordered in the following sense:
the convex hull of V is equal to the polygon enclosed by the path
v1 -> v2 -> ... -> vn -> v1.
"""
function get_Dynamic(id::String, domain::Domain_euclidean, verbose::Bool=false)
    verb(verbose, "Creating dynamic $id...")

    if (id == "translation")

        res = MathFunction{Dynamic}(id, x::RdPoint -> [Velocity(RdPoint(x.p .+ 1.0), 1.0, 1.0)])

    else
        throw(ErrorException("unknown dynamic id '$id'"))
    end

    verb(verbose, "... done.")
    return res
end

export get_Scalar
function get_Scalar(id::String, domain::Domain_euclidean, verbose::Bool=false)
    verb(verbose, "Creating scalar function $id...")

    if (id == "norm")

        res = MathFunction{Scalar}(id, x::RdPoint -> norm(x.p))

    elseif (id == "norm_truncated")

        res = MathFunction{Scalar}(id, x::RdPoint -> min.(norm(x.p),1.0))

    elseif (id == "norm_to_center")

        thecenter = RdPoint(0.5 .* (domain.bounds[1].p .+ domain.bounds[2].p))
        res = MathFunction{Scalar}(id, x::RdPoint -> get_distance(domain, x, thecenter))

    elseif (id == "sinus")

        mul = 8.0
        res = MathFunction{Scalar}(id, x::RdPoint -> sum(sin.(mul .* x.p)) / domain.dim)

    else
        throw(ErrorException("unknown scalar id '$id'"))
    end

    verb(verbose, "... done.")
    return res
end