export get_ellipses
function get_ellipses(id::String, verbose::Bool=false)
    verb(verbose, "Creating $id...")

    if (id == "ellcube")
        domain = Ellipses(id, 1.0, 1.0, false)

    elseif (id == "ellplane")
        domain = Ellipses(id, 1.0, 0.0, false)

    elseif (id == "ellctedet")
        domain = Ellipses(id, 1.5, sqrt(2)*1.5, true) 

    else
        throw(ErrorException("unknown domain id '$id'"))
    end

    verb(verbose, "... done.")
    return domain
end

export get_Scalar
function get_Scalar(id::String, domain::Ellipses, verbose::Bool=false)
    verb(verbose, "Creating scalar function $id...")
    
    if (id == "one")  

        res = MathFunction{Scalar}(id, x::EPoint -> 1.0)
    
    elseif (id == "norm")

        res = MathFunction{Scalar}(id, x::EPoint -> get_distance(domain, EPoint(0.0,0.0,0.0), x))

    elseif (id == "distToSquare")

        res = MathFunction{Scalar}(id, x::EPoint -> get_distance(domain, x, EPoint(x.c11^2, x.c12*(x.c11+x.c22), x.c22^2)))

    elseif (id == "eigmax")

        res = MathFunction{Scalar}(id, x::EPoint -> min(2.0, (x.c11+x.c22)/2 + sqrt(abs((x.c11+x.c22)^2/4 - x.det))))

    else
        throw(ErrorException("unknown scalar id '$id'"))
    end

    verb(verbose, "... done.")
    return res
end


export get_Dynamic
function get_Dynamic(id::String, domain::Ellipses, verbose::Bool=false)
    verb(verbose, "Creating dynamic $id...")

    if (id == "eikonal")

        res = MathFunction{Dynamic}(id, 
            x::EPoint -> [
                Velocity(EPoint(x.alpha+2*domain.alga,x.beta,x.gamma), 1.0, 1.0),
                Velocity(EPoint(x.alpha-2*domain.alga,x.beta,x.gamma), 1.0, 1.0),
                Velocity(EPoint(x.alpha,x.beta,x.gamma+2*domain.alga), 1.0, 1.0),
                Velocity(EPoint(x.alpha,x.beta,x.gamma-2*domain.alga), 1.0, 1.0),
                # 
                Velocity(EPoint(x.alpha+2*domain.alga,x.beta,x.gamma+2*domain.alga), 1.0, 1.0),
                Velocity(EPoint(x.alpha+2*domain.alga,x.beta,x.gamma-2*domain.alga), 1.0, 1.0),
                Velocity(EPoint(x.alpha-2*domain.alga,x.beta,x.gamma+2*domain.alga), 1.0, 1.0),
                Velocity(EPoint(x.alpha-2*domain.alga,x.beta,x.gamma-2*domain.alga), 1.0, 1.0),
                # 
            ]
        )

    elseif (id == "robust")
        """
            Comes from the control problem over the trajectories
                ̇xₜ = A(u(t)) xₜ, 
            where u(⋅) ∈ L⁰(0,T;U), U = {-1,0,1} ⊂ ℝ², and 
                A : u ↦ u R - (1-|u|) I,
            R being the generator of the rotation of angle π/2 and I the identity. 
            In words, 
                if u = -1, we turn clockwise, 
                if u =  0, we go towards 0, 
                if u =  1, we turn anticlockwise.

            The dynamic over the set of matrices is described as 
                dG/dt 
                = dyn_{u(t)}(G)
                = (Gₜ)^{-1/2} A(u(t))ᵗ Gₜ^{1/2} + Gₜ^{1/2} A(u(t)) Gₜ^{-1/2}, 
            with the same controls. 

            Note: if B = [a b; b c], and R = [0 1; -1 0], then
            B R B^{-1} + B^{-1} Rᵗ B = 1/det(B) * [-2b(a+c) a²-c²; a²-c² 2b(a+c)].
            Moreover, B^2 = [a^2+b^2  b(a+c); b(a+c)  b^2+c^2]. Hence if B^2 = X,
            there holds 
                X^{1/2} R X^{-1/2} + X^{-1/2} Rᵗ X^{1/2}
                = 1/sqrt(det(X)) * [-2x.c12  x.c11-x.c22; x.c11-x.c22 2x.c12].

        """

        function robustcall(x::EPoint)
            sqrtx = sqrt([x.c11 x.c12; x.c12 x.c22])
            normR = sqrt(2.0 * ((x.c11+x.c22)^2 / x.det - 4.0)) # trace of the square of the rotating velocity
            return [
                Velocity(EPoint(sqrtx * exp(1.0 / sqrt(x.det) .* [-2*x.c12  x.c11-x.c22; x.c11-x.c22 2*x.c12]) * sqrtx), normR, 1.0), 
                Velocity(EPoint([exp(-2.0)*x.c11, exp(-2.0)*x.c12, exp(-2.0)*x.c22]), 2*sqrt(2.0), 1.0), # square root of the trace  
                Velocity(EPoint(sqrtx * exp(1.0 / sqrt(x.det) .* [2*x.c12  x.c22-x.c11; x.c22-x.c11 -2*x.c12]) * sqrtx), normR, 1.0)
            ]
        end

        res = MathFunction{Dynamic}(id, x::EPoint -> robustcall(x))

    elseif (id == "robustHam")
        """
            Same control problem as in the above, but now
            A has a null trace. We take the famous case 
            U = {-1,1}, and 
            A(u) = [u 0; 0 -u].
        """
        function robustHamcall(x::EPoint)
            # very basic implementation
            sqrtx = sqrt([x.c11 x.c12; x.c12 x.c22])
            isqrtx = inv(sqrtx)
            Vis = [
                sqrtx * [u 0.0; 0.0 -u] * isqrtx + isqrtx * [u 0.0; 0.0 -u] * sqrtx for u in [-1.0, 1.0]
            ]
            normVis = [
                sqrt(Vi[1,1]^2 + 2*Vi[1,2]^2 + Vi[2,2]^2) for Vi in Vis
            ]
            return [
                Velocity(EPoint(sqrtx * exp(2.0 * Vi / (normVi + 1.0*(normVi <= 1e-7))) * sqrtx), normVi, 1.0) for (Vi,normVi) in zip(Vis,normVis)
            ]
        end

        res = MathFunction{Dynamic}(id, robustHamcall)

    elseif (id == "robustRot")
        """
            Same control problem as in the above, but now
            A has a null trace. We take the famous case 
            U = {-1,0,1}, and 
            A(u) = [0 1; u 0].

            Note: if B = [a b; b c], then 
            B A(u) B^{-1} + B^{-1} A(u)ᵗ B
            = 1/det B [2b(uc-a)  uc²-(1+u)b²+a²; uc²-(1+u)b²+a² 2b(a-uc)].
        """
        function robustRotcall(x::EPoint)
            sqrtx = sqrt([x.c11 x.c12; x.c12 x.c22])

            # very basic 
            isqrtx = inv(sqrtx)
            Vis = [
                sqrtx * [0.0 1.0; u 0.0] * isqrtx + isqrtx * [0.0 u; 1.0 0.0] * sqrtx for u in [-1.0, 0.0, 1.0]
            ]
            return [ 
                Velocity(EPoint(sqrtx * exp(2.0 * Vi) * sqrtx), normVi, 1.0) for (Vi,normVi) in zip(Vis,normVis)
            ]
        end

        res = MathFunction{Dynamic}(id, robustRotcall)

    elseif (id == "geodesic")

        function call_geodesic(x::EPoint)
            coeff1 = cosh(5.0 * sqrt((x.c11+x.c22)^2/4 - x.det))
            coeff2 = sinh(5.0 * sqrt((x.c11+x.c22)^2/4 - x.det)) / sqrt((x.c11+x.c22)^2/4 - x.det)
            exp5AA = EPoint(exp(5.0 * (x.c11+x.c22)/2) .* [
                (coeff1 - coeff2*(x.c11+x.c22)/2) + coeff2 * x.c11, # coordinate [1,1] of the exp
                coeff2 * x.c12,                                     # coordinate [1,2] of the exp
                (coeff1 - coeff2*(x.c11+x.c22)/2) + coeff2 * x.c22  # coordinate [2,2] of the exp
            ])
            return [Velocity(exp5AA, t, 1.0) for t in LinRange(-2.0,2.0,7)]
        end
        res = MathFunction{Dynamic}(id, call_geodesic)

    else
        throw(ErrorException("unknown dynamic id '$id'"))
    end

    verb(verbose, "... done.")
    return res
end

export get_Analytical
function get_Analytical(domain::Ellipses, dynamicid::String, costid::String, T::Float64, probe::Bool=false, verbose::Bool=false)
    verb(verbose, "Computing analytical V... ")

    id = domain.id * "_" * dynamicid * "_" * costid

    if (dynamicid == "eikonal" && costid == "norm") 

        if probe 
            return true 
        else
            res = MathFunction{Scalar}(id, (t::Float64, x::EPoint) -> max(0.0, get_distance(domain, EPoint(0.0,0.0,0.0), x) - (T-t)))
        end

    elseif (dynamicid == "robust" && costid == "eigmax") 

        if probe 
            return true 
        else
            res = MathFunction{Scalar}(id, (t::Float64, x::EPoint) -> exp(-2*(T-t)) * ((x.c11+x.c22)/2 + sqrt((x.c11+x.c22)^2/4 - x.det)))
        end

    elseif (dynamicid == "robustHam" && costid == "eigmax") 

        if probe 
            return true 
        else
            function rHamCall(t::Float64, x::EPoint)
                tswitch = min(T-t, 0.25 * log(max(x.c11,x.c22) / min(x.c11,x.c22)))
                trace = exp(-2.0*tswitch) * max(x.c11,x.c22) + exp(2.0*tswitch) * min(x.c11,x.c22)
                return min(2.0, trace/2.0 + sqrt(abs((trace^2)/4.0 - x.det)))
            end

            res = MathFunction{Scalar}(id, rHamCall)
        end
        

    else
        if probe 
            return false 
        else
            throw(ErrorException("unknown analytical solution for domain $(domain.id), dynamic $dynamicid and cost $costid"))
        end
    end

    verb(verbose, "... done.")
    return res 
end

export get_analytical_yu
function get_analytical_yu(domain::Ellipses, dynamic::MathFunction{Dynamic}, costid::String, T::Float64, N::Int64, xx::EPoint, verbose::Bool=false)
    verb(verbose, "Computing analytical (y,u)... ")

    if (dynamic.id == "eikonal" && costid == "norm")

        y = Vector{EPoint}(undef, N+1)
        y[end] = EPoint(0.0,0.0,0.0)
        vv = Velocity(xx, 1.0, 1.0)
        for n in N:-1:1
            y[n] = get_exponential(domain, vv, y[n+1], T/N)
        end

        u = Vector{Velocity}(undef, N)
        for n in 1:N
            dyns = dynamic.call(y[n])
            iu = argmin([get_distance(domain, y[n+1], get_exponential(domain, dyn, y[n], T/N)) for dyn in dyns])
            u[n] = dyns[iu]
        end
        
    elseif (dynamic.id == "robust" && costid == "eigmax") 

        y = [EPoint(exp(-2*t) .* [xx.c11, xx.c12, xx.c22]) for t in LinRange(0.0,T,N+1)]
        u = [Velocity(EPoint([exp(-2.0)*x.c11, exp(-2.0)*x.c12, exp(-2.0)*x.c22]), 1.0, 1.0) for x in y[1:end-1]]

    else
        throw(ErrorException("unknown analytical solution for domain $(domain.id), dynamic $dynamicid and cost $costid"))
    end

    verb(verbose, "... done.")
    return y, u
end