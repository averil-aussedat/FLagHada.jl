export expV 
expV(beta::Float64) = [cosh(beta/sqrt(2)) sinh(beta/sqrt(2)); sinh(beta/sqrt(2)) cosh(beta/sqrt(2))]

export get_distance 
"""
    Compute the distance 
    ```math
        \\delta_2(A, B) = \\left( \\sum_{j=1}^2 \\log^2(\\lambda_j (A^{-1} B)) \\right)^{1/2}.
    ```
    To compute the eigenvalues of A^{-1} B, we use the dimension 2, in which they are given by
    ```math
        \\frac{\\text{Trace}(C) \\pm \\sqrt{\\text{Trace}^2(C) - 4 \\det C}}{2}.
    ```
    Since the det of `A^{-1} B` is given by `det(B) / det(A)`, it suffices to compute trace(A^{-1} B).
    Again by explicit computations, in the space of symmetric matrices, 
    ```math 
        \\Trace(A^{-1} B) = \\frac{A_{22} B_{11} - 2 A_{12} B_{12} + A_{11} B_{22}}{\\det A}.
    ```
    Plugging this in the above, we get 
    ```math
        \\lambda_{\\pm} = \\frac{S \\pm \\sqrt{S^2 - 4 \\det A \\det B}}{2 \\det A}, 
        \\qquad \\text{where} \\qquad 
        S = A_{22} B_{11} - 2 A_{12} B_{12} + A_{11} B_{22}.
    ```
"""
function get_distance(domain::Ellipses, A1::EPoint, A2::EPoint)::Float64
    if abs(A1.alpha-A2.alpha) + abs(A1.beta-A2.beta) + abs(A1.gamma-A2.gamma) <= 1e-7
        return 0.0
    else
        # optimized expression in the case of positive matrices 
        S = A1.c11 * A2.c22 - 2 * A1.c12 * A2.c12 + A1.c22 * A2.c11 
        return sqrt(log((S - sqrt(abs(S^2 - 4*A1.det*A2.det))) / (2.0*A1.det))^2 + log((S + sqrt(abs(S^2 - 4*A1.det*A2.det))) / (2.0*A1.det))^2)
    end
end

export get_exponential_slow
""" 
    For testing purposes
"""
function get_exponential_slow(domain::Ellipses, velocity::Velocity, p::EPoint, h::Float64)
    dist = get_distance(domain, p, velocity.target)
    h = get_kappa_time(velocity, h, dist)
    if dist <= 1e-8
        return p 
    else
        sqrtE = sqrt([p.c11 p.c12; p.c12 p.c22]) 
        isqrtE = inv(sqrtE)
        return EPoint(sqrtE * (isqrtE * [velocity.target.c11 velocity.target.c12; velocity.target.c12 velocity.target.c22] * isqrtE)^(velocity.scale * h / dist) * sqrtE)
    end
end

export get_exponential 
"""
    Return the evaluation at time ``h`` of geodesic linking ``p`` to ``velocity.target``. 
    This function being critical, it it as optimized as I can.
"""
function get_exponential(domain::Ellipses, velocity::Velocity, p::EPoint, h::Float64)
    dist = get_distance(domain, p, velocity.target)
    # print("Changing h = ", h)
    h = get_kappa_time(velocity, h, dist)
    # println(" into ", h)
    if dist <= 1e-8
        return p 
    else 
        h = velocity.scale * h / dist
        if abs(velocity.target.beta - p.beta) <= 1e-7 # moving along the image of a cartesian plane through the isometry A -> e^{beta/2 V} A e^{beta/2 V}
            return EPoint((1-h) * p.alpha + h * velocity.target.alpha, p.beta, (1-h) * p.gamma + h * velocity.target.gamma)
        else
            ### Denote B = p^{-1/2} * target * p^{-1/2}. Then
            rho0 = exp(0.5*(p.alpha + p.gamma)) * cosh(0.5*(p.alpha - p.gamma))
            rho1 = exp(0.5*(velocity.target.alpha + velocity.target.gamma)) * cosh(0.5*(velocity.target.alpha - velocity.target.gamma))
            tau0 = tanh(0.5*(p.alpha - p.gamma))
            tau1 = tanh(0.5*(velocity.target.alpha - velocity.target.gamma))
            traceB = 2.0 * rho1 / (rho0 * (1.0 - tau0^2)) * (cosh((p.beta - velocity.target.beta)/sqrt(2)) - tau0 * tau1)
            deterB = velocity.target.det / p.det 
            hdiffl = sqrt(traceB^2/4 - deterB) # half of the difference between eigenvalues 

            # the eigenvalues of B are (traceB +- sqrt(traceB^2 - 4 * deterB))/2
            # computing the coefficients kI and kB of the decomposition log(B) = kI * I + kB * B 
            if abs(traceB^2 - 4 * deterB) <= 1e-8 # if both eigenvalues are "the same"
                # then the unique eigenvalue is traceB/2 
                kI = log(traceB/2.0) - 1.0
                kB = 2.0 / traceB
            else 
                lBp = (traceB + sqrt(traceB^2 - 4 * deterB)) / 2.0
                lBm = (traceB - sqrt(traceB^2 - 4 * deterB)) / 2.0
                kI = (log(lBm) * lBp - lBm * log(lBp)) / (lBp - lBm)
                kB = (log(lBp) - log(lBm)) / (lBp - lBm)
            end

            # using now the explicit formula for the exponential 
            coeffp = exp(h * (kI + kB * traceB/2)) * (cosh(h * kB * hdiffl) - traceB/2 * sinh(h * kB * hdiffl) / hdiffl)
            coefft = exp(h * (kI + kB * traceB/2)) *  sinh(h * kB * hdiffl) / hdiffl

            # println("coeffp : ", coeffp, ", coefft : ", coefft)

            return EPoint([
                coeffp * p.c11 + coefft * velocity.target.c11, 
                coeffp * p.c12 + coefft * velocity.target.c12, 
                coeffp * p.c22 + coefft * velocity.target.c22 
            ])
        end
    end
end