function get_distance_to_segment(point, aa, bb)
    if ((point .- aa)' * (bb .- aa) <= 0) || ((point .- bb)' * (aa .- bb) <= 0)
        # one of two extremals 
        return min(sqrt((point[1]-aa[1])^2+(point[2]-aa[2])^2), sqrt((point[1]-bb[1])^2+(point[2]-bb[2])^2))
    else
        # projection on the line 
        return sqrt((point[1]-aa[1])^2+(point[2]-aa[2])^2 - ((point.-aa)' * (bb .- aa))^2 / ((aa[1]-bb[1])^2 + (aa[2]-bb[2])^2))
    end
end

function get_distance_to_2dhull(hull, point)
    if point in VPolygon(hull)
        return 0.0 
    else
        return minimum([get_distance_to_segment(point, hull[i], hull[1+(i%length(hull))]) for i in eachindex(hull)])
    end
end

"""
    Return h such that 
        the gradient flow of scale * d(⋅,x₀) issued from x at time h
    coincides with 
        the gradient flow of scale * κ∘d(⋅,x₀) issed from x at time t. 
    Inputs:
        * vv a velocity  
        * t the allowed running time 
        * truedist the distance between x and x₀.

    Here κ : R⁺ → R⁺ is equal to 
        x - thresh/2        if x >= thresh 
        x^2/(2*thresh)      if x < thresh.

    Let y_h = d(x_h,x₀). In the first system, one has 
        d/dh y_h = - scale if y_h > 0, and 0 otherwise. 
    The solution is max(0,d(x,x₀) - scale * h).
    In the second system, one has 
        d/dt y_t = - scale if y_t > thresh, and - scale * y_t / thresh otherwise. 
    The solution is 
        y_t = d(x,x₀) - scale * t if this quantity is >= thresh, 
    and 
        y_t = exp(- scale/thresh * (t - max(0,d(x,x₀)-thresh)/scale)) min(d(x,x₀),thresh) otherwise.
        
    In the case where d(x,x₀) - scale * t >= thresh, both t and h are equal. 
    Otherwise, one wants 
        d(x,x₀) - scale * h = exp(- scale/thresh * (t - max(0,d(x,x₀)-thresh)/scale)) min(d(x,x₀),thresh)
                            = exp(- scale/thresh *  t + max(0,d(x,x₀)/thresh  - 1.0)) min(d(x,x₀),thresh)
"""
function get_kappa_time(vv::Velocity, t::Float64, truedist::Float64)
    if (abs(vv.scale) <= 1e-7) || (abs(vv.thresh) <= 1e-7) || (truedist - vv.scale * t >= vv.thresh)
        return t
    else 
        return (truedist - exp(- vv.scale/vv.thresh * t + max(0.0,truedist/vv.thresh - 1.0))*min(truedist,vv.thresh)) / vv.scale 
    end
end