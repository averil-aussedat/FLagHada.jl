using LinearAlgebra
using Plots
using Polynomials
using Colors

export verb
function verb(verbose::Bool, theverbs...) # I is proud 
    if verbose
        if string(stacktrace()[2].func) != "top-level scope"
            println("[$(string(stacktrace()[2].func))] ", theverbs...)
        else
            println(theverbs...)
        end
    end
end

export numtostr
function numtostr(r, n::Int64)
    if 'e' in "$r"
        r = 0.0
    end
    thestr = "$r"[1:min(end,n)]
    return thestr * "0"^(n-length(thestr))
end

export plot_errors
"""
$(SIGNATURES)

Simple loglog display of series of errors.
- listeh : vector of steps 
- errors : series to plot. One colums = one series
- names : vector of tags

"""
function plot_errors(listeh::Vector{Float64}, errors::VecOrMat{Float64}, names::Vector{String}, verbose::Bool=false)
    verb(verbose, "Creating error plot... ")

    p = plot(xaxis=:log, yaxis=:log,xlabel="step",ylabel="error",legend=:outerbottom)
    abserr = abs.(errors)
    col = ["#B8255F", "#FF9933", "#AFB83B", "#6ACCBC", "#4073FF", "#AF38EB", "#808080", "#CCAC93"]
    mk = [:circle, :diamond, :xcross, :star5, :rect, :ltriangle, :pentagon]
    abserr[abserr .<= 1e-14] .= 1e-14
    for i=1:length(errors[1,:])
        if (maximum(abserr[:,i]) < 1e-10)
            verb(verbose, "Skipping series $i, identically null")
        else
            poly = fit(log10.(listeh), log10.(abserr[:,i]), 1)
            plot!(p, listeh, abserr[:,i], shape=mk[i], ls=:solid, lc=col[i], mc=col[i], lw=1.0, label=names[i])
            poly_approx = 10 .^(poly.(log10.(listeh)))
            plot!(p, listeh, poly_approx, shape=mk[i], ls=:dash,  lc=col[i], mc=col[i], lw=1.0, label="a = " * "$(poly.coeffs[2])"[1:5])
        end
    end

    verb(verbose, "... done.")
    return p
end

# needs LinRegOutliers
# export plot_errors_robust
# """
# $(SIGNATURES)

# Robust loglog display of series of errors.
# - listeh : vector of steps 
# - errors : series to plot. One colums = one series
# - names : vector of tags

# """
# function plot_errors_robust(listeh::Vector{Float64}, errors::VecOrMat{Float64}, names::Vector{String}, verbose::Bool=false)
#     verb(verbose, "Creating error plot... ")

#     xticks = 10 .^(round.(LinRange(minimum(log10.(listeh)),maximum(log10.(listeh)),5),digits=2))
#     p = plot(xaxis=:log, yaxis=:log, xlabel="step", ylabel="error", legend=:outerbottom, xticks=xticks)
#     abserr = abs.(errors) # 2e, 3e, 4e, 5e, 
#     col1 = ["#6ACCBC", "#4073FF", "#FF9933", "#B8255F", "#D5E04C"]
#     col2 = ["#2B9483", "#4A1167", "#A7611B", "#61092C", "#7F8813"]
#     mk = [:circle, :diamond, :xcross, :star5, :rect, :ltriangle, :pentagon]
#     abserr[abserr .<= 1e-14] .= 1e-14
#     for i=1:length(errors[1,:])
#         if (maximum(abserr[:,i]) < 1e-10)
#             verb(verbose, "Skipping series $i, identically null")
#         else
#             # outliers-robust linear regression. Doc https://github.com/jbytecode/LinRegOutliers/blob/master/examples.md
#             reg = createRegressionSetting(@formula(le ~ lh), DataFrame((lh = log10.(listeh), le=log10.(abserr[:,i]))))
#             thelts = lts(reg)
#             poly = thelts["betas"]
#             poly_approx = 10 .^(poly[2] .* log10.(listeh) .+ poly[1])

#             # if outliers are those who do not fit the window of the regression 
#             # extension = maximum(poly_approx) - minimum(poly_approx)
#             # topline    = maximum(poly_approx) + 0.3 * extension
#             # bottomline = minimum(poly_approx) - 0.3 * extension
#             # inliers = [j for j in 1:length(listeh) if bottomline <= abserr[j,i] <= topline]
#             # outliers= [j for j in 1:length(listeh) if j ∉ inliers]
            
#             outliers = thelts["outliers"]
#             inliers = [j for j in 1:length(listeh) if !(j in outliers)]
#             bottomline = minimum(abserr[inliers,i])

#             # println("inliers : ", inliers)
#             # println("outliers : ", outliers)

#             scatter!(p, listeh[inliers], abserr[inliers,i], shape=mk[i], ls=:solid, lc=col1[i], mc=col1[i], msc=col2[i], lw=1.0, label=names[i], ms=3)
#             plot!(p, listeh[inliers], poly_approx[inliers], ls=:dash, lc=col2[i], mc=col2[i], lw=2.0, label="a = " * "$(poly[2])"[1:min(end,5)]) # shape=mk[i], , ms=2
#             if length(outliers) > 0
#                 verb(verbose, "Plotting ", length(outliers), " outlier" * "s"^(length(outliers) > 1) * " on series $i")
#                 scatter!(p, listeh[outliers], fill(bottomline,length(outliers)), shape=mk[i], mc="red", msc="darkred", lw=1.0, ms=3, label=:none)
#             end
#         end
#     end

#     verb(verbose, "... done.")
#     return p
# end

# inspired from https://discourse.julialang.org/t/quiver-plot-arrowhead-sizing/81789/4
# as multiplies the branches of the arrow head 
# ec larger = pointier arrow 
function arrow0!(p, x, y, u, v; as=0.07, lw=1, lc=:black, la=1, ec=3.0)
    v1, v2 = [u;v], [-v;u]
    v4 = as * (ec*v1 + v2)/sqrt(1.0 + ec^2)
    v5 = as * (ec*v1 - v2)/sqrt(1.0 + ec^2)
    plot!(p, [x,x+u], [y,y+v], lw=lw, lc=lc, la=la)
    plot!(p, [x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lw=lw, lc=lc, la=la)
    plot!(p, [x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lw=lw, lc=lc, la=la)
end

function arrow0!(p, x, y, z, u, v, w; as=0.07, lw=1, lc=:black, la=1, ec=3.0)
    v1, v2 = [u;v], [-v;u]
    v4 = as * (ec*v1 + v2)/sqrt(1.0 + ec^2)
    v5 = as * (ec*v1 - v2)/sqrt(1.0 + ec^2)
    plot!(p, [x,x+u], [y,y+v], [z,z+w], lw=lw, lc=lc, la=la)
    plot!(p, [x+u,x+u-v5[1]], [y+v,y+v-v5[2]], [z+w,z+w], lw=lw, lc=lc, la=la)
    plot!(p, [x+u,x+u-v4[1]], [y+v,y+v-v4[2]], [z+w,z+w], lw=lw, lc=lc, la=la)
end