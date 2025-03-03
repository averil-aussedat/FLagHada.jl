"""
    Cartesian step in (alpha, gamma) with nag divisions 
    Then non-cartesian step in beta in order to keep a mesh size of `step`
"""
function get_mesh(domain::Ellipses, step::Float64, verbose::Bool=false)
    verb(verbose, "Meshing $(domain.id) with step $step...")

    # looking for a vector beta respecting the step constraints 
    gmax = cosh(domain.alga)^2 * cosh(sqrt(2)*domain.beta) - sinh(domain.alga)^2
    max_dist_betas = sqrt(log(gmax + sqrt(gmax^2 - 1.0))^2 + log(gmax - sqrt(gmax^2 - 1.0))^2)
    nb = ceil(Int, max_dist_betas / step)
    betas = collect(LinRange(-domain.beta, domain.beta, nb))

    # Now we know what the total amount of points will be
    if domain.ctedet 
        nag = ceil(Int, 2.0 * sqrt(2) * domain.alga / step)
        npoints = nag * nb 

        # First meshing the cartesian part
        cartpart = [[exp(aa) 0.0; 0.0 exp(-aa)] for aa in LinRange(-domain.alga, domain.alga, nag)]        
    else
        nag = ceil(Int, 2.0 * domain.alga / step)
        npoints = nag^2 * nb 

        # First meshing the cartesian part
        cartpart = [[exp(aa) 0.0; 0.0 exp(gg)] for aa in LinRange(-domain.alga, domain.alga, nag) for gg in LinRange(-domain.alga, domain.alga, nag)]
    end

    # then propagating 
    eV(beta) = [cosh(beta/(2*sqrt(2))) sinh(beta/(2*sqrt(2))); sinh(beta/(2*sqrt(2))) cosh(beta/(2*sqrt(2)))]
    points = [EPoint(eV(beta) * AA * eV(beta)) for beta in betas for AA in cartpart]

    verb(verbose, "... done with $npoints points.")
    return EllMesh(npoints, step, points, nag, nb, betas)
end

"""
    Return the index of a mesh.step-approximation of the closest point
"""
function search_in_mesh(domain::Ellipses, mesh::EllMesh, x::EPoint)

    nalpha = iprojparamFromParam(domain.alga, x.alpha, mesh.nag)
    nbeta = iprojparamFromParam(domain.beta, x.beta, mesh.nb)

    if domain.ctedet 
        return (nbeta-1) * mesh.nag + nalpha 
    else 
        ngamma = iprojparamFromParam(domain.alga, x.gamma, mesh.nag)
        return (nbeta-1) * mesh.nag^2 + (nalpha-1) * mesh.nag + ngamma
    end
end

function imeshFromiABG(ctedet::Bool, ia::Int64, ib::Int64, ig::Int64, nag::Int64)::Int64
    if ctedet # we discard ig
        return 1 + (ib-1)*nag + (ia-1)
    else
        return 1 + (ib-1)*nag^2 + (ia-1)*nag + (ig-1)
    end
end

function paramFromiparam(width::Float64, index::Int64, total::Int64)
    return width*(2.0*(index-1)/(total-1)-1.0)
end

function iprojparamFromParam(width::Float64, param::Float64, total::Int64)
    if width ≈ 0.0
        return 1
    else
        return max(1, min(total, 1 + round(Int, (total-1) * 0.5*(param/width + 1.0))))
    end
end

function otheriproj(width::Float64, param::Float64, total::Int64, iprojparam::Int64)
    if iprojparam <= 1
        return 2
    elseif iprojparam >= total 
        return total-1
    else
        return max(1, min(total, iprojparam + 1 - 2 * signbit(param - width*(2.0*(iprojparam-1)/(total-1)-1.0))))
    end
end

function ibetaFromImesh(ctedet::Bool, imesh::Int64, nag::Int64)
    if ctedet 
        return 1 + floor(Int, (imesh - 1) / nag)
    else
        return 1 + floor(Int, (imesh - 1) / nag^2)
    end
end

function ialphaFromImesh(ctedet::Bool, imesh::Int64, nag::Int64)
    if ctedet 
        return imesh - floor(Int,(imesh-1)/nag)*nag
    else
        return 1 + floor(Int, (imesh - floor(Int,(imesh - 1) / nag^2)*nag^2 - 1) / nag)
    end
end

function igammaFromImesh(imesh::Int64, nag::Int64) # never called if ctedet
    return imesh - floor(Int,(imesh - 1) / nag^2)*nag^2 - floor(Int, (imesh - floor(Int,(imesh - 1) / nag^2)*nag^2 - 1) / nag) * nag
end

# projection of x over [-a,a]
function proj(x::Float64, a::Float64)
    if -a <= x <= a 
        return x 
    else
        return sign(x) * a
    end
end

"""
    Based on a lot of assumptions.
"""
function get_mesh_convex(domain::Ellipses, mesh::EllMesh, points::Vector{EPoint}, verbose::Bool=false)
    verb(verbose, "Getting the mesh convex...")
    
    res = [search_in_mesh(domain, mesh, point) for point in points]
    ibetares = sort(union([ibetaFromImesh(domain.ctedet, re, mesh.nag) for re in res])) # indexes of the betas 
    minibetares = minimum(ibetares)
    # for each betares, what are the indexes of the points belonging to that strata 
    pByBetas = [[ip for ip in 1:length(res) if ibetaFromImesh(domain.ctedet, res[ip],mesh.nag) == betares] for betares in ibetares]
    if domain.ctedet
        # the indexation with beta as the outer index implies 
        # that within each { beta = betares[i] }, 
        # alpha is increasing with ip 
        pByBetas = [union(minimum(bloci),maximum(bloci)) for bloci in pByBetas]
    end
    extremals_by_stratae = [Vector{Int64}() for _ in minibetares:maximum(ibetares)]

    h = 1e-3

    # getting the intersection of the convex with each strata of constant beta 
    for (ibloci, bloci) in enumerate(pByBetas)
        for xi in bloci # index in res
            # adding the current point as an extremal in its strata
            ibxi = ibetaFromImesh(domain.ctedet,res[xi],mesh.nag) # index in mesh.betas
            push!(extremals_by_stratae[ibxi-minibetares+1], res[xi]) # index in mesh.points
            for blocj in pByBetas[ibloci+1:end]
                for xj in blocj 
                    ibxj = ibetaFromImesh(domain.ctedet,res[xj],mesh.nag)
                    # if there is at least one intermediate strata between both 
                    # we look for the intersection of said intermediate strata 
                    # with the geodesic linking points[ix] and points[xj]
                    # 
                    #  pts[ix]
                    #          `⋅.
                    #              `⋅.
                    #                  `⋅.
                    #                  our target
                    #                           `⋅.
                    #                               `⋅.
                    #                                  pts[xj]
                    # 
                    # --|-------- ⋯ ------|-------- ⋯ ---|----> beta
                    # beta[ibxi]      beta[ibxk]      beta[ibxj]
                    #
                    # println("Entering here, ibxi = ", ibxi, ", ibxj = ", ibxj)
                    if ibxi+1 <= ibxj-1
                        vv = Velocity(points[xj], 1.0, 0.0)
                    end
                    for ibxk in ibxi+1:ibxj-1
                        # a coarse Newton method to find a zero of
                        # beta(geod(t)) - beta_obj
                        ceinture = 0; bretelles = 20; goon = true 
                        tn = (mesh.betas[ibxk] - mesh.betas[ibxi]) / (mesh.betas[ibxj] - mesh.betas[ibxi])
                        tnp1 = tn
                        # println("Init : ", tn)

                        while (goon && ceinture <= bretelles)
                            current_beta = get_exponential(domain, vv, points[xi], tn).beta
                            # println("\tpassing")
                            beta_prime = (get_exponential(domain, vv, points[xi], tn+h).beta - current_beta) / h
                            # println("\ttrespassing")
                            tnp1 = tn - proj((current_beta - mesh.betas[ibxk]) / beta_prime, 0.5) # ahem
                            # tnp1 = max(0.0,min(1.0,tn - (current_beta - mesh.betas[ibxk]) / beta_prime))
                            # println("\ttnp1 : ", tnp1, ", beta : ", current_beta, " for target ", mesh.betas[ibxk])
                            if abs(current_beta - mesh.betas[ibxk]) <= 1e-6
                                goon = false 
                            end
                            tn = tnp1
                            ceinture += 1
                        end
                        if ceinture >= bretelles 
                            throw(ErrorException("Newton got out, last tn : $tnp1"))
                        end

                        # println("End\n")
                        # println("Pushing ", search_in_mesh(domain, mesh, get_exponential(domain, vv, points[xi], tnp1)), " in strata $ibxk")

                        push!(extremals_by_stratae[ibxk-minibetares+1], 
                                search_in_mesh(domain, mesh, get_exponential(domain, vv, points[xi], tnp1)))
                    end
                end
            end
        end
    end

    # at this point, extremals_by_stratae contains, for each beta attained by the convex, 
    # a set of points whose convex hull is (?) the intersection of the convex with the mesh and the strata 
    # only stays to do a convex hull in the cartesian part :D

    if domain.ctedet 
        # cartesian part is isometric to a line 
        for bloci in extremals_by_stratae
            append!(res, minimum(bloci):maximum(bloci))
        end
        res = union(res)
    else
        # cartesian part is isometric to a 2D plane
        for (ibeta, extremals) in zip(minibetares:maximum(ibetares), extremals_by_stratae)
            extremals = union(extremals)
            if length(extremals) <= 1
                res = union(res, extremals)
            else 
                indexes = [[ialphaFromImesh(domain.ctedet, extrem, mesh.nag), igammaFromImesh(extrem, mesh.nag)] for extrem in extremals]
                mina = minimum([ii[1] for ii in indexes])
                ming = minimum([ii[2] for ii in indexes])
                maxa = maximum([ii[1] for ii in indexes])
                maxg = maximum([ii[2] for ii in indexes])
                hull = convex_hull(indexes)
                res = union(res, [
                        imeshFromiABG(domain.ctedet,ia,ibeta,ig,mesh.nag) 
                            for ia in mina:maxa 
                                for ig in ming:maxg 
                                    if get_distance_to_2dhull(hull, [ia,ig]) <= 0.5 # convex in integers
                    ]
                )
            end
        end
    end

    verb(verbose, "... done.")
    return res
end

export get_mesh_interpolant 
"""
    Return a list of tuples (index, weight)
    where point ≈ sum weight * mesh.points[index].
    
"""
function get_mesh_interpolant(domain::Ellipses, mesh::EllMesh, point::Point, verbose::Bool=false)

    if domain.ctedet 
        throw(ErrorException("get_mesh_interpolant pas encore implémentééé"))
    end

    verb(verbose, "Interpolating point (abg) = ", point.alpha, ", ", point.beta, ", ", point.gamma)

    if domain.alga ≈ 0.0

        if domain.beta ≈ 0.0 # single point domain. Better treat all the cases...
            return [(1,1.0)]
        else 
            nbeta = iprojparamFromParam(domain.beta, point.beta, mesh.nb)
            mbeta = otheriproj(domain.beta, point.beta, mesh.nb, nbeta)
            distn = get_distance(domain, point, mesh.points[nbeta]) # using the fact that indexation 
            distm = get_distance(domain, point, mesh.points[mbeta]) # can only be 1:nbeta in this case
            return [(nbeta, distm/(distm+distn)), (mbeta, distn/(distm+distn))]

        end

    else
        nalpha = iprojparamFromParam(domain.alga, point.alpha, mesh.nag) 
        ngamma = iprojparamFromParam(domain.alga, point.gamma, mesh.nag)
        malpha = otheriproj(domain.alga, point.alpha, mesh.nag, nalpha)
        mgamma = otheriproj(domain.alga, point.gamma, mesh.nag, ngamma)

        verb(verbose, "nalpha : ", nalpha, ", malpha : ", malpha, ", ngamma : ", ngamma, ", mgamma : ", mgamma)

        if domain.beta ≈ 0.0
            # flat euclidian space 

            dista0 = get_distance(domain, point, EPoint(paramFromiparam(domain.alga,nalpha,mesh.nag),point.beta,point.gamma)) # moving alpha, one direction
            dista1 = get_distance(domain, point, EPoint(paramFromiparam(domain.alga,malpha,mesh.nag),point.beta,point.gamma)) # moving alpha, other direction
            distg0 = get_distance(domain, point, EPoint(point.alpha,point.beta,paramFromiparam(domain.alga,ngamma,mesh.nag))) # moving gamma, one direction
            distg1 = get_distance(domain, point, EPoint(point.alpha,point.beta,paramFromiparam(domain.alga,mgamma,mesh.nag))) # moving gamma, other direction

            verb(verbose, "distances : ", dista0, ", ", dista1, ", ", distg0, ", ", distg1)

            res = [(imeshFromiABG(domain.ctedet, ia,1,ig,mesh.nag), numiaig/((dista0+dista1)*(distg0+distg1))) for (ia,ig,numiaig) in [
                    (nalpha,ngamma,dista1*distg1), 
                    (nalpha,mgamma,dista1*distg0), 
                    (malpha,ngamma,dista0*distg1), 
                    (malpha,mgamma,dista0*distg0)
                ]
            ]
            verb(verbose, "res : ", res)

            return [(imeshFromiABG(domain.ctedet, ia,1,ig,mesh.nag), numiaig/((dista0+dista1)*(distg0+distg1))) for (ia,ig,numiaig) in [
                    (nalpha,ngamma,dista1*distg1), 
                    (nalpha,mgamma,dista1*distg0), 
                    (malpha,ngamma,dista0*distg1), 
                    (malpha,mgamma,dista0*distg0)
                ]
            ]

        else
            # full 3D space, assuming Euclidian for testing 

            nbeta = iprojparamFromParam(domain.beta, point.beta, mesh.nb)
            mbeta = otheriproj(domain.beta, point.beta, mesh.nb, nbeta)
            
            dista0 = get_distance(domain, point, EPoint(paramFromiparam(domain.alga,nalpha,mesh.nag),point.beta,point.gamma)) # moving alpha, one direction
            dista1 = get_distance(domain, point, EPoint(paramFromiparam(domain.alga,malpha,mesh.nag),point.beta,point.gamma)) # moving alpha, other direction
            distb0 = get_distance(domain, point, EPoint(point.alpha,paramFromiparam(domain.beta,nbeta,mesh.nb), point.gamma)) # moving beta, one direction
            distb1 = get_distance(domain, point, EPoint(point.alpha,paramFromiparam(domain.beta,mbeta,mesh.nb), point.gamma)) # moving beta, other direction
            distg0 = get_distance(domain, point, EPoint(point.alpha,point.beta,paramFromiparam(domain.alga,ngamma,mesh.nag))) # moving gamma, one direction
            distg1 = get_distance(domain, point, EPoint(point.alpha,point.beta,paramFromiparam(domain.alga,ngamma,mesh.nag))) # moving gamma, other direction

            return [(imeshFromiABG(domain.ctedet, ia,ib,ig,mesh.nag), numiaibig/((dista0+dista1)*(distb0+distb1)*(distg0+distg1))) for (ia,ib,ig,numiaibig) in [
                    (nalpha,nbeta,ngamma,dista1*distb1*distg1), 
                    (nalpha,nbeta,mgamma,dista1*distb1*distg0), 
                    (nalpha,mbeta,ngamma,dista1*distb0*distg1), 
                    (nalpha,mbeta,mgamma,dista1*distb0*distg0), 
                    (malpha,nbeta,ngamma,dista0*distb1*distg1), 
                    (malpha,nbeta,mgamma,dista0*distb1*distg0), 
                    (malpha,mbeta,ngamma,dista0*distb0*distg1), 
                    (malpha,mbeta,mgamma,dista0*distb0*distg0)
                ]
            ]

        end

    end
end