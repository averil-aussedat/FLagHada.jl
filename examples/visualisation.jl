function create_and_save_pics(testcase::TestCase, testcaseid::String, foldername::String, domain::Domain, plotter, dynamic::MathFunction{Dynamic}, 
                        mesh::Mesh, hatVV::Matrix{Float64}, hatV, pictype, iN::Int64, N::Int64, dx::Float64, verbose::Bool)

    if testcase.domaintype == "network"
        if (pictype == "value" || (typeof(pictype) == Vector{String} && "value" in pictype))
            if mesh.npoints <= 1000
                submesh = collect(1:mesh.npoints)
            else
                submesh = mesh.juncind ∪ collect(1:round(Int,mesh.npoints/1000):mesh.npoints) 
                verb(verbose, "Plotting on a submesh of ", length(submesh), " points")
            end
            p = plot_function(domain, plotter, mesh, hatVV[:,1], submesh, verbose)
            if testcaseid == "eikonale_tripod"
                zlims = [-0.1,1.1]
            elseif testcaseid == "sncf_intricate"
                zlims = [-0.1,6.1]
            else 
                zlims = [-0.1,1.1]
            end
            plot!(p, zlims=zlims, dpi=300)
            savefig(p, "data/$foldername/$(testcase.runid)_V$(iN)"); verb(verbose, "Saving $(testcase.runid)_V$(iN)")
        end
        if (pictype == "feedback" || (typeof(pictype) == Vector{String} && "feedback" in pictype))
            q = plot_feedback_map(domain, plotter, dynamic, hatV, testcase.T/N, testcase.T, verbose)
            display(q)
            savefig(q, "data/$foldername/$(testcase.runid)_FB$(iN)"); verb(verbose, "Saving $(testcase.runid)_FB$(iN)")
        end

    elseif testcase.domaintype == "ellipses"
        if (pictype == "value" || (typeof(pictype) == Vector{String} && "value" in pictype))
            p = plot_function_ctedet(domain, mesh, hatVV[:,1], verbose)
            plot!(p, zlims=[0.9,2.01], xlabel="α = -γ", ylabel="β", dpi=300, camera=(50,20), colorbar=false, grid=true)
            savefig(p, "data/$foldername/$(testcase.runid)_V$(iN)"); verb(verbose, "Saving $(testcase.runid)_V$(iN)")
        end
        if (pictype == "feedback" || (typeof(pictype) == Vector{String} && "feedback" in pictype))
            q = plot_feedback_map_ctedet(domain, dynamic, hatV, testcase.T/N, testcase.T, verbose)
            display(q)
            savefig(q, "data/$foldername/$(testcase.runid)_FB$(iN)"); verb(verbose, "Saving $(testcase.runid)_FB$(iN)")
        end

    elseif testcase.domaintype == "gluing"
        if (pictype == "value" || (typeof(pictype) == Vector{String} && "value" in pictype))
            verb(verbose, "Plotting value...")
            if domain.id == "planeline"
                p = plot(camera=(-120,30), legend=false, zlim=[-0.05, 1.0], aspect_ratio=1.0)
                ms = max(1.5, dx*36)
            else
                p = plot_gluing(domain, plotter)
                ms = max(0.8,dx*7.5)
            end
            plot_values!(p, domain, plotter, mesh, hatVV[:,1], ms)
            plot!(p, title="", dpi=400)
            savefig(p, "data/$foldername/$(testcase.runid)_V$(iN)"); verb(verbose, "Saving $(testcase.runid)_V$(iN)")
            verb(verbose, "... done.")
        end
        if (pictype == "feedback" || (typeof(pictype) == Vector{String} && "feedback" in pictype))
            q = plot_gluing(domain, plotter)
            plot!(q,dpi=400)
            plot_feedback_map!(q, domain, plotter, dynamic, hatV, testcase.T/N, testcase.T, verbose)
            savefig(q, "data/$foldername/$(testcase.runid)_FB$(iN)"); verb(verbose, "Saving $(testcase.runid)_FB$(iN)")
        end
    else 
        throw(ErrorException("what is domain type $(testcase.domaintype) exactly"))
    end

end

function compute_errors_controls!(testcase::TestCase, domain::Domain, dynamic::MathFunction{Dynamic}, hatV, iN::Int64, N::Int64, 
                                  initpoints::Vector{<:Point}, save_errors::Bool, has_analytical::Bool, verbose::Bool)

    worserr = -1.0; worsrerr = -1.0
    meanerr = 0.0;  meanrerr = 0.0; cardmeanrerr = 0

    for initp in initpoints
        haty, hatu = get_control(domain, dynamic, hatV, testcase.T, N, initp, false)

        if save_errors && has_analytical
            err = abs(terminalcost.call(haty[end]) - analytical.call(0.0, initp))
            worserr = max(worserr, err)
            meanerr += err
            if abs(analytical.call(0.0, initp)) > 1e-7
                rerr = abs(terminalcost.call(haty[end]) - analytical.call(0.0, initp)) / abs(analytical.call(0.0, initp))
                worsrerr = max(worsrerr, rerr)
                meanrerr += rerr
                cardmeanrerr += 1
            else 
                println("Warning: renormalized loss of optimality not computed")
            end
            verb(verbose, "Starting from ", pt_to_str(initp), ", getting to ", pt_to_str(haty[end]), ", err : ", err)
        end
    end
    
    testcase.erroruSL_lossw[iN] = worserr 
    testcase.erroruSL_lossm[iN] = meanerr / length(initpoints)
    testcase.erroruSL_percw[iN] = worsrerr 

    if cardmeanrerr > 0
        testcase.erroruSL_percm[iN] = meanrerr / cardmeanrerr
    end
    verb(verbose, "Loss of optimality : worst ", testcase.erroruSL_lossw[iN], ", mean ", testcase.erroruSL_lossm[iN], 
                      ", renormalized : worst ", testcase.erroruSL_percw[iN], ", mean ", testcase.erroruSL_percm[iN])

end