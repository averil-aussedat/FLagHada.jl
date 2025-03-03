using FLagHada
using Parameters # with_kw
using DocStringExtensions
using ProgressBars
using Plots
using DelimitedFiles

include("testcases.jl")
include("visualisation.jl")

println("Welcome in main.")

#####################################################
# Parameters 
#####################################################

VSL = true 
uSL = true
save_errors = true
save_values = true
create_pics = true

### Networks

# testcaseid = "eikonale_tripod"
# testcaseid = "eikonale_intricate"
testcaseid = "sncf_intricate"
# testcaseid = "tripod_chattering"

### Ellipses 

# testcaseid = "eikonale_ellcube"
# testcaseid = "eikonale_ellplane"
# testcaseid = "robustHam" # A(u) = [u 0; 0 -u], u ∈ {-1,1}
# testcaseid = "robustRot" # A(u) = [0 1; u 0], u ∈ {-1,0,1} 

### Gluings

# testcaseid = "eikonale_gluing"
# testcaseid = "eikonale_planeline"
# testcaseid = "windrobot_planeline"
# testcaseid = "eikonale_ISS"
# testcaseid = "windrobot_ISS"

verbose = true

picruns = [1,20,50,100]
# picruns = [1]
# pictype = "value"
# pictype = "feedback"
# pictype = "opttraj"
pictype = ["value", "feedback"]

#####################################################
# Initialization 
#####################################################

testcase = get_TestCase(testcaseid)

if testcase.domaintype == "network"
    domain, plotter = get_network_and_plotter(testcase.domainid)
elseif testcase.domaintype == "ellipses"
    domain = get_ellipses(testcase.domainid, verbose)
elseif testcase.domaintype == "gluing"
    domain = get_gluing(testcase.domainid, verbose)
    plotter = get_gluing_plotter(testcase.domainid, verbose)
else 
    throw(ErrorException("what is domain type $(testcase.domaintype) exactly"))
end
dynamic = get_Dynamic(testcase.dynamicid, domain)
terminalcost = get_Scalar(testcase.costid, domain)
if save_errors
    has_analytical = get_Analytical(domain, testcase.dynamicid, testcase.costid, testcase.T, true)
    if has_analytical
        analytical = get_Analytical(domain, testcase.dynamicid, testcase.costid, testcase.T, false, true)
    end
end
initpoints = get_init_points(testcaseid, domain)
foldername = join(split(testcase.runid,"_")[1:2],"_")
if !isdir("data/$foldername")
    mkdir("data/$foldername")
end
if !isdir("data/$foldername/hatV")
    mkdir("data/$foldername/hatV")
end

#####################################################
# Iterations
#####################################################

for (iN,(N,dx)) in enumerate(zip(testcase.NN, testcase.dxs))

    if (length(testcase.NN) >= 4)
        verb(verbose, "\n########### Run $iN / $(length(testcase.NN)) ###########")
    end
    mesh = get_mesh(domain, dx, verbose)
    testcase.npoints[iN] = mesh.npoints

    if save_errors && has_analytical 
        verb(verbose, "Creating analytical VV for N = $N...")
        VV = hcat([[analytical.call((n-1)*testcase.T/N,mesh.points[j]) for j in 1:mesh.npoints] for n in 1:N+1]...)
        verb(verbose, "... done.")
    end 

    ##########################
    # SemiLag part
    ##########################

    if VSL 

        if testcase.scheme == "SL"
            testcase.exec_time[iN] = @elapsed hatVV = semiLag(domain, dynamic, terminalcost, testcase.T, N, mesh, verbose)
        elseif testcase.scheme == "SLinterp"
            testcase.exec_time[iN] = @elapsed hatVV = semiLagInterp(domain, dynamic, terminalcost, testcase.T, N, mesh, verbose)
        else
            throw(ErrorException("unknown method $(testcase.method)"))
        end

        if (uSL || create_pics)
            if testcase.buildV == "closestN"
                # closest neighbour. 
                function hatV(t::Float64, x::Point)
                    tn = 1+round(Int64, N*t/testcase.T)
                    ix = search_in_mesh(domain, mesh, x)
                    # println("ix : ", ix)
                    # println("At time $tn, choosing ", pt_to_str(mesh.points[ix]))
                    return hatVV[ix,tn]
                end
            
            else 
                throw(ErrorException("unknown way to build V: $(testcase.buildV)"))
            end
        end
        
        if save_errors && has_analytical
            verb(verbose, "Summing errors for N = $N...")
            @views testcase.errorVSL_glob[iN] = maximum(abs.(hatVV .- VV)) / maximum(abs.(VV))
            @views testcase.errorVSL_init[iN] = maximum(abs.(hatVV[:,1] .- VV[:,1])) / maximum(abs.(VV[:,1]))
            verb(verbose, "Errors on V: global ", testcase.errorVSL_glob[iN], ", initial ", testcase.errorVSL_init[iN])

            # tmax = argmax([maximum(abs.(hatVV[:,tn] .- VV[:,tn])) for tn in 1:N])
            # println("Time of the maximal error : $tmax / $N")

            verb(verbose, "... done.")
        end

        # if mesh.npoints <= 1000
        #     submesh = collect(1:mesh.npoints)
        # else
        #     submesh = mesh.juncind ∪ collect(1:round(Int,mesh.npoints/1000):mesh.npoints) 
        #     verb(verbose, "Plotting on a submesh of ", length(submesh), " points")
        # end
        # # tmax = 1
        # p = plot_function(domain, plotter, mesh, hatVV[:,tmax], submesh, verbose)
        # p = plot_function(domain, plotter, mesh, VV[:,tmax], submesh, verbose)
        # p = plot_function(domain, plotter, mesh, hatVV[:,tmax] .- VV[:,tmax], submesh, verbose)
        # display(p)
        # display(plot(p,q,layout=[1,1]))

        # p = plot_gluing(domain, plotter)
        # q = plot_gluing(domain, plotter)
        # plot_values!(p, domain, plotter, mesh, hatVV[:,tmax], 1.5)
        # plot_values!(q, domain, plotter, mesh, VV[:,tmax], 1.5, true)
        # plot_values!(p, domain, plotter, mesh, hatVV[:,tmax] .- VV[:,tmax], 1.5)
        # display(p)
        # display(plot(p,q,layout=[1,1]))

        if create_pics && iN ∈ picruns
            create_and_save_pics(testcase, testcaseid, foldername, domain, plotter, dynamic, mesh, hatVV, hatV, pictype, iN, N, dx, verbose)
        end

        if save_values
            writedlm("data/$foldername/hatV/$(testcase.runid)_hatV$iN", hatVV)
        end

        if uSL 
            compute_errors_controls!(testcase, domain, dynamic, hatV, iN, N, initpoints, save_errors, has_analytical, verbose)
        end

    end
end

#####################################################
# Diagnostics
#####################################################

save_testcase(testcase, verbose)
if (save_errors && has_analytical && length(testcase.NN) >= 5)
    run_diagnostics(testcase, VSL, uSL, verbose)
end

println("Bye")