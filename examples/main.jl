using FLagHada
using Parameters # with_kw
using DocStringExtensions
using ProgressBars
using Plots
using DelimitedFiles
using JLD2

include("testcases.jl")
include("visualisation.jl")

println("Welcome in main.")

# To the passerby: 
# This is the full code used in our work on CAT(0) spaces. 
# Feel free to look at it: if you are actually interested into running it (who knows), 
# you may prefer to start with the demo test cases.

# How to run it: 
# uncomment the testcaseid that you want 
# in the REPL: 
#   ] activate .
#   using FLagHada
#   include("examples/main.jl") 

# Warnings: 
# Some error plots may be empty if the error is of machine precision
# In particular for the demo examples...
# The first run takes longer.

#####################################################
# Parameters 
#####################################################

VSL = true # semi-lagrangian on the value function 
uSL = true # compute errors on the controls as well
save_errors = true
save_values = true
create_pics = true

### Demo test cases (should run in reasonable time)

testcaseid = "demo_network" 
# testcaseid = "demo_ellipses"
# testcaseid = "demo_gluing"

### Networks

# testcaseid = "eikonale_tripod"
# testcaseid = "eikonale_intricate"
# testcaseid = "sncf_intricate"
# testcaseid = "tripod_chattering"

### Ellipses 

# testcaseid = "eikonale_ellcube"
# testcaseid = "eikonale_ellplane"
# testcaseid = "robustHam" # A(u) = [u 0; 0 -u], u ∈ {-1,1}

### Gluings

# testcaseid = "eikonale_gluing"
# testcaseid = "eikonale_planeline"
# testcaseid = "windrobot_planeline"
# testcaseid = "eikonale_ISS"
# testcaseid = "windrobot_ISS"

verbose = true

# index at which we make pictures
picruns = [1,2,3,4,5] # for demos
# picruns = [1,20,50,100] 
# picruns = [1]
# pictype = "value"
# pictype = "feedback"
# pictype = "opttraj"
pictype = ["value", "feedback"] # can be "value", "feedback" and "opttraj". !! No guarantee to be implemented, check it yourself

#####################################################
# Initialization 
#####################################################

testcase = get_TestCase(testcaseid)

if testcase.domaintype == "network"
    domain, plotter = get_network_and_plotter(testcase.domainid)
elseif testcase.domaintype == "ellipses"
    domain = get_ellipses(testcase.domainid, verbose)
    plotter = "not used"
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

    if VSL 

        if testcase.scheme == "SL"
            testcase.exec_time[iN] = @elapsed hatVV = semiLag(domain, dynamic, terminalcost, testcase.T, N, mesh, verbose)
        else
            throw(ErrorException("unknown method $(testcase.method)"))
        end

        if (uSL || create_pics)
            if testcase.buildV == "closestN"
                # closest neighbour. 
                function hatV(t::Float64, x::Point)
                    tn = 1+round(Int64, N*t/testcase.T)
                    ix = search_in_mesh(domain, mesh, x)
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
            verb(verbose, "... done.")
        end

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