using Dates 
using JLD2 # format of saved files 

@with_kw struct TestCase 
    # immutable part
    runid::String
    domaintype::String
    domainid::String 
    dynamicid::String 
    costid::String
    scheme::String
    buildV::String
    T::Float64 
    dxs::Vector{Float64}
    NN::Vector{Int64}

    # mutable part
    "Number of mesh points"
    npoints::Vector{Int64}
    "Supremum on time-space of the error" 
    errorVSL_glob::Vector{Float64}
    "Supremum on space, initial time" 
    errorVSL_init::Vector{Float64}
    "Loss of optimality, worst case"
    erroruSL_lossw::Vector{Float64}
    "Loss of optimality, mean"
    erroruSL_lossm::Vector{Float64}
    "Renormalized loss of optimality, worst case" 
    erroruSL_percw::Vector{Float64}
    "Renormalized loss of optimality, mean" 
    erroruSL_percm::Vector{Float64}
    "Execution time of the runs"
    exec_time::Vector{Float64}

    function TestCase(domaintype::String, domainid::String, dynamicid::String, costid::String, scheme::String, buildV::String,
                T::Float64, init::Float64, last::Float64, number::Int64, exponent::Int64, fact::Float64=1.0)
        dxs = init * exp.(log(last/init) .* (1:number) ./ number)
        NN = [ceil(Int64, fact*T*(1.0/dx)^(1/exponent)) for dx in dxs]
        runid = domainid[1:min(end,5)] *"_"* dynamicid[1:min(end,5)] *"_"* costid[1:min(end,5)] *"_T"* replace(numtostr(T,5),"."=>"-") *"_"*
            "Nmin$(NN[1])_Nmax$(NN[end])_" * Dates.format(now(), "dd-mm-yy_HHhMM")
        return new(runid, domaintype, domainid, dynamicid, costid, scheme, buildV, T, dxs, NN, [1 for _ in 1:length(NN)], [fill(42.0,length(NN)) for _ in 1:7]...)
    end
end

"""
    Bank of test cases.
"""
function get_TestCase(testcaseid::String)

    ##############################
    # Demos 
    ##############################

    # params: 
    #   domaintype      "network", "ellipses" or "gluing" (feel free to add some, requirements in FLagHada.jl)
    #   domainid        passed to the constructor of the specific type of domain 
    #   dynamicid       passed to get_Dynamic
    #   costid          passed to get_Scalar
    #   scheme          always "SL" for now
    #   buildV          always "closestN" for now
    #   T               terminal time, a positive float 
    #   init            largest  space step for the space mesh, a positive float
    #   last            smallest space step for the space mesh, a positive float
    #   number          number of space meshes. Step will be distributed log-uniformly
    #   exponent        2. For the CFL condition
    #   fact (opt)      default 1.0, larger = more time steps = smaller Δt

    if (testcaseid == "demo_network")

        return TestCase("network", "tripod", "eikonal", "norm", "SL", "closestN", 0.8, 0.1, 0.01, 5, 2)

    elseif (testcaseid == "demo_ellipses")

        return TestCase("ellipses", "ellctedet", "robustHam", "eigmax", "SL", "closestN", 1.0, 0.7, 0.05, 5, 2)
    
    elseif (testcaseid == "demo_gluing")

        return TestCase("gluing", "planeline", "eikonal", "norm", "SL", "closestN", 1.0, 0.1, 0.01, 5, 2)

    ##############################
    # Networks
    ##############################
    
    elseif (testcaseid == "eikonale_tripod")

        return TestCase("network", "tripod", "eikonal", "norm", "SL", "closestN", 0.5, 0.05, 0.000005, 100, 2) # big simulation

    elseif (testcaseid == "eikonale_intricate")
        
        return TestCase("network", "intricate", "eikonal", "norm", "SL", "closestN", 6.0, 0.08, 0.01, 10, 2)

    elseif (testcaseid == "sncf_intricate")
        
        return TestCase("network", "intricate", "sncf", "dist_south", "SL", "closestN", 6.0, 0.08, 0.00005, 100, 2) # big simulation

    elseif (testcaseid == "tripod_chattering")

        return TestCase("network", "tripod", "tripod_twoaims", "norm", "SL", "closestN", 2.0, 0.1, 0.001, 500, 2)

    ##############################
    # Ellipses
    ##############################

    elseif (testcaseid == "eikonale_ellcube")

        return TestCase("ellipses", "ellcube", "eikonal", "norm", "SL", "closestN", 0.5, 0.7, 0.1, 5, 2)

    elseif (testcaseid == "eikonale_ellplane")

        return TestCase("ellipses", "ellplane", "eikonal", "norm", "SL", "closestN", 1.0, 0.1, 0.01, 5, 2, 4.0)

    elseif (testcaseid == "robust")

        return TestCase("ellipses", "ellcube", "robust", "eigmax", "SL", "closestN", 1.0, 0.7, 0.03, 10, 2) 

    elseif (testcaseid == "robustHam")

        return TestCase("ellipses", "ellctedet", "robustHam", "eigmax", "SL", "closestN", 0.3, 0.08, 0.005, 100, 2, 6.0) # big simulation

    ##############################
    # Gluings
    ##############################

    elseif (testcaseid == "eikonale_gluing")

        return TestCase("gluing", "2dsquare", "eikonal", "norm", "SL", "closestN", 0.7, 0.1, 0.001, 10, 2, 2.0)

    elseif (testcaseid == "eikonale_planeline")

        return TestCase("gluing", "planeline", "eikonal", "norm", "SL", "closestN", 1.0, 0.012, 0.001, 100, 2)

    elseif (testcaseid == "windrobot_planeline")

        return TestCase("gluing", "planeline", "windrobot", "distmeteors", "SL", "closestN", 1.8, 0.2, 0.002, 100, 2, 1.5) # big simulation

    elseif (testcaseid == "eikonale_ISS")

        return TestCase("gluing", "ISS", "eikonal", "norm", "SL", "closestN", 3.0, 1.0, 0.03, 10, 2)

    elseif (testcaseid == "windrobot_ISS")

        return TestCase("gluing", "ISS", "windrobot", "distmeteors", "SL", "closestN", 2.5, 0.2, 0.0045, 100, 2) # big simulation

    else 
        throw(ErrorException("unknown test case $testcaseid"))
    end
end

"""
    Return the initial points of trajectories 
    for the numerical approximation 
"""
function get_init_points(testcaseid::String, domain::Domain)

    ##############################
    # Demos 
    ##############################

    if (testcaseid == "demo_network")

        return domain.juncpoints

    elseif (testcaseid == "demo_ellipses")

        return [EPoint(domain.alga,domain.beta,domain.alga)]

    elseif (testcaseid == "demo_gluing")

        return [GPoint(2,[1.0])]

    ##############################
    # Networks
    ##############################

    elseif (testcaseid in ["eikonale_tripod", "eikonale_intricate", "tripod_chattering", "sncf_intricate"])

        return domain.juncpoints

    ##############################
    # Ellipses
    ##############################

    elseif (testcaseid in ["eikonale_ellcube", "eikonale_ellplane"])

        return [EPoint(domain.alga,domain.beta,domain.alga)]

    elseif (testcaseid == "robustHam")

        ABG = [rand(3) .- 0.5 for _ in 1:10]
        for abg in ABG
            abg[3] = - abg[1]
        end
        return [EPoint(abg...) for abg in ABG]

    ##############################
    # Gluings
    ##############################

    elseif testcaseid == "eikonale_gluing"

        return [GPoint(1,ones(domain.comps[1].dim))]

    elseif testcaseid == "eikonale_planeline"

        return [GPoint(2,[1.0])]

    elseif testcaseid == "windrobot_planeline"

        return [GPoint(1,[0.0,4.0]), GPoint(1,[0.0,0.0]), GPoint(1,[4.0,0.0]), GPoint(1,[4.0,4.0]), 
                GPoint(2,[4.0]), GPoint(1,[2.0, 2.0]), GPoint(1,[4.0, 2.0]), GPoint(2,[2.0])]

    elseif testcaseid in ["eikonale_ISS", "windrobot_ISS"]

        return [GPoint(1, [0.0, 3.0]), GPoint(2, [0.0, 3.0]), GPoint(3, [0.0, 3.0]), GPoint(4, [0.0, 3.0]),
                GPoint(9, [1.0, 0.0]), GPoint(10, [0.0, 1.5, 0.0])]

    else 
        throw(ErrorException("unknown test case id $testcaseid"))
    end

end

function save_testcase(testcase::TestCase, verbose::Bool=false; suffix::String="")
    verb(verbose, "Saving testcase in $(testcase.runid)...")
    foldername = join(split(testcase.runid,"_")[1:2],"_")
    # julia-friendly format 
    save_object("data/$foldername/" * testcase.runid * suffix * ".jld2", testcase) 

    # [whatever-is-not-julia]-friendly format
    open("data/$foldername/" * testcase.runid * suffix * ".txt", "w") do textfile
        # column names 
        write(textfile, "runid"); write(textfile, ("\t"))
        write(textfile, "domaintype"); write(textfile, ("\t"))
        write(textfile, "domainid"); write(textfile, ("\t"))
        write(textfile, "dynamicid"); write(textfile, ("\t"))
        write(textfile, "costid"); write(textfile, ("\t"))
        write(textfile, "scheme"); write(textfile, ("\t"))
        write(textfile, "buildV"); write(textfile, ("\t"))
        write(textfile, "T"); write(textfile, ("\t"))
        write(textfile, "dx"); write(textfile, ("\t"))
        write(textfile, "N"); write(textfile, ("\t"))
        write(textfile, "npoints"); write(textfile, ("\t"))
        write(textfile, "errorVSLglob"); write(textfile, ("\t"))
        write(textfile, "errorVSLinit"); write(textfile, ("\t"))
        write(textfile, "erroruSLlossw"); write(textfile, ("\t"))
        write(textfile, "erroruSLlossm"); write(textfile, ("\t"))
        write(textfile, "erroruSLpercw"); write(textfile, ("\t"))
        write(textfile, "erroruSLpercm"); write(textfile, ("\t"))
        write(textfile, "exectime"); write(textfile, ("\n"))
        # contents
        for iN in 1:length(testcase.NN)
            write(textfile, testcase.runid); write(textfile, ("\t"))
            write(textfile, testcase.domaintype); write(textfile, ("\t"))
            write(textfile, testcase.domainid); write(textfile, ("\t"))
            write(textfile, testcase.dynamicid); write(textfile, ("\t"))
            write(textfile, testcase.costid); write(textfile, ("\t"))
            write(textfile, testcase.scheme); write(textfile, ("\t"))
            write(textfile, testcase.buildV); write(textfile, ("\t"))
            write(textfile, "$(testcase.T)"); write(textfile, ("\t"))
            write(textfile, "$(testcase.dxs[iN])"); write(textfile, ("\t"))
            write(textfile, "$(testcase.NN[iN])"); write(textfile, ("\t"))
            write(textfile, "$(testcase.npoints[iN])"); write(textfile, ("\t"))
            write(textfile, "$(testcase.errorVSL_glob[iN])"); write(textfile, ("\t"))
            write(textfile, "$(testcase.errorVSL_init[iN])"); write(textfile, ("\t"))
            write(textfile, "$(testcase.erroruSL_lossw[iN])"); write(textfile, ("\t"))
            write(textfile, "$(testcase.erroruSL_lossm[iN])"); write(textfile, ("\t"))
            write(textfile, "$(testcase.erroruSL_percw[iN])"); write(textfile, ("\t"))
            write(textfile, "$(testcase.erroruSL_percm[iN])"); write(textfile, ("\t"))
            write(textfile, "$(testcase.exec_time[iN])"); write(textfile, ("\n"))
        end
    end # closes textfile 
    verb(verbose, "... done.")
end

function load_testcase(runid::String, verbose::Bool=false)
    verb(verbose, "Loading " * runid * "...")
    foldername = join(split(runid,"_")[1:2],"_")
    testcase = load_object("data/$foldername/" * runid * ".jld2")
    verb(verbose, "... done.")
    return testcase
end

# side comment: in the simulations, I use plot_errors_robust instead of plot_errors
# which can treat outliers. The code is provided in interface.jl, I just do not want to impose LinRegOutliers
function run_diagnostics(testcase::TestCase, VSL::Bool, uSL::Bool, verbose::Bool=false; suffix::String="")
    foldername = join(split(testcase.runid,"_")[1:2],"_")
    # Errors on the value function 
    if VSL
        verb(verbose, "Plotting error on the value function for test case $(testcase.runid)...")
        savefig(plot_errors(testcase.dxs, testcase.errorVSL_glob, ["Relative global error on V w.r.t. the space step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_errVSL_dx" * suffix)
        savefig(plot_errors(testcase.T./testcase.NN, testcase.errorVSL_glob, ["Relative global error on V w.r.t. the time step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_errVSL_dt" * suffix)
        verb(verbose, "... done.")
    end

    if uSL
        verb(verbose, "Plotting error on the optimal trajectory for test case $(testcase.runid)...")
        # Errors on the optimal trajectories
        savefig(plot_errors(testcase.T./testcase.NN, testcase.erroruSL_lossw, ["Worst-case loss of optimality w.r.t. the time step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_erruSLw_dt" * suffix)
        savefig(plot_errors(testcase.T./testcase.NN, testcase.erroruSL_percw, ["Worst-case renormalized loss of optimality w.r.t. the time step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_rerruSLw_dt" * suffix)
        savefig(plot_errors(testcase.T./testcase.NN, testcase.erroruSL_lossm, ["Mean loss of optimality w.r.t. the time step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_erruSLm_dt" * suffix)
        savefig(plot_errors(testcase.T./testcase.NN, testcase.erroruSL_percm, ["Mean renormalized loss of optimality w.r.t. the time step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_rerruSLm_dt" * suffix)
        savefig(plot_errors(testcase.dxs,            testcase.erroruSL_lossw, ["Worst-case loss of optimality w.r.t. the space step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_erruSLw_dx" * suffix)
        savefig(plot_errors(testcase.dxs,            testcase.erroruSL_percw, ["Worst-case renormalized loss of optimality w.r.t. the space step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_rerruSLw_dx" * suffix)
        savefig(plot_errors(testcase.dxs,            testcase.erroruSL_lossm, ["Mean loss of optimality w.r.t. the space step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_erruSLm_dx" * suffix)
        savefig(plot_errors(testcase.dxs,            testcase.erroruSL_percm, ["Mean renormalized loss of optimality w.r.t. the space step, $(testcase.scheme) scheme"], verbose), "data/$foldername/$(testcase.runid)_rerruSLm_dx" * suffix)
        verb(verbose, "... done.")
    end
end