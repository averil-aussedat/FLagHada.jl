using Plots

@testset "domains" begin 
    
    # domainid = "1dsegment"
    # domainid = "2dsquare"
    # domainid = "3dcube"
    domainid = "ISS"

    domain = get_gluing(domainid)
    plotter = get_gluing_plotter(domainid)
    p = plot_gluing(domain, plotter)

    p1 = GPoint(1, [-0.25, -2.0])
    # p2 = GPoint(3, [0.25, 2.0])
    p2 = GPoint(10, [-0.25, 1.0, 0.0])
    vv = Velocity(p2, 1.0, 1.0)
    Ng = 40
    geodesic = [get_exponential(domain, vv, p1, h) for h in LinRange(0.0, get_distance(domain, p1, p2), Ng)]
    # println("geodesic : ", geodesic)
    for i in 1:Ng
        println("point $i : ", pt_to_str(geodesic[i]))
        if i < Ng 
            println("\tdistance to next : ", get_distance(domain, geodesic[i], geodesic[i+1]))
        end
    end
    plot_points!(p, domain, plotter, geodesic, ms=1)
    plot!(p, dpi=400)
    savefig(p, "data/ISS")

    # mesh = get_mesh(domain, 0.45)
    # plot_points!(p, domain, plotter, mesh.points, ms=1)

    # # pt = GPoint(1, [0.01, 2.0])
    # pt = GPoint(5, [0.7])
    # ipt = search_in_mesh(domain, mesh, pt)
    # plot_points!(p, domain, plotter, [pt], color="red")
    # plot_points!(p, domain, plotter, [mesh.points[ipt]], color="green")

    # pts = [
    #     GPoint(1, [-0.5, -2.0]),
    #     GPoint(1, [0.5, -2.0]),
    #     GPoint(1, [-0.5, 1.0]),
    #     # GPoint(5, [0.0]),
    #     GPoint(4, [0.5, -3.0]),
    #     GPoint(4, [0.5, 3.0]),
    # ]
    # ipts = get_mesh_convex(domain, mesh, pts)
    # println("type ipts : ", typeof(ipts))
    # println("ipts : ", ipts)
    # plot_points!(p, domain, plotter, pts, color="red")
    # plot_points!(p, domain, plotter, [mesh.points[ipt] for ipt in ipts], color="green")

    display(p)
end
