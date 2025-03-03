# @testset "convex hulls" begin
#     domain, plotter = get_network_and_plotter("tripod")

#     # 1 point 
#     @test get_convex(domain, [domain.juncpoints[1]]) == Network_convex([domain.juncpoints[1]])
#     # 2 points, no problem 
#     @test get_convex(domain, [domain.juncpoints[1],domain.juncpoints[2]]) == Network_convex([domain.juncpoints[1],domain.juncpoints[2]])
#     # 2 points + trash
#     @test get_convex(domain, [domain.juncpoints[1],domain.juncpoints[2],EdgeP(1,2,0.5,0.5)]) == Network_convex([domain.juncpoints[1],domain.juncpoints[2]])
#     # 3 edge points, no problem
#     points = [EdgeP(1,2,0.5,0.5),EdgeP(1,3,0.5,0.5),EdgeP(1,4,0.5,0.5)]
#     solution = Network_convex(points)
#     @test get_convex(domain, points) == solution 
#     # 3 junc points, no problem
#     points = [domain.juncpoints[2],domain.juncpoints[3],domain.juncpoints[4]]
#     solution = Network_convex(points)
#     @test get_convex(domain, points) == solution 
#     # A lot of points
#     @test get_convex(domain, union(points, [domain.juncpoints[1],EdgeP(1,2,0.5,0.5), EdgeP(1,2,0.4,0.6), EdgeP(1,4,0.7,0.3)])) == solution
# end

# @testset "distances to convex hulls" begin 
#     domain, plotter = get_network_and_plotter("tripod")

#     # 1 point 
#     @test get_dist_to_convex(domain, Network_convex([domain.juncpoints[1]]), domain.juncpoints[1]) ≈ 0.0
#     @test get_dist_to_convex(domain, Network_convex([domain.juncpoints[1]]), EdgeP(1, 2, 0.5,0.5)) ≈ 0.5
#     @test get_dist_to_convex(domain, Network_convex([domain.juncpoints[1]]), domain.juncpoints[2]) ≈ 1.0

#     # 2 points 
#     convex = Network_convex([EdgeP(1,4,0.76,0.24), EdgeP(1,4,0.9,0.1)])
#     @test get_dist_to_convex(domain, convex, EdgeP(1,4,0.8,0.2)) ≈ 0.0
#     @test get_dist_to_convex(domain, convex, domain.juncpoints[4]) ≈ 0.1
#     @test get_dist_to_convex(domain, convex, domain.juncpoints[1]) ≈ 0.76
#     @test get_dist_to_convex(domain, convex, domain.juncpoints[2]) ≈ 1.76

#     # more points 
#     convex = Network_convex([EdgeP(1,2,0.5,0.5),EdgeP(1,3,0.5,0.5),EdgeP(1,4,0.5,0.5)])
#     @test get_dist_to_convex(domain, convex, domain.juncpoints[1]) ≈ 0.0
#     @test get_dist_to_convex(domain, convex, EdgeP(1,2,0.4,0.6)) ≈ 0.0
#     @test get_dist_to_convex(domain, convex, EdgeP(1,2,0.6,0.4)) ≈ 0.1
#     @test get_dist_to_convex(domain, convex, domain.juncpoints[2]) ≈ 0.5
# end

# @testset "directions" begin 
#     domain, plotter = get_network_and_plotter("tripod")

#     # Junction directions 
#     xx = domain.juncpoints[1] # central junction 
#     @test get_direction(domain, xx, domain.juncpoints[2]) == 2
#     @test get_direction(domain, xx, domain.juncpoints[3]) == 3
#     @test get_direction(domain, xx, domain.juncpoints[4]) == 4
#     @test get_direction(domain, xx, EdgeP(1,2,0.5,0.5)) == 2
#     @test get_direction(domain, xx, EdgeP(1,3,0.5,0.5)) == 3
#     @test get_direction(domain, xx, EdgeP(1,4,0.5,0.5)) == 4
#     @test get_direction(domain, xx, EdgeP(2,1,0.5,0.5)) == 2
#     @test get_direction(domain, xx, EdgeP(3,1,0.5,0.5)) == 3
#     @test get_direction(domain, xx, EdgeP(4,1,0.5,0.5)) == 4

#     xx = domain.juncpoints[2]
#     @test get_direction(domain, xx, domain.juncpoints[1]) == 1
#     @test get_direction(domain, xx, domain.juncpoints[3]) == 1
#     @test get_direction(domain, xx, domain.juncpoints[4]) == 1
#     @test get_direction(domain, xx, EdgeP(1,2,0.5,0.5)) == 1
#     @test get_direction(domain, xx, EdgeP(1,3,0.5,0.5)) == 1
#     @test get_direction(domain, xx, EdgeP(1,4,0.5,0.5)) == 1
#     @test get_direction(domain, xx, EdgeP(2,1,0.5,0.5)) == 1
#     @test get_direction(domain, xx, EdgeP(3,1,0.5,0.5)) == 1
#     @test get_direction(domain, xx, EdgeP(4,1,0.5,0.5)) == 1

#     # Edge directions 
#     xx = EdgeP(1,2,0.5,0.5)
#     @test get_direction(domain, xx, domain.juncpoints[1]) == 1
#     @test get_direction(domain, xx, domain.juncpoints[2]) == 2
#     @test get_direction(domain, xx, domain.juncpoints[3]) == 1
#     @test get_direction(domain, xx, EdgeP(1,2,0.75,0.25)) == 2
#     @test get_direction(domain, xx, EdgeP(1,2,0.25,0.75)) == 1
#     @test get_direction(domain, xx, EdgeP(1,3,0.75,0.25)) == 1
#     @test get_direction(domain, xx, EdgeP(2,1,0.75,0.25)) == 1 
#     @test get_direction(domain, xx, EdgeP(2,1,0.25,0.75)) == 2
#     @test get_direction(domain, xx, EdgeP(3,1,0.75,0.25)) == 1
# end

# @testset "mesh" begin 
#     domain, plotter = get_network_and_plotter("tripod")

#     mesh = get_mesh(domain, 2.0) # so large that no interior points 
#     # junction test 
#     @test all([mesh.points[search_in_mesh(domain, mesh, x)] == x for x in domain.juncpoints])
#     # edge point test 
#     @test all([mesh.points[search_in_mesh(domain, mesh, x)] == x for x in mesh.points if typeof(x) == EdgeP])
#     # other points test 
#     @test search_in_mesh(domain, mesh, EdgeP(1,2,0.2,0.8)) == 1 
#     @test search_in_mesh(domain, mesh, EdgeP(1,2,0.5,0.5)) in [1,2] 
#     @test search_in_mesh(domain, mesh, EdgeP(1,2,0.55,0.45)) == 2 

#     mesh = get_mesh(domain, 0.06) # now interior points 
#     # junction test 
#     @test all([mesh.points[search_in_mesh(domain, mesh, x)] == x for x in domain.juncpoints])
#     # edge point test 
#     @test all([mesh.points[search_in_mesh(domain, mesh, x)] == x for x in mesh.points if typeof(x) == EdgeP])
#     # other points test 
#     thepoint = EdgeP(1,2,pi/10,1.0-pi/10)
#     thefound = search_in_mesh(domain, mesh, thepoint)
#     println("Point : ", thepoint, ", found : ", mesh.points[thefound])
#     @test get_distance(domain, thepoint, mesh.points[thefound]) <= minimum([get_distance(domain, thepoint, x) for x in mesh.points])
#     thepoint = EdgeP(4,1,pi/10,1.0-pi/10)
#     thefound = search_in_mesh(domain, mesh, thepoint)
#     @test get_distance(domain, thepoint, mesh.points[thefound]) <= minimum([get_distance(domain, thepoint, x) for x in mesh.points])
# end

@testset "mesh_conv" begin 
    domain, plotter = get_network_and_plotter("tripod")

    mesh = get_mesh(domain, 0.02)
    # all mesh test 
    @test sort(get_mesh_convex(domain, mesh, domain.juncpoints)) == 1:mesh.npoints 
    # one edge test
    @test sort(get_mesh_convex(domain, mesh, [domain.juncpoints[1], domain.juncpoints[4]])) == union(1, 102:151, 154)
    # one side of edge test 
    @test sort(get_mesh_convex(domain, mesh, [EdgeP(1,4,0.5,0.5), domain.juncpoints[4]])) == union(126:151, 154)
    # strictly inner test 
    @test sort(get_mesh_convex(domain, mesh, [EdgeP(1,4,0.5,0.5), EdgeP(1,4,0.75,0.25)])) == union(126:139)
    # around one junction test 
    @test sort(get_mesh_convex(domain, mesh, [EdgeP(1,2,0.5,0.5), EdgeP(1,3,0.5,0.5), EdgeP(1,4,0.5,0.5)])) == union(1:26, 52:76, 102:126)

end