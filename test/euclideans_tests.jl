@testset "convex hulls" begin
    # 1D test
    domain = Domain_euclidean("segment", 1, [RdPoint([0.0]), RdPoint([1.0])])
    @test get_convex(domain, [RdPoint([0.2]), RdPoint([0.3]), RdPoint([0.8])]) == Euclidean_convex_1D(0.2, 0.8)

    # 2D test, one point
    domain = Domain_euclidean("square", 2, [RdPoint([0.0,0.0]), RdPoint([1.0,1.0])])
    @test get_convex(domain, [RdPoint([0.1, 0.1])]) == Euclidean_convex_2D_singlepoint(RdPoint([0.1, 0.1]))

    # 2D test, two points
    @test get_convex(domain, [RdPoint([0.1, 0.1]), RdPoint([0.1, 0.8])]) == Euclidean_convex_2D_segment(RdPoint([0.1, 0.1]), RdPoint([0.1, 0.8]))

    # 2D test, three points without problems 
    points = [RdPoint([0.1, 0.1]), RdPoint([0.1, 0.8]), RdPoint([0.8, 0.8])]
    solution = Euclidean_convex_2D_fullset([RdPoint([0.1, 0.1]), RdPoint([0.8, 0.8]), RdPoint([0.1, 0.8])])
    @test get_convex(domain, points) == solution 

    # 2D test, three points, adding useless inner points
    @test get_convex(domain, union(points, [RdPoint([0.2, 0.7]), RdPoint([0.7, 0.78]), RdPoint([0.5, 0.7])])) == solution 

    # 2D test, more points 
    points = [RdPoint([0.5,0.5] .+ 0.4 * [cos(pi*1.1+2*pi*k/5),sin(pi*1.1+2*pi*k/5)]) for k in 0:4]
    solution = Euclidean_convex_2D_fullset(points)
    innerpoints = [RdPoint([0.5,0.5] .+ 0.2 * s * [1.0,1.0]) for s in [-1.0,1.0]]
    @test get_convex(domain, union(points, innerpoints)) == solution

    # 2D test, strictly aligned points 
    points = [RdPoint([0.1,0.1]), RdPoint([0.2,0.2]), RdPoint([0.5,0.5]), RdPoint([0.8,0.8])]
    solution = Euclidean_convex_2D_fullset([points[1],points[end]])
    @test get_convex(domain, points) == solution

    # 2D test, non-strictly aligned points 
    points = [RdPoint([0.1, 0.1]), RdPoint([0.8, 0.1]), RdPoint([0.8, 0.8]), RdPoint([0.1, 0.8])]
    solution = Euclidean_convex_2D_fullset(points)
    innerpoints = [RdPoint([0.399999999,0.4])]
    @test get_convex(domain, union(points, innerpoints)) == solution
end

# @testset "distances to convex hulls" begin 
#     # 1D test
#     domain = Domain_euclidean("segment", 1, [RdPoint([0.0]), RdPoint([1.0])])
#     convex = Euclidean_convex_1D(0.2, 0.8)
#     @test get_dist_to_convex(domain, convex, RdPoint([0.2])) ≈ 0.0
#     @test get_dist_to_convex(domain, convex, RdPoint([0.2])) ≈ 0.0
#     @test get_dist_to_convex(domain, convex, RdPoint([0.0])) ≈ 0.2
#     @test get_dist_to_convex(domain, convex, RdPoint([1.0])) ≈ 0.2

#     # 2D test, one point
#     domain = Domain_euclidean("square", 2, [RdPoint([0.0,0.0]), RdPoint([1.0,1.0])])
#     convex = Euclidean_convex_2D_singlepoint(RdPoint([0.4, 0.8]))
#     @test get_dist_to_convex(domain, convex, RdPoint([0.4, 0.8])) ≈ 0.0
#     @test get_dist_to_convex(domain, convex, RdPoint([0.4, 0.4])) ≈ 0.4

#     # 2D test, two points
#     convex = Euclidean_convex_2D_segment(RdPoint([0.1, 0.1]), RdPoint([0.1, 0.8]))
#     @test get_dist_to_convex(domain, convex, RdPoint([0.1, 0.5])) ≈ 0.0
#     @test get_dist_to_convex(domain, convex, RdPoint([0.8, 0.8])) ≈ 0.7

#     # 2D test, more points 
#     points = [RdPoint([0.5,0.5] .+ 0.4 * [cos(pi*1.1+2*pi*k/5),sin(pi*1.1+2*pi*k/5)]) for k in 0:4]
#     convex = Euclidean_convex_2D_fullset(points)
#     @test get_dist_to_convex(domain, convex, RdPoint([0.5, 0.5])) ≈ 0.0
#     @test get_dist_to_convex(domain, convex, RdPoint([0.5,0.5] .+ 0.6 * [cos(pi*1.1),sin(pi*1.1)])) ≈ 0.2
# end