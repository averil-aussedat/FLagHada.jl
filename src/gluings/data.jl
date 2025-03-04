export get_gluing
function get_gluing(id::String, verbose::Bool=false)
    verb(verbose, "Creating $id...")

    if id in ["1dsegment", "2dsquare", "3dcube"]

        if id == "1dsegment"
            dim = 1
        elseif id == "2dsquare"
            dim = 2
        else
            dim = 3
        end

        # id, ncomps, comps, get_crossing
        domain = Gluing(id, [Domain_euclidean("sub" * id, dim, [RdPoint(zeros(dim)), RdPoint(ones(dim))])], (p1::GPoint, p2::GPoint) -> [])

    elseif id == "planeline"

        #   _________
        #  |         |
        #  |         |_________
        #  |      1  |       2
        #  |_________|
        #  

        domain = Gluing(id, [
            Domain_euclidean("plane", 2, [RdPoint([0.0,0.0]),RdPoint([4.0,4.0])]), 
            Domain_euclidean("line", 1, [RdPoint([0.0]),RdPoint([4.0])]), 
        ], (p1::GPoint, p2::GPoint) -> 
            if (p1.comp == p2.comp) && (norm(p1.p .- p2.p) <= 1e-7)
                if p1.comp == 1 && p1.p ≈ [4.0,2.0]
                    return [p1, GPoint(2,[0.0])]
                elseif p1.comp == 2 && p1.p ≈ [0.0]
                    return [p1, GPoint(1,[4.0,2.0])]
                else
                    return []
                end
            
            elseif p1.comp < p2.comp 
                return [GPoint(1,[4.0,2.0]), GPoint(2,[0.0])]

            elseif p1.comp > p2.comp
                return [GPoint(2,[0.0]), GPoint(1,[4.0,2.0])]

            else
                return []
            end
        )

    elseif id == "ISS"

        #    _____    _____                 _____    _____
        #   |     |  |     |    _______    |     |  |     |
        #   |     |  |     |   |       |   |     |  |     |
        #   |     |  |     |   |       |   |     |  |     | 
        #   |     |  |     |   |       |   |     |  |     |
        #   |     |  |     |   |    10 |   |     |  |     |
        #   |     |  |     |   |___ ___|   |     |  |     |
        #   |     |  |     |       |8      |     |  |     |
        #   |     |__|     |_______|_______|     |__|     |
        #   |     | 6|     |  _____|_____ 5|     | 7|     |
        #   |     |  |     | |___________| |     |  |     |
        #   |     |  |     |            9  |     |  |     |
        #   |     |  |     |               |     |  |     |
        #   |     |  |     |               |     |  |     |
        #   |     |  |     |               |     |  |     |
        #   |  1  |  |  2  |               |  3  |  |  4  |
        #   |_____|  |_____|               |_____|  |_____|
        # 

        solarpanels = [Domain_euclidean("panel$i", 2, [RdPoint([-0.5, -3.0]), RdPoint([0.5, 3.0])]) for i in 1:4] 
        unity = Domain_euclidean("unity", 1, [RdPoint([-1.5]), RdPoint([1.5])])
        S5 = Domain_euclidean("S5", 1, [RdPoint([0.0]), RdPoint([0.5])])
        P5 = Domain_euclidean("P5", 1, [RdPoint([0.0]), RdPoint([0.5])])
        S0 = Domain_euclidean("S0", 1, [RdPoint([-0.5]), RdPoint([1.0])])
        harmony = Domain_euclidean("harmony", 2, [RdPoint([-1.0,0.0]), RdPoint([1.0,0.5])])
        zarya = Domain_euclidean("zarya", 3, [RdPoint([-0.5,0.0,-0.05]), RdPoint([0.5,1.5,0.05])])

        function get_crossing_ISS(p1::GPoint, p2::GPoint)::Vector{GPoint}

            if (p1.comp == p2.comp) && (norm(p1.p .- p2.p) <= 1e-7) # to treat separately

                if (p1.comp == 1 && p1.p ≈ [0.5,0.0]) || (p1.comp == 6 && p1.p ≈ [0.0])
                    return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0])]
                elseif (p1.comp == 6 && p1.p ≈ [0.5]) || (p1.comp == 2 && p1.p ≈ [-0.5,0.0])
                    return [GPoint(6,[0.5]), GPoint(2,[-0.5,0.0])]
                elseif (p1.comp == 2 && p1.p ≈ [0.5,0.0]) || (p1.p == 5 && p1.p ≈ [-1.5])
                    return [GPoint(2,[0.5,0.0]), GPoint(5,[-1.5])]
                elseif (p1.comp == 5 && p1.p ≈ [1.5]) || (p1.comp == 3 && p1.p ≈ [-0.5,0.0])
                    return [GPoint(5,[1.5]), GPoint(3,[-0.5,0.0])]
                elseif (p1.comp == 3 && p1.p ≈ [0.5,0.0]) || (p1.comp == 7 && p1.p ≈ [0.0])
                    return [GPoint(3,[0.5,0.0]), GPoint(7,[0.0])]
                elseif (p1.comp == 7 && p1.p ≈ [0.5]) || (p1.comp == 4 && p1.p ≈ [-0.5,0.0])
                    return [GPoint(7, [0.5]), GPoint(4,[-0.5,0.0])]
                elseif (p1.comp == 5 && p1.p ≈ [0.0]) || (p1.comp == 8 && p1.p ≈ [0.0])
                    return [GPoint(5,[0.0]), GPoint(8,[0.0])]
                elseif (p1.comp == 8 && p1.p ≈ [-0.5]) || (p1.comp == 9 && p1.p ≈ [0.0,0.5])
                    return [GPoint(8,[-0.5]), GPoint(9,[0.0,0.5])]
                elseif (p1.comp == 8 && p1.p ≈ [1.0]) || (p1.comp == 10 && p1.p ≈ [0.0,0.0,0.0])
                    return [GPoint(8,[1.0]), GPoint(10,[0.0,0.0,0.0])]
                else 
                    return [] # not an intersection point
                end


            elseif p1.comp == p2.comp 
                return Vector{GPoint}()

            elseif p1.comp < p2.comp 

                if p1.comp == 1

                    if p2.comp == 2
                        return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0]), GPoint(6,[0.5]), GPoint(2,[-0.5,0.0])]
                    elseif p2.comp == 3
                        return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0]), GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), 
                                GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[1.5]), GPoint(3,[-0.5,0.0])]
                    elseif p2.comp == 4
                        return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0]), GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), 
                                GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[1.5]), GPoint(3,[-0.5,0.0]), 
                                GPoint(3,[0.5,0.0]), GPoint(7,[0.0]), GPoint(7,[0.5]), GPoint(4,[-0.5,0.0])]
                    elseif p2.comp == 5
                        return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0]), GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), 
                                GPoint(2,[0.5,0.0]), GPoint(5,[-1.5])]
                    elseif p2.comp == 6
                        return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0])]
                    elseif p2.comp == 7
                        return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0]), GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), 
                                GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[1.5]), GPoint(3,[-0.5,0.0]), 
                                GPoint(3,[0.5,0.0]), GPoint(7,[0.0])]
                    elseif p2.comp == 8
                        return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0]), GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), 
                                GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[0.0]), GPoint(8,[0.0])]
                    elseif p2.comp == 9
                        return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0]), GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), 
                                GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[0.0]), GPoint(8,[0.0]), 
                                GPoint(8,[-0.5]), GPoint(9,[0.0,0.5])]
                    elseif p2.comp == 10
                        return [GPoint(1,[0.5,0.0]), GPoint(6,[0.0]), GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), 
                                GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[0.0]), GPoint(8,[0.0]), 
                                GPoint(8,[1.0]), GPoint(10, [0.0,0.0,0.0])]
                    else 
                        throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                    end 

                elseif p1.comp == 2

                    if p2.comp == 3
                        return [GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[1.5]), GPoint(3,[-0.5,0.0])]
                    elseif p2.comp == 4
                        return [GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[1.5]), GPoint(3,[-0.5,0.0]), 
                                GPoint(3,[0.5,0.0]), GPoint(7,[0.0]), GPoint(7,[0.5]), GPoint(4,[-0.5,0.0])
                                ]
                    elseif p2.comp == 5
                        return [GPoint(2,[0.5,0.0]), GPoint(5,[-1.5])]
                    elseif p2.comp == 6
                        return [GPoint(2,[-0.5,0.0]), GPoint(6,[0.5])]
                    elseif p2.comp == 7
                        return [GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[1.5]), GPoint(3,[-0.5,0.0]), 
                                GPoint(3,[0.5,0.0]), GPoint(7,[0.0])]
                    elseif p2.comp == 8
                        return [GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[0.0]), GPoint(8,[0.0])]
                    elseif p2.comp == 9
                        return [GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[0.0]), GPoint(8,[0.0]), 
                                GPoint(8,[-0.5]), GPoint(9,[0.0,0.5])]
                    elseif p2.comp == 10
                        return [GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), GPoint(5,[0.0]), GPoint(8,[0.0]), 
                                GPoint(8,[1.0]), GPoint(10, [0.0,0.0,0.0])]
                    else 
                        throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                    end 

                elseif p1.comp == 3

                    if p2.comp == 4
                        return [GPoint(3,[0.5,0.0]), GPoint(7,[0.0]), GPoint(7,[0.5]), GPoint(4,[-0.5,0.0])
                                ]
                    elseif p2.comp == 5
                        return [GPoint(3,[-0.5,0.0]), GPoint(5,[1.5])]
                    elseif p2.comp == 6
                        return [GPoint(3,[-0.5,0.0]), GPoint(5,[1.5]), GPoint(5,[-1.5]), GPoint(2,[0.5,0.0]), 
                                GPoint(2,[-0.5,0.0]), GPoint(6,[0.5])]
                    elseif p2.comp == 7
                        return [GPoint(3,[0.5,0.0]), GPoint(7,[0.0])]
                    elseif p2.comp == 8
                        return [GPoint(3,[-0.5,0.0]), GPoint(5,[1.5]), GPoint(5,[0.0]), GPoint(8,[0.0])]
                    elseif p2.comp == 9
                        return [GPoint(3,[-0.5,0.0]), GPoint(5,[1.5]), GPoint(5,[0.0]), GPoint(8,[0.0]), 
                                GPoint(8,[-0.5]), GPoint(9,[0.0,0.5])]
                    elseif p2.comp == 10
                        return [GPoint(3,[-0.5,0.0]), GPoint(5,[1.5]), GPoint(5,[0.0]), GPoint(8,[0.0]), 
                                GPoint(8,[1.0]), GPoint(10, [0.0,0.0,0.0])]
                    else 
                        throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                    end 

                elseif p1.comp == 4

                    if p2.comp == 5
                        return [GPoint(4,[-0.5,0.0]), GPoint(7,[0.5]), GPoint(7,[0.0]), GPoint(3,[0.5,0.0]), 
                                GPoint(3, [-0.5,0.0]), GPoint(5,[1.5])]
                    elseif p2.comp == 6
                        return [GPoint(4,[-0.5,0.0]), GPoint(7,[0.5]), GPoint(7,[0.0]), GPoint(3,[0.5,0.0]), 
                                GPoint(3,[-0.5,0.0]), GPoint(5,[1.5]), GPoint(5,[-1.5]), GPoint(2,[0.5,0.0]), 
                                GPoint(2,[-0.5,0.0]), GPoint(6,[0.5])]
                    elseif p2.comp == 7
                        return [GPoint(4,[-0.5,0.0]), GPoint(7,[0.5])]
                    elseif p2.comp == 8
                        return [GPoint(4,[-0.5,0.0]), GPoint(7,[0.5]), GPoint(7,[0.0]), GPoint(3,[0.5,0.0]), 
                                GPoint(3, [-0.5,0.0]), GPoint(5,[1.5]), GPoint(5,[0.0]), GPoint(8,[0.0])]
                    elseif p2.comp == 9
                        return [GPoint(4,[-0.5,0.0]), GPoint(7,[0.5]), GPoint(7,[0.0]), GPoint(3,[0.5,0.0]), 
                                GPoint(3, [-0.5,0.0]), GPoint(5,[1.5]), GPoint(5,[0.0]), GPoint(8,[0.0]), 
                                GPoint(8,[-0.5]), GPoint(9,[0.0,0.5])]
                    elseif p2.comp == 10
                        return [GPoint(4,[-0.5,0.0]), GPoint(7,[0.5]), GPoint(7,[0.0]), GPoint(3,[0.5,0.0]), 
                                GPoint(3, [-0.5,0.0]), GPoint(5,[1.5]), GPoint(5,[0.0]), GPoint(8,[0.0]), 
                                GPoint(8,[1.0]), GPoint(10, [0.0,0.0,0.0])]
                    else 
                        throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                    end 

                elseif p1.comp == 5

                    if p2.comp == 6
                        return [GPoint(5,[-1.5]), GPoint(2,[0.5,0.0]), GPoint(2,[-0.5,0.0]), GPoint(6,[0.5])]
                    elseif p2.comp == 7
                        return [GPoint(5,[1.5]), GPoint(3,[-0.5,0.0]), GPoint(3, [0.5,0.0]), GPoint(7,[0.0])]
                    elseif p2.comp == 8
                        return [GPoint(5,[0.0]), GPoint(8,[0.0])]
                    elseif p2.comp == 9
                        return [GPoint(5,[0.0]), GPoint(8,[0.0]), GPoint(8,[-0.5]), GPoint(9,[0.0,0.5])]
                    elseif p2.comp == 10
                        return [GPoint(5,[0.0]), GPoint(8,[0.0]), GPoint(8,[1.0]), GPoint(10, [0.0,0.0,0.0])]
                    else 
                        throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                    end 

                elseif p1.comp == 6

                    if p2.comp == 7
                        return [GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), 
                                GPoint(5,[1.5]), GPoint(3,[-0.5,0.0]), GPoint(3,[0.5,0.0]), GPoint(7,[0.0])]
                    elseif p2.comp == 8
                        return [GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), 
                                GPoint(5,[0.0]), GPoint(8,[0.0])]
                    elseif p2.comp == 9
                        return [GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), 
                                GPoint(5,[0.0]), GPoint(8,[0.0]), GPoint(8,[-0.5]), GPoint(9,[0.0,0.5])]
                    elseif p2.comp == 10
                        return [GPoint(6,[0.5]), GPoint(2,[-0.5,0.0]), GPoint(2,[0.5,0.0]), GPoint(5,[-1.5]), 
                                GPoint(5,[0.0]), GPoint(8,[0.0]), GPoint(8,[1.0]), GPoint(10, [0.0,0.0,0.0])]
                    else 
                        throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                    end

                elseif p1.comp == 7

                    if p2.comp == 8
                        return [GPoint(7,[0.0]), GPoint(3,[0.5,0.0]), GPoint(3,[-0.5,0.0]), GPoint(5,[1.5]), 
                                GPoint(5,[0.0]), GPoint(8,[0.0])]
                    elseif p2.comp == 9
                        return [GPoint(7,[0.0]), GPoint(3,[0.5,0.0]), GPoint(3,[-0.5,0.0]), GPoint(5,[1.5]), 
                                GPoint(5,[0.0]), GPoint(8,[0.0]), GPoint(8,[-0.5]), GPoint(9,[0.0,0.5])]
                    elseif p2.comp == 10
                        return [GPoint(7,[0.0]), GPoint(3,[0.5,0.0]), GPoint(3,[-0.5,0.0]), GPoint(5,[1.5]), 
                                GPoint(5,[0.0]), GPoint(8,[0.0]), GPoint(8,[1.0]), GPoint(10, [0.0,0.0,0.0])]
                    else 
                        throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                    end

                elseif p1.comp == 8

                    if p2.comp == 9
                        return [GPoint(8,[-0.5]), GPoint(9,[0.0,0.5])]
                    elseif p2.comp == 10
                        return [GPoint(8,[1.0]), GPoint(10, [0.0,0.0,0.0])]
                    else 
                        throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                    end

                elseif p1.comp == 9

                    if p2.comp == 10
                        return [GPoint(9,[0.0,0.5]), GPoint(8,[-0.5]), GPoint(8,[1.0]), GPoint(10, [0.0,0.0,0.0])]
                    else 
                        throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                    end

                else
                    throw(ErrorException("component pair ($(p1.comp),$(p2.comp)) unknown"))
                end

            else # p1.comp > p2.comp 
                return get_crossing_ISS(p2, p1)[end:-1:1]
            end
        end # end get_crossing

        domain = Gluing(id, solarpanels ∪ [unity, S5, P5, S0, harmony, zarya], get_crossing_ISS)

    else
        throw(ErrorException("unknown domain id '$id'"))
    end

    verb(verbose, "... done.")
    return domain
end

export get_gluing_plotter 
function get_gluing_plotter(id::String, verbose::Bool=false)
    verb(verbose, "Creating plotter for $id...")

    if id == "1dsegment"
        
        plotter = Gluing_plotter(2, (domain::Gluing, comp::Int, p::EuPoint) -> [p.p[1],0.0])

    elseif id == "2dsquare"
        
        plotter = Gluing_plotter(2, (domain::Gluing, comp::Int, p::EuPoint) -> p.p)
        
    elseif id == "3dcube"

        plotter = Gluing_plotter(3, (domain::Gluing, comp::Int, p::EuPoint) -> p.p)

    elseif id == "planeline"

        plotter = Gluing_plotter(2, (domain::Gluing, comp::Int, p::EuPoint) ->
            if comp == 1
                return p.p 
            else
                return [4.0 + p.p[1], 2.0]
            end
        )

    elseif id == "ISS"

        plotter = Gluing_plotter(3, (domain::Gluing, comp::Int, p::EuPoint) ->
            if comp == 1
                return [p.p[1]-3.5,p.p[2],0.0]
            elseif comp == 2
                return [p.p[1]-2.0,p.p[2],0.0]
            elseif comp == 3
                return [p.p[1]+2.0,p.p[2],0.0]
            elseif comp == 4
                return [p.p[1]+3.5,p.p[2],0.0]
            elseif comp == 5
                return [p.p[1],0.0,0.0]
            elseif comp == 6
                return [p.p[1]-3.0,0.0,0.0]
            elseif comp == 7
                return [p.p[1]+2.5,0.0,0.0]
            elseif comp == 8
                return [0.0,p.p[1],0.0]
            elseif comp == 9
                return [p.p[1],p.p[2]-1.0,0.0]
            else
                return [p.p[1],p.p[2]+1.0,p.p[3]]
            end
        )

    end

    verb(verbose, "... done.")
    return plotter
end

export get_Dynamic
function get_Dynamic(id::String, domain::Gluing, verbose::Bool=false)
    verb(verbose, "Creating dynamic $id...")

    if (id == "eikonal" && domain.id == "1dsegment")

        h = 10.0
        return MathFunction{Dynamic}(id, x::GPoint -> [Velocity(GPoint(1,[-h]), 1.0, 1.0), Velocity(GPoint(1,[1.0+h]), 1.0, 1.0)])

    elseif (id == "eikonal" && domain.id == "2dsquare")

        h = 10.0
        return MathFunction{Dynamic}(id, x::GPoint -> [
            Velocity(GPoint(1,[x.p[1]-h,x.p[2]]), 1.0, 1.0),
            Velocity(GPoint(1,[x.p[1]+h,x.p[2]]), 1.0, 1.0),
            Velocity(GPoint(1,[x.p[1],x.p[2]-h]), 1.0, 1.0),
            Velocity(GPoint(1,[x.p[1],x.p[2]+h]), 1.0, 1.0),
        ])

    elseif (id == "eikonal" && domain.id == "3dcube")

        h = 10.0
        return MathFunction{Dynamic}(id, x::GPoint -> [
            Velocity(GPoint(1,[x.p[1]-h,x.p[2],x.p[3]]), 1.0, 1.0),
            Velocity(GPoint(1,[x.p[1]+h,x.p[2],x.p[3]]), 1.0, 1.0),
            Velocity(GPoint(1,[x.p[1],x.p[2]-h,x.p[3]]), 1.0, 1.0),
            Velocity(GPoint(1,[x.p[1],x.p[2]+h,x.p[3]]), 1.0, 1.0),
            Velocity(GPoint(1,[x.p[1],x.p[2],x.p[3]-h]), 1.0, 1.0),
            Velocity(GPoint(1,[x.p[1],x.p[2],x.p[3]+h]), 1.0, 1.0),
        ])

    elseif (id == "eikonal" && domain.id == "planeline")

        h = 10.0
        return MathFunction{Dynamic}(id, x::GPoint -> 
            if x.comp == 1
                if sum(abs.(x.p .- [4.0,2.0])) <= 1e-7
                    mag = 1.0
                else
                    mag = norm(x.p .- [4.0,2.0]) / sum(abs.(x.p .- [4.0,2.0]))
                end
                return [
                    Velocity(GPoint(1,[x.p[1]-h,x.p[2]]), 1.0, 1.0),
                    Velocity(GPoint(1,[x.p[1]+h,x.p[2]]), 1.0, 1.0),
                    Velocity(GPoint(1,[x.p[1],x.p[2]-h]), 1.0, 1.0),
                    Velocity(GPoint(1,[x.p[1],x.p[2]+h]), 1.0, 1.0),
                    Velocity(GPoint(2,[h]),mag, 1.0)
                ]
            else 
                return [
                    Velocity(GPoint(1,[-h,0.5]), 1.0, 1.0), 
                    Velocity(GPoint(1,[1.0, h]), 1.0, 1.0),
                    Velocity(GPoint(1,[1.0,-h]), 1.0, 1.0),
                    Velocity(GPoint(2,[h]),1.0, 1.0)
                ]
            end
        )

    elseif (id == "windrobot" && domain.id == "planeline")

        h = 10.0
        wind = 2.0
        return MathFunction{Dynamic}(id, x::GPoint -> 
            if x.comp == 1
                if sum(abs.(x.p .- [4.0,2.0])) <= 1e-7 # avoid dangerous divisions 
                    mag = wind
                else
                    mag = norm(x.p .- [4.0,2.0]) / ((4.0 - x.p[1])/wind + abs(x.p[2] - 2.0))
                end
                return [
                    Velocity(GPoint(1,[0.0,x.p[2]]), 1.0, 1.0),     # left
                    Velocity(GPoint(1,[4.0,x.p[2]]), wind, 1.0),    # right
                    Velocity(GPoint(1,[x.p[1],0.0]), 1.0, 1.0),     # down
                    Velocity(GPoint(1,[x.p[1],4.0]), 1.0, 1.0),     # up 
                    Velocity(GPoint(2,[4.0]),mag, 1.0)              # jump on the segment 
                ]

            else 
                return [
                    Velocity(GPoint(1,[0.0,2.0]), 1.0, 1.0), 
                    Velocity(GPoint(2,[4.0]),wind, 1.0)
                ]
            end
        )

    elseif (id == "eikonal" && domain.id == "ISS")

        h = 10.0

        res = MathFunction{Dynamic}(id, x::GPoint -> 
            if x.comp ∈ [1,2,3,4]
                return [
                    Velocity(GPoint(x.comp,[x.p[1],x.p[2]+h]),1.0,1.0),      # up
                    Velocity(GPoint(x.comp,[x.p[1],x.p[2]-h]),1.0,1.0),      # down 
                    Velocity(GPoint(x.comp,[x.p[1]-h,x.p[2]]),1.0,1.0),      # left 
                    Velocity(GPoint(x.comp,[x.p[1],x.p[2]+h]),1.0,1.0),      # right 
                ] ∪ [
                    Velocity(GPoint(4,[h,0.0]),1.0,1.0) for _ in [1] if x.comp in [1,2,3]   # left jump
                ] ∪ [
                    Velocity(GPoint(1,[-h,0.0]),1.0,1.0) for _ in [1] if x.comp in [2,3,4]  # right jump
                ] 
                
            elseif x.comp ∈ [5,6,7]
                return [Velocity(GPoint(1,[-h,0.0]),1.0,1.0), Velocity(GPoint(4,[h,0.0]),1.0,1.0)] 
                
            elseif x.comp == 8
                return [Velocity(GPoint(9,[0.0,-h]),1.0,1.0), Velocity(GPoint(10,[0.0,h,0.0]),1.0,1.0)]
            elseif x.comp == 9
                return [
                    Velocity(GPoint(9,[-h,0.0]),1.0,1.0),       # up
                    Velocity(GPoint(9,[ h,0.0]),1.0,1.0),       # down 
                    Velocity(GPoint(9,[0.0,-h]),1.0,1.0),       # left 
                    Velocity(GPoint(9,[0.0, h]),1.0,1.0),       # right 
                    Velocity(GPoint(10,[0.0,h,0.0]),1.0,1.0)    # jump up
                ]
            else 
                return [
                    Velocity(GPoint(9,[0.0,-h]),1.0,1.0)        # jump down
                ] ∪ [
                    Velocity(GPoint(10,x.p .+ a*b),1.0,1.0) for a ∈ [-1.0,1.0] for b in [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] # sides 
                ]
            end
        )

    elseif (id == "windrobot" && domain.id == "ISS")

        wind = 2.0
        thresh = 0.25

        euclideannorm(v1::Float64, v2::Float64) = sqrt(v1^2 + v2^2)
        distordednorm(v1::Float64, v2::Float64) = max(v1,0.0)/wind + max(-v1,0.0) + abs(v2)

        res = MathFunction{Dynamic}(id, x::GPoint -> 
            if x.comp ∈ [1,2,3,4]
                if abs(x.p[2]) <= 1e-7 # avoid dangerous divisions
                    return [
                        Velocity(GPoint(x.comp,[x.p[1], 3.0]),1.0,thresh),       # up 
                        Velocity(GPoint(x.comp,[x.p[1],-3.0]),1.0,thresh),       # down 
                        Velocity(GPoint(1,[-0.5,0.0]),1.0,thresh),               # left
                        Velocity(GPoint(4,[ 0.5,0.0]),wind,thresh)               # right
                    ]
                else
                    locres = [
                        Velocity(GPoint(x.comp,[x.p[1], 3.0]),1.0, thresh),    # up
                        Velocity(GPoint(x.comp,[x.p[1],-3.0]),1.0, thresh),    # down 
                        Velocity(GPoint(x.comp,[-0.5,x.p[2]]),1.0, thresh),    # left
                        Velocity(GPoint(x.comp,[ 0.5,x.p[2]]),wind,thresh),    # right 
                    ]
                    if x.comp ∈ [2,3,4]
                        locres = vcat(locres, [
                            Velocity(GPoint(1,[-0.5,0.0]), euclideannorm(-0.5-x.p[1],-x.p[2]) / (abs(-0.5-x.p[1])+abs(x.p[2])), thresh) # jump left
                        ])
                    end
                    if x.comp ∈ [1,2,3]
                        locres = vcat(locres, [
                            Velocity(GPoint(4,[ 0.5,0.0]), euclideannorm(0.5-x.p[1],-x.p[2]) / distordednorm(0.5-x.p[1],-x.p[2]), thresh) # jump right
                        ])
                    end
                    return locres
                end
            elseif x.comp ∈ [5,6,7]
                return [Velocity(GPoint(1,[-0.5,0.0]),1.0,thresh), Velocity(GPoint(4,[0.5,0.0]),wind, thresh)]
            elseif x.comp == 8
                return [Velocity(GPoint(1,[-0.5,0.0]),1.0,thresh), Velocity(GPoint(4,[0.5,0.0]),max(1.0,wind*(1.0-abs(x.p[1]))), thresh)]
            else # x.comp ∈ [9,10]
                return [Velocity(GPoint(1,[-0.5,0.0]),1.0,thresh)] # arbitrary right-left point
            end
        )

    else
        throw(ErrorException("unknown dynamic id '$id'"))
    end

    verb(verbose, "... done.")
    return res
end

export get_Scalar
function get_Scalar(id::String, domain::Gluing, verbose::Bool=false)
    verb(verbose, "Creating scalar function $id...")

    if (id == "norm" && domain.id ∈ ["1dsegment","2dsquare","3dcube"])

        res = MathFunction{Scalar}(id, x::GPoint -> norm(x.p))

    elseif (id == "norm" && domain.id == "planeline")

        res = MathFunction{Scalar}(id, x::GPoint -> 
            if x.comp == 1
                return norm(x.p)
            else
                return sqrt(4.0^2 + 2.0^2) + x.p[1]        
            end
        )

    elseif (id == "normend" && domain.id == "planeline")

        res = MathFunction{Scalar}(id, x::GPoint -> 
            if x.comp == 1
                return 4.0 + norm(x.p .- [4.0,2.0])
            else
                return 4.0 - x.p[1]
            end
        )

    elseif (id == "distmeteors" && domain.id == "planeline")

        res = MathFunction{Scalar}(id, x::GPoint -> min(
            get_distance(domain, GPoint(1,[1.0,3.0]), x), 
            get_distance(domain, GPoint(2,[3.0]), x)
        ))

    elseif (id == "norm" && domain.id == "ISS")

        res = MathFunction{Scalar}(id, x::GPoint -> get_distance(domain, x, GPoint(5,[0.0])))

    elseif (id == "distmeteors" && domain.id == "ISS")

        res = MathFunction{Scalar}(id, x::GPoint -> minimum([get_distance(domain, x, GPoint(c,[a,b])) for (c,a,b) in [(1,0.0,-2.0),(2,0.0,2.0),(3,0.0,-1.0)]]))
    
    else
        throw(ErrorException("unknown scalar id '$id'"))
    end

    verb(verbose, "... done.")
    return res
end

export get_Analytical
"""
    Return the exact solution as a MathFunction{Scalar}.
    When probe = true, only return the existence of a solution.
"""
function get_Analytical(domain::Gluing, dynamicid::String, costid::String, T::Float64, probe::Bool=false, verbose::Bool=false)
    verb(verbose, "Computing analytical V for domain $(domain.id), dynamic $dynamicid and cost $costid... ")

    id = domain.id * "_" * dynamicid * "_" * costid

    if (dynamicid == "eikonal" && costid == "norm" && domain.id ∈ ["1dsegment","2dsquare","3dcube"])

        if probe 
            return true 
        else
            if domain.comp[1].dim == 1
                res = MathFunction{Scalar}(id, (t::Float64, x::GPoint) -> max(0.0, abs(x.p[1]) - (T-t)))
            elseif domain.comp[1].dim == 2
                res = MathFunction{Scalar}(id, (t::Float64, x::GPoint) -> 
                    if t >= T - 1e-7
                        return norm(x.p)
                    elseif sum(abs.(x.p)) <= T-t + 1e-7
                        return 0.0
                    else
                        point = [0.0,0.0]
                        aa = x.p .- [T-t,0.0]
                        bb = x.p .- [0.0,T-t]
                        if ((point .- aa)' * (bb .- aa) <= 0) || ((point .- bb)' * (aa .- bb) <= 0)
                            # one of two extremals 
                            return min(sqrt((point[1]-aa[1])^2+(point[2]-aa[2])^2), sqrt((point[1]-bb[1])^2+(point[2]-bb[2])^2))
                        else
                            # projection on the line 
                            return sqrt((point[1]-aa[1])^2+(point[2]-aa[2])^2 - ((point.-aa)' * (bb .- aa))^2 / ((aa[1]-bb[1])^2 + (aa[2]-bb[2])^2))
                        end
                        # return max(0.0, norm(x.p)*(1.0 - (T-t)/sum(abs.(x.p))))
                    end
                )
            else
                # not tested 
                res = MathFunction{Scalar}(id, (t::Float64, x::GPoint) -> 
                    if t >= T - 1e-7
                        return norm(x.p)
                    elseif sum(abs.(x.p)) <= T-t + 1e-7
                        return 0.0
                    else
                        """
                            p1 = x - [T-t,0,0], p2 = x - [0,T-t,0], p3 = x - [0,0,T-t]
                            α [1,1,1] ∈ Vect(triangle)
                            λ p1 + μ p2 + (1-λ-μ) p3 = -α [1,1,1]

                            Sachant que p1-p3 = (T-t) [-1,0,1] et p2-p3 = (T-t) [0,-1,1],
                            on a 
                                q1 ≔ [-1,0,1] ∧ [1,1,1] = [-1,2,-1], 
                                q2 ≔ [0,-1,1] ∧ [1,1,1] = [-2,1,1],
                            et 

                            λ (p1-p3) ⋅ q1 + μ (p2-p3) ⋅ q1 = - p3 ⋅ q1 = 0 + μ (T-t) (-3) 
                            λ (p1-p3) ⋅ q2 + μ (p2-p3) ⋅ q2 = - p3 ⋅ q2 = λ (T-t) 3 

                            donc 
                                λ = - p3 ⋅ q2 / ( 3(T-t)) 
                                μ = - p3 ⋅ q1 / (-3(T-t))
                                proj = λ (x-a) + μ (x-b) + (1-λ-μ) (x-c) = - (λ a + μ b + (1-λ-μ) c)
                                     = 1/3 (-p3⋅q2 [1,0,0] + p3⋅q1 [0,1,0] + (1+p3⋅q2-p3⋅q1) [0,0,1])
                                     = [p3 ⋅ q2 / 3, p3 ⋅ q1 / 3, (1+p3⋅q2-p3⋅q1) / 3]
                                     = - (T-t) [λ, μ, 1-λ-μ]
                        """
                        λ = - (x.p .- [0.0,0.0,T-t])' * [-2,1,1 ] / ( 3.0 * (T-t))
                        μ = - (x.p .- [0.0,0.0,T-t])' * [-1,2,-1] / (-3.0 * (T-t))
                        if λ >= 0.0 && μ >= 0.0 && λ + μ <= 1.0
                            # barycentric coordinates : we are in the interior, 
                            # the distance is the norm of the projection 
                            return sqrt(λ^2 + μ^2 + (1-λ-μ)^2) * (T-t)
                        else
                            # minimum reached at one of the edges 
                            return sqrt(minimum([(x.p[1]-(T-t))^2+x.p[2]^2+x.p[3]^2, x.p[1]^2+(x.p[2]-(T-t))^2+x.p[3]^2, x.p[1]^2+x.p[2]^2+(x.p[3]-(T-t))^2]))
                        end
                    end
                )
            end
        end

    elseif (dynamicid == "eikonal" && costid == "norm" && domain.id == "planeline")

        if probe 
            return true 
        else
            res = MathFunction{Scalar}(id, (t::Float64, x::GPoint) -> 
                if t >= T - 1e-7
                    if x.comp == 1
                        return norm(x.p)
                    else
                        return sqrt(4.0^2 + 2.0^2) + x.p[1]
                    end
                elseif x.comp == 1
                    if abs(x.p[1] - x.p[2]) > T - t # not reaching the diagonal before the end of time
                        if x.p[2] > x.p[1]
                            return sqrt(x.p[1]^2 + (x.p[2] - (T-t))^2)
                        else
                            return sqrt((x.p[1]-(T-t))^2 + x.p[2]^2)
                        end
                    else # reaching the diagonal 
                        return max(0.0, sqrt(2) * minimum(x.p) - (T-t-abs(x.p[1] - x.p[2]))/sqrt(2))
                    end
                else # x.comp = 2
                    if x.p[1] >= T-t # stay within component 2
                        return sqrt(4.0^2 + 2.0^2) + x.p[1] - (T-t)
                    else
                        if 2.0 >= T-t-x.p[1] # not reaching the diagonal before the end of time  
                            return sqrt((4.0 - (T-t-x.p[1]))^2 + 2.0^2)
                        else  # reaching the diagonal 
                            return max(0.0, sqrt(2.0^2 + 2.0^2) - (T-t-x.p[1]-2.0)/sqrt(2)) # times √2 because norm of a constant vector in dim 2 
                        end
                    end
                end
            )
        end

    elseif (dynamicid == "eikonal" && costid == "normend" && domain.id == "planeline")

        if probe 
            return true 
        else
            res = MathFunction{Scalar}(id, (t::Float64, x::GPoint) -> 
                if t >= T - 1e-7
                    if x.comp == 1
                        return 4.0 + norm(x.p .- [4.0,2.0])
                    else
                        return 4.0 - x.p[1]
                    end
                elseif x.comp == 1
                    #    __________________
                    #   |                  |
                    #   |               \  |
                    #   |                \ |
                    #   |                 \|____
                    #   |                 /|
                    #   |                / |
                    #   |               /  |
                    #   |__________________|
                    # 
                    if abs(abs(x.p[1]-4.0) - abs(x.p[2]-2.0)) > T-t # not reaching the diagonal 
                        # endpoint of the characteristic: x.p + (T-t) * (1.0*(|x.p[2]| < |x.p[1]|), -sign(x.p[2]-2.0)*(|x.p[1]| < |x.p[2]|))
                        return 1.0 + sqrt((x.p[1]+(T-t)*(abs(x.p[2]) < abs(x.p[1])) - 4.0)^2 + (x.p[2]-(T-t)*sign(x.p[2]-2.0)*(abs(x.p[2]) > abs(x.p[1])) - 2.0)^2)
                    elseif abs(abs(x.p[1]-4.0) - abs(x.p[2]-2.0)) + min(abs(x.p[1]-4.0),abs(x.p[2]-2.0))/sqrt(2) > T-t # not reaching the junction 
                        # endpoint of the characteristic: (-α,sign(x.p[2]-2.0)*α) + (T-t-||x.p[1]-1.0|-|x.p[2]-0.5||)/√2 * (target-proj)/|target-proj|, 
                        #     where α = min(|x.p[1]-4.0|,|x.p[2]-2.0|)
                        α = min(abs(x.p[1]-4.0),abs(x.p[2]-2.0))
                        proj = [-α, sign(x.p[2]-2.0)*α]
                        endpoint = proj + (T-t-abs(abs(x.p[1]-4.0) - abs(x.p[2]-2.0)))/√2 * [1.0, -sign(x.p[2]-2.0)]
                        return 4.0 + norm(endpoint .- [4.0,2.0])
                    else # reaching the junction 
                        α = min(abs(x.p[1]-4.0),abs(x.p[2]-2.0))
                        proj = [-α, sign(x.p[2]-2.0)*α]
                        return max(0.0, 4.0 - (T - (t + abs(abs(x.p[1]-4.0) - abs(x.p[2]-2.0)) + norm(proj .- [4.0,2.0]) * sqrt(2))))
                    end
                else # x.comp = 2
                    return max(0.0, 4.0 - x.p[1] - (T-t))
                end
            )
        end

    elseif (dynamicid == "windrobot" && costid == "distmeteors" && domain.id == "planeline")

        if probe 
            return true 
        else 

            x1 = [1.0,3.0]
            x2 = [3.0]
            wind = 2.0
            # returning endpoint, target_is_reached, surplus_of_time
            function move_in_planeline(x1::Float64, x2::Float64, tar1::Float64, tar2::Float64, allowedtime::Float64)

                if x1 >= tar1 # moving with speed 1, no wind in action 
                    if allowedtime <= abs(abs(x1-tar1) - abs(x2-tar2)) # not reaching the diagonal 
                        return [
                            x1 - allowedtime * (abs(x1-tar1) > abs(x2-tar2)), 
                            x2 - allowedtime * sign(x2-tar2)*(abs(x1-tar1) < abs(x2-tar2))
                        ], false, 0.0
                    # distance = √2 * min(abs(x1-tar1), abs(x2-tar2))
                    # speed = 1/√2
                    # time = distance / speed = 2 * min(abs(x1-tar1), abs(x2-tar2))
                    elseif allowedtime <= abs(abs(x1-tar1) - abs(x2-tar2)) + 2.0 * min(abs(x1-tar1), abs(x2-tar2)) # reaching the diagonal but not the junction 
                        timetodiag = abs(abs(x1-tar1) - abs(x2-tar2))
                        return [
                            x1 -  timetodiag * (abs(x1-tar1) > abs(x2-tar2)) - (allowedtime-timetodiag) / 2.0, 
                            x2 - (timetodiag * (abs(x1-tar1) < abs(x2-tar2)) + (allowedtime-timetodiag) / 2.0) * sign(x2-tar2)
                        ], false, 0.0
                    else # reaching the junction 
                        return [tar1, tar2], true, allowedtime - abs(abs(x1-tar1) - abs(x2-tar2)) - 2.0 * min(abs(x1-tar1), abs(x2-tar2))
                    end
                else # x1-tar1 < 0, wind in action 
                    # the action of the wind changes the inclination of the diagonal 
                    # now described by abs(x2-tar2) = wind * abs(x1-tar1)

                    if abs(x2 - tar2) <= wind * abs(x1 - tar1) # moving horizontally first
                        timetodiag = abs((tar1-x1)/wind - abs(x2 - tar2)/wind^2)
                        if allowedtime <= timetodiag
                            return [x1 + allowedtime * wind, x2], false, 0.0
                        end # implicit else 
                        projondiag = [x1 + timetodiag * wind, x2]
                    else # moving vertically first 
                        timetodiag = abs(x2 - tar2) - wind * abs(x1-tar1)
                        if allowedtime <= timetodiag
                            return [x1, x2 + sign(tar2-x2) * allowedtime], false, 0.0
                        end # implicit else 
                        projondiag = [x1, x2 + sign(tar2-x2) * timetodiag]
                    end
                    # so allowedtime > timetodiag 
                    speedondiag = wind / sqrt(1.0 + wind^2)
                    timetojunc = timetodiag + norm(projondiag .- [tar1,tar2]) / speedondiag
                    if allowedtime <= timetojunc
                        return [
                            projondiag[1] + (allowedtime-timetodiag) * speedondiag/sqrt(1+wind^2), 
                            projondiag[2] - (allowedtime-timetodiag) * speedondiag/sqrt(1+wind^2) * sign(x2-tar2) * wind
                        ], false, 0.0
                    else # reaching the junction 
                        return [tar1, tar2], true, allowedtime - timetojunc
                    end
                end
            end

            function V1pl(t::Float64, x::GPoint) # target C1=>[1,3]
                if t >= T - 1e-7
                    if x.comp == 1
                        return sqrt((x.p[1]-x1[1])^2 + (x.p[2]-x1[2])^2)
                    else 
                        return sqrt((4.0-x1[1])^2 + (2.0-x1[2])^2) + x.p[1]
                    end

                elseif x.comp == 1
                    endpoint, target_is_reached, surplus_of_time = move_in_planeline(x.p[1], x.p[2], x1[1], x1[2], T-t)
                    if target_is_reached
                        return 0.0
                    else
                        return sqrt((endpoint[1]-x1[1])^2 + (endpoint[2]-x1[2])^2)
                    end

                else 
                    if x.p[1] > T-t # we don't reach the panel 1 before the end of the time 
                        return sqrt((4.0-x1[1])^2 + (2.0-x1[2])^2) + x.p[1] - (T-t)
                    else
                        return V1pl(t + x.p[1], GPoint(1,[4.0,2.0]))
                    end
                end
            end # function V1

            function V2pl(t::Float64, x::GPoint) # target C2=>[3.0]
                if t >= T - 1e-7
                    if x.comp == 1
                        return x2[1] + sqrt((x.p[1]-4.0)^2 + (x.p[2]-2.0)^2)
                    else 
                        return abs(x2[1] - x.p[1])
                    end

                elseif x.comp == 1
                    endpoint, target_is_reached, surplus_of_time = move_in_planeline(x.p[1], x.p[2], 4.0, 2.0, T-t)
                    if target_is_reached
                        return V2pl(T-surplus_of_time, GPoint(2,[0.0]))
                    else
                        return x2[1] + sqrt((endpoint[1]-4.0)^2 + (endpoint[2]-2.0)^2)
                    end

                else 
                    if x.p[1] >= x2[1]
                        return max(0.0, abs(x2[1]-x.p[1]) - (T-t))
                    else
                        return max(0.0, abs(x2[1]-x.p[1]) - wind*(T-t))
                    end
                end
            end # function V2

            res = MathFunction{Scalar}(id, (t::Float64, x::GPoint) -> min(V1pl(t,x), V2pl(t,x)))

        end

    elseif (dynamicid == "eikonal" && costid == "norm" && domain.id == "ISS")

        if probe 
            return true 
        else
            res = MathFunction{Scalar}(id, (t::Float64, x::GPoint) -> max(0.0, get_distance(domain, GPoint(5,[0.0]), x) - (T-t)))
        end

    elseif (dynamicid == "windrobot" && costid == "distmeteors" && domain.id == "ISS")

        if probe 
            return true 
        else 

            wind = 2.0

            # returning endpoint, target_is_reached, surplus_of_time
            function move_in_panel(x1::Float64, x2::Float64, tar1::Float64, tar2::Float64, allowedtime::Float64)
                
                if x1 >= tar1 # moving with speed 1, no wind in action 
                    if allowedtime <= abs(abs(x1-tar1) - abs(x2-tar2)) # not reaching the diagonal 
                        return [
                            x1 - allowedtime * (abs(x1-tar1) > abs(x2-tar2)), 
                            x2 - allowedtime * sign(x2-tar2)*(abs(x1-tar1) < abs(x2-tar2))
                        ], false, 0.0
                    # distance = √2 * min(abs(x1-tar1), abs(x2-tar2))
                    # speed = 1/√2
                    # time = distance / speed = 2 * min(abs(x1-tar1), abs(x2-tar2))
                    elseif allowedtime <= abs(abs(x1-tar1) - abs(x2-tar2)) + 2.0 * min(abs(x1-tar1), abs(x2-tar2)) # reaching the diagonal but not the junction 
                        timetodiag = abs(abs(x1-tar1) - abs(x2-tar2))
                        return [
                            x1 -  timetodiag * (abs(x1-tar1) > abs(x2-tar2)) - (allowedtime-timetodiag) / 2.0, 
                            x2 - (timetodiag * (abs(x1-tar1) < abs(x2-tar2)) + (allowedtime-timetodiag) / 2.0) * sign(x2-tar2)
                        ], false, 0.0
                    else # reaching the junction 
                        return [tar1, tar2], true, allowedtime - abs(abs(x1-tar1) - abs(x2-tar2)) - 2.0 * min(abs(x1-tar1), abs(x2-tar2))
                    end
                else # x1-tar1 < 0, wind in action 
                    # the action of the wind changes the inclination of the diagonal 
                    # now described by abs(x2-tar2) = wind * abs(x1-tar1)

                    if abs(x2 - tar2) <= wind * abs(x1 - tar1) # moving horizontally first
                        timetodiag = abs((tar1-x1)/wind - abs(x2 - tar2)/wind^2)
                        if allowedtime <= timetodiag
                            return [x1 + allowedtime * wind, x2], false, 0.0
                        end # implicit else 
                        projondiag = [x1 + timetodiag * wind, x2]
                    else # moving vertically first 
                        timetodiag = abs(x2 - tar2) - wind * abs(x1-tar1)
                        if allowedtime <= timetodiag
                            return [x1, x2 + sign(tar2-x2) * allowedtime], false, 0.0
                        end # implicit else 
                        projondiag = [x1, x2 + sign(tar2-x2) * timetodiag]
                    end
                    # so allowedtime > timetodiag 
                    speedondiag = wind / sqrt(1.0 + wind^2)
                    timetojunc = timetodiag + norm(projondiag .- [tar1,tar2]) / speedondiag
                    if allowedtime <= timetojunc
                        return [
                            projondiag[1] + (allowedtime-timetodiag) * speedondiag/sqrt(1+wind^2), 
                            projondiag[2] - (allowedtime-timetodiag) * speedondiag/sqrt(1+wind^2) * sign(x2-tar2) * wind
                        ], false, 0.0
                    else # reaching the junction 
                        return [tar1, tar2], true, allowedtime - timetojunc
                    end
                end
            end

            # ̇yₜ = -2(1-yₜ) => yₜ = e^{2(t-s)} yₛ + e^{2(t-s)} ∫_τ=s^t e^{-2(τ-s)} (-2) dτ
            # = e^{2(t-s)} yₛ - 2 e^{2t} ∫_τ=s^t e^{-2τ} dτ
            # = e^{2(t-s)} yₛ + e^{2t} (e^{-2t} - e^{-2s})
            # = e^{2(t-s)} (yₛ-1) + 1.
            # So, starting from |x| < 1/2 at time t0, and letting y := |x|, the point 
            # 0 is reached when 
            # 0 = e^{2(t-t0)} (|x|-1) + 1
            # e^{2(t-t0)} = 1 / (1-|x|)
            # 2(t-t0) = log(1/(1-|x|)) = - log(1-|x|)
            # t = t0 - 1/2 log(1-|x|)

            # returning endpoint, junction_is_reached, surplus_of_time
            function move_in_unity(x::Float64, allowedtime::Float64)
                if x > 1/2
                    if allowedtime <= x-1/2
                        # moving at speed 1 before the acceleration zone 
                        return x-allowedtime, false, 0.0
                    else 
                        # jumping at the beginning of the acceleration zone 
                        allowedtime -= x-1/2
                        x = 1/2
                    end
                end
                # by now, |x| <= 1/2
                timetojunc = - 0.5 * log(1.0 - abs(x))
                if allowedtime >= timetojunc
                    return 0.0, true, allowedtime - timetojunc
                else
                    res = sign(x) * (exp(2*allowedtime)*(abs(x)-1.0) + 1.0)
                    return res, false, 0.0
                end
            end

            function V1(t::Float64, x::GPoint) # target C1=>[0.0,-2.0]
                if t >= T - 1e-7
                    if x.comp == 1
                        return sqrt(x.p[1]^2 + (x.p[2]-(-2.0))^2)
                    elseif x.comp == 6
                        return sqrt(2.0^2 + 0.5^2) + x.p[1]
                    elseif x.comp == 2
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + sqrt((x.p[1]-(-0.5))^2 + x.p[2]^2)
                    elseif x.comp == 5
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + x.p[1]-(-1.5)
                    elseif x.comp == 3
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 3.0 + sqrt((x.p[1]-(-0.5))^2 + x.p[2]^2)
                    elseif x.comp == 7
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 3.0 + 1.0 + x.p[1]
                    elseif x.comp == 4
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 3.0 + 1.0 + 0.5 + sqrt((x.p[1]-(-0.5))^2 + x.p[2]^2)
                    elseif x.comp == 8
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 1.5 + abs(x.p[1]) 
                    elseif x.comp == 9
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 1.5 + 0.5 + sqrt(x.p[1]^2 + (x.p[2]-0.5)^2)
                    else 
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 1.5 + 1.0 + norm(x.p)
                    end

                elseif x.comp == 1
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], 0.0, -2.0, T-t)
                    if target_is_reached
                        return 0.0
                    else
                        return sqrt(endpoint[1]^2 + (endpoint[2] - (-2.0))^2)
                    end

                elseif x.comp == 6
                    if x.p[1] > T-t # we don't reach the panel 1 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + x.p[1] - (T-t)
                    else
                        return V1(t + x.p[1], GPoint(1,[0.5,0.0]))
                    end

                elseif x.comp == 2
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], -0.5, 0.0, T-t)
                    if target_is_reached
                        return V1(T-surplus_of_time, GPoint(6,[0.5]))
                    else
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + sqrt((-0.5-endpoint[1])^2 + endpoint[2]^2)
                    end

                elseif x.comp == 5
                    if x.p[1]-(-1.5) > T-t # we don't reach the panel 2 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + (x.p[1]+1.5) - (T-t)
                    else
                        return V1(t + x.p[1] + 1.5, GPoint(2,[0.5,0.0]))
                    end

                elseif x.comp == 3
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], -0.5, 0.0, T-t)
                    if target_is_reached
                        return V1(T-surplus_of_time, GPoint(5,[1.5]))
                    else
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 3.0 + sqrt((-0.5-endpoint[1])^2 + endpoint[2]^2)
                    end

                elseif x.comp == 7
                    if x.p[1] > T-t # we don't reach the panel 3 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 3.0 + 1.0 + x.p[1] - (T-t)
                    else
                        return V1(t + x.p[1], GPoint(3,[0.5,0.0]))
                    end

                elseif x.comp == 4
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], -0.5, 0.0, T-t)
                    if target_is_reached
                        return V1(T-surplus_of_time, GPoint(7,[0.5]))
                    else
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 3.0 + 1.0 + 0.5 + sqrt((-0.5-endpoint[1])^2 + endpoint[2]^2)
                    end

                elseif x.comp == 8
                    endpoint, junction_is_reached, surplus_of_time = move_in_unity(x.p[1], T-t)
                    if junction_is_reached
                        return V1(T-surplus_of_time, GPoint(5,[0.0]))
                    else 
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 1.5 + abs(endpoint)
                    end

                elseif x.comp == 9 # Harmony: we are only allowed to move at speed one towards Unity (8)
                    if sqrt(x.p[1]^2 + (x.p[2]-0.5)^2) > T-t # we don't reach 8 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 1.5 + 0.5 + sqrt(x.p[1]^2 + (x.p[2]-0.5)^2) - (T-t)
                    else
                        return V1(t + sqrt(x.p[1]^2 + (x.p[2]-0.5)^2), GPoint(8,[-0.5]))
                    end

                else # x.comp = 10, Zarya: we are only allowed to move at speed one towards Unity (8)
                    if norm(x.p) > T-t # we don't reach 8 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + 1.0 + 1.5 + 1.0 + norm(x.p) - (T-t)
                    else
                        return V1(t + norm(x.p), GPoint(8,[1.0]))
                    end
                end # distinction x.comp
            end # function V1 

            function V2(t::Float64, x::GPoint) # target C2=>[0.0,2.0]
                if t >= T 
                    if x.comp == 2
                        return sqrt(x.p[1]^2 + (x.p[2]-2.0)^2)
                    elseif x.comp == 6
                        return sqrt(2.0^2 + 0.5^2) + (0.5 - x.p[1])
                    elseif x.comp == 1
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + sqrt((x.p[1]-0.5)^2 + x.p[2]^2)
                    elseif x.comp == 5
                        return sqrt(2.0^2 + 0.5^2) + x.p[1]-(-1.5)
                    elseif x.comp == 3
                        return sqrt(2.0^2 + 0.5^2) + 3.0 + sqrt((x.p[1]-(-0.5))^2 + x.p[2]^2)
                    elseif x.comp == 7
                        return sqrt(2.0^2 + 0.5^2) + 3.0 + 1.0 + x.p[1]
                    elseif x.comp == 4
                        return sqrt(2.0^2 + 0.5^2) + 3.0 + 1.0 + 0.5 + sqrt((x.p[1]-(-0.5))^2 + x.p[2]^2)
                    elseif x.comp == 8
                        return sqrt(2.0^2 + 0.5^2) + 1.5 + abs(x.p[1])
                    elseif x.comp == 9
                        return sqrt(2.0^2 + 0.5^2) + 1.5 + 0.5 + sqrt(x.p[1]^2 + (x.p[2]-0.5)^2)
                    else 
                        return sqrt(2.0^2 + 0.5^2) + 1.5 + 1.0 + norm(x.p)
                    end

                elseif x.comp == 2
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], 0.0, 2.0, T-t)
                    if target_is_reached
                        return 0.0
                    else
                        return sqrt(endpoint[1]^2 + (endpoint[2]-2.0)^2)
                    end

                elseif x.comp == 6
                    if 0.5-x.p[1] > (T-t)*wind # we don't reach the panel 2 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + 0.5-x.p[1] - (T-t)*wind
                    else
                        return V2(t + (0.5-x.p[1])/wind, GPoint(2,[-0.5,0.0]))
                    end

                elseif x.comp == 1
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], 0.5, 0.0, T-t)
                    if target_is_reached
                        return return V2(T-surplus_of_time, GPoint(6,[0.0]))
                    else
                        return sqrt(2.0^2 + 0.5^2) + 0.5 + sqrt((endpoint[1]-0.5)^2 + endpoint[2]^2)
                    end

                elseif x.comp == 5
                    if x.p[1]-(-1.5) > T-t # we don't reach the panel 2 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + (x.p[1]+1.5) - (T-t)
                    else
                        return V2(t + x.p[1] + 1.5, GPoint(2,[0.5,0.0]))
                    end

                elseif x.comp == 3
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], -0.5, 0.0, T-t)
                    if target_is_reached
                        return return V2(T-surplus_of_time, GPoint(5,[1.5]))
                    else
                        return sqrt(2.0^2 + 0.5^2) + 3.0 + sqrt((endpoint[1]+0.5)^2 + endpoint[2]^2)
                    end

                elseif x.comp == 7
                    if x.p[1] > T-t # we don't reach the panel 3 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + 3.0 + 1.0 + x.p[1] - (T-t)
                    else
                        return V2(t + x.p[1], GPoint(3,[0.5,0.0]))
                    end

                elseif x.comp == 4
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], -0.5, 0.0, T-t)
                    if target_is_reached
                        return return V2(T-surplus_of_time, GPoint(7,[0.5]))
                    else
                        return sqrt(2.0^2 + 0.5^2) + 3.0 + 1.0 + 0.5 + sqrt((endpoint[1]+0.5)^2 + endpoint[2]^2)
                    end

                elseif x.comp == 8
                    endpoint, junction_is_reached, surplus_of_time = move_in_unity(x.p[1], T-t)
                    if junction_is_reached
                        return V2(T-surplus_of_time, GPoint(5,[0.0]))
                    else 
                        return sqrt(2.0^2 + 0.5^2) + 1.5 + abs(endpoint)
                    end

                elseif x.comp == 9 # Harmony: we are only allowed to move at speed one towards Unity (8)
                    if sqrt(x.p[1]^2 + (x.p[2]-0.5)^2) > T-t # we don't reach 8 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + 1.5 + 0.5 + sqrt(x.p[1]^2 + (x.p[2]-0.5)^2) - (T-t)
                    else
                        return V2(t + sqrt(x.p[1]^2 + (x.p[2]-0.5)^2), GPoint(8,[-0.5]))
                    end

                else # x.comp = 10, Zarya: we are only allowed to move at speed one towards Unity (8)
                    if norm(x.p) > T-t # we don't reach 8 before the end of the time 
                        return sqrt(2.0^2 + 0.5^2) + 1.5 + 1.0 + norm(x.p) - (T-t)
                    else
                        return V2(t + norm(x.p), GPoint(8,[1.0]))
                    end
                end # distinction x.comp
            end # function V2

            function V3(t::Float64, x::GPoint) # target C3=>[0.0,-1.0]
                if t >= T 
                    if x.comp == 3
                        return sqrt(x.p[1]^2 + (x.p[2]-(-1.0))^2)
                    elseif x.comp == 7
                        return sqrt(1.0^2 + 0.5^2) + x.p[1]
                    elseif x.comp == 4
                        return sqrt(1.0^2 + 0.5^2) + 0.5 + sqrt((x.p[1]-(-0.5))^2 + x.p[2]^2)
                    elseif x.comp == 5
                        return sqrt(1.0^2 + 0.5^2) + 1.5-x.p[1]
                    elseif x.comp == 2
                        return sqrt(1.0^2 + 0.5^2) + 3.0 + sqrt((x.p[1]-0.5)^2 + x.p[2]^2)
                    elseif x.comp == 6
                        return sqrt(1.0^2 + 0.5^2) + 3.0 + 1.0 + 0.5-x.p[1]
                    elseif x.comp == 1
                        return sqrt(1.0^2 + 0.5^2) + 3.0 + 1.0 + 0.5 + sqrt((x.p[1]-0.5)^2 + x.p[2]^2)
                    elseif x.comp == 8
                        return sqrt(1.0^2 + 0.5^2) + 1.5 + abs(x.p[1])
                    elseif x.comp == 9
                        return sqrt(1.0^2 + 0.5^2) + 1.5 + 0.5 + sqrt(x.p[1]^2 + (x.p[2]-0.5)^2)
                    else 
                        return sqrt(1.0^2 + 0.5^2) + 1.5 + 1.0 + norm(x.p)
                    end

                elseif x.comp == 3
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], 0.0, -1.0, T-t)
                    if target_is_reached
                        return 0.0
                    else
                        return sqrt(endpoint[1]^2 + (endpoint[2]-(-1.0))^2)
                    end

                elseif x.comp == 7
                    if T - t <= x.p[1]
                        return sqrt(1.0^2 + 0.5^2) + x.p[1] - (T-t)
                    else 
                        return V3(t+x.p[1], GPoint(3,[0.5,0.0]))
                    end

                elseif x.comp == 4
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], -0.5, 0.0, T-t)
                    if target_is_reached
                        return V3(T-surplus_of_time, GPoint(7,[0.5]))
                    else
                        return sqrt(1.0^2 + 0.5^2) + 0.5 + sqrt((endpoint[1]-(-0.5))^2 + endpoint[2]^2)
                    end

                elseif x.comp == 5
                    if 1.5 - x.p[1] > wind*(T-t) # we don't reach the panel 3 before the end of the time 
                        return sqrt(1.0^2 + 0.5^2) + (1.5-x.p[1]) - (T-t)*wind
                    else
                        return V3(t + (1.5-x.p[1])/wind, GPoint(3,[-0.5,0.0]))
                    end

                elseif x.comp == 2
                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], 0.5, 0.0, T-t)
                    if target_is_reached
                        return V3(T-surplus_of_time, GPoint(5,[-1.5]))
                    else
                        return sqrt(1.0^2 + 0.5^2) + 3.0 + sqrt((endpoint[1]-0.5)^2 + endpoint[2]^2)
                    end

                elseif x.comp == 6
                    if 0.5-x.p[1] > (T-t)*wind # we don't reach the panel 2 before the end of the time 
                        return sqrt(1.0^2 + 0.5^2) + 3.0 + 1.0 + 0.5-x.p[1] - (T-t)*wind
                    else
                        return V3(t + (0.5-x.p[1])/wind, GPoint(2,[-0.5,0.0]))
                    end

                elseif x.comp == 1

                    endpoint, target_is_reached, surplus_of_time = move_in_panel(x.p[1], x.p[2], 0.5, 0.0, T-t)
                    if target_is_reached
                        return V3(T-surplus_of_time, GPoint(6,[0.0]))
                    else
                        return sqrt(1.0^2 + 0.5^2) + 3.0 + 1.0 + 0.5 + sqrt((endpoint[1]-0.5)^2 + endpoint[2]^2)
                    end

                elseif x.comp == 8
                    endpoint, junction_is_reached, surplus_of_time = move_in_unity(x.p[1], T-t)
                    if junction_is_reached
                        return V3(T-surplus_of_time, GPoint(5,[0.0]))
                    else 
                        return sqrt(1.0^2 + 0.5^2) + 1.5 + abs(endpoint)
                    end

                elseif x.comp == 9 # Harmony: we are only allowed to move at speed one towards Unity (8)
                    if sqrt(x.p[1]^2 + (x.p[2]-0.5)^2) > T-t # we don't reach 8 before the end of the time 
                        return sqrt(1.0^2 + 0.5^2) + 1.5 + 0.5 + sqrt(x.p[1]^2 + (x.p[2]-0.5)^2) - (T-t)
                    else
                        return V3(t + sqrt(x.p[1]^2 + (x.p[2]-0.5)^2), GPoint(8,[-0.5]))
                    end

                else # x.comp = 10, Zarya: we are only allowed to move at speed one towards Unity (8)
                    if norm(x.p) > T-t # we don't reach 8 before the end of the time 
                        return sqrt(1.0^2 + 0.5^2) + 1.5 + 1.0 + norm(x.p) - (T-t)
                    else
                        return V3(t + norm(x.p), GPoint(8,[1.0]))
                    end
                end # distinction x.comp
            end # function V3
            
            res = MathFunction{Scalar}(id, (t::Float64, x::GPoint) -> min(V1(t,x), min(V2(t,x),V3(t,x))))

        end

    else
        if probe 
            return false 
        else
            throw(ErrorException("unknown analytical solution for domain $(domain.id), dynamic $dynamicid and cost $costid"))
        end
    end

    verb(verbose, "... done.")
    return res 
end

export get_analytical_yu
function get_analytical_yu(domain::Gluing, dynamic::MathFunction{Dynamic}, costid::String, T::Float64, N::Int64, xx::GPoint, verbose::Bool=false)
    verb(verbose, "Computing analytical (y,u)... ")

    if (dynamic.id == "eikonal" && costid == "norm")

        y = Vector{GPoint}(undef, N+1)
        y[end] = GPoint(5,[0.0])
        vv = Velocity(xx, 1.0, 1.0)
        for n in N:-1:1
            y[n] = get_exponential(domain, vv, y[n+1], T/N)
        end

        u = Vector{Velocity}(undef, N)
        for n in 1:N
            dyns = dynamic.call(y[n])
            iu = argmin([get_distance(domain, y[n+1], get_exponential(domain, dyn, y[n], T/N)) for dyn in dyns])
            u[n] = dyns[iu]
        end

    else
        throw(ErrorException("unknown analytical solution for domain $(domain.id), dynamic $(dynamic.id) and cost $costid"))
    end

    verb(verbose, "... done.")
    return y, u
end
