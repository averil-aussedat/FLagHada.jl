"""
Small function to search shortest path in a network
embedded in an euclidian space.
"""
function get_distances(juncpoints::Vector{JuncP}, plotter::Network_plotter, verbose::Bool=false)
    verb(verbose, "Computing distances...")

    njunctions = length(plotter.junccoords)
    nit_max = 1000

    distances = Inf .* ones(njunctions, njunctions)

    for (ip, p) in zip(1:njunctions, juncpoints)
        # the algorithm starts from p
        distances[ip,ip] = 0.0
        toverify = [(ip,jp) for jp in p.neighbours]

        verb(verbose, "ip $ip, init toverify: ", toverify)

        it = 0
        while ((length(toverify) > 0) && (it <= nit_max)) # ceinture, bretelles
            toverify_next = []

            verb(verbose, "ip $ip, \tit $it, toverify beginning : ", toverify)

            for (istart, istop) in toverify
                thedist = norm(plotter.junccoords[istop] - plotter.junccoords[istart])
                if (distances[ip,istop] > distances[ip,istart] + thedist)
                    for neigh in juncpoints[istop].neighbours
                        push!(toverify_next, (istop, neigh))
                    end
                    distances[ip,istop] = distances[ip,istart] + thedist
                end
            end

            toverify = copy(toverify_next)
            verb(verbose, "ip $ip, \tit $it, toverify end : ", toverify)
            it += 1
        end
        if (it >= nit_max)
            throw(ErrorException("$(nit_max) reached iterations for ip = $ip"))
        end

        verb(verbose, "ip $ip, terminal line: ", distances[ip,:], "\n")
    end

    verb(verbose, "... done.")
    return distances
end

export get_network_and_plotter
"""
    Return a network domain and a "plotter" of 2-dimensional coordinates for graphs.
"""
function get_network_and_plotter(id::String, verbose::Bool=false)
    verb(verbose, "Creating $id...")

    if (id == "segment")

        diameter = 1.0
        distances = [
            0.0         diameter; 
            diameter    0.0
        ]
        domain = Network(id, 2, [JuncP(1,[2]), JuncP(2,[1])], diameter, Inf, distances)
        plotter = Network_plotter([diameter/2 * [-1.0,0.0], diameter/2 * [1.0, 0.0]])

    elseif (id == "linear")

        junccoords = [-1.0 -0.5 0.0 1.0 3.0]
        npoints = length(junccoords)
        diameter = junccoords[end] - junccoords[1]
        distances = zeros(npoints, npoints)
        for (i1,j1) in enumerate(junccoords)
            for (i2,j2) in enumerate(junccoords)
                distances[i1,i2] = abs(j2 - j1)
            end
        end
        domain = Network(id, npoints, [JuncP(i,[j for j in [i-1,i+1] if (j>0 && j<=npoints)]) for i in 1:npoints], diameter, Inf, distances)
        plotter = Network_plotter(vec([[x,x] for (ii,x) in enumerate(junccoords)]))

    elseif (id == "tripod")

        #      2
        #      |
        #      |
        #      1
        #     / \
        #    /   \
        #   3     4

        diameter = 1.0
        distances = [
            0.0         diameter    diameter    diameter; 
            diameter    0.0         2*diameter  2*diameter;
            diameter    2*diameter  0.0         2*diameter;
            diameter    2*diameter  2*diameter  0.0
        ]
        domain = Network(id, 4, [JuncP(1,[2,3,4]), JuncP(2,[1]), JuncP(3,[1]), JuncP(4,[1])], 2*diameter, Inf, distances)
        plotter = Network_plotter(vcat([[0.0,0.0]], [diameter .* [cos(th),sin(th)] for th in [π/2, π/2+2π/3, π/2 + 4π/3]]))

    elseif (id == "CAT0")

        junccoords = [
            [ 0.0, 0.0],[-1.0, 0.2],[-1.8,-0.1],[-2.0, 1.2],[-0.5, 1.5],
            [ 0.1,-1.0],[ 2.3,-1.5],[ 1.0,-1.8],[-0.5,-1.9],[-1.0,-1.6],
            [-2.0,-2.0],[-2.0,-1.1],[-1.2,-0.8],[-0.5,-0.9],[ 2.0, 1.5],
            [ 2.5, 1.8],[ 2.2, 0.3],[ 2.4,-1.0],[ 1.0,-0.2],[ 0.8, 0.5],
            [ 1.4, 1.0],[ 2.8, 0.0],[-0.9,-0.5]]
        plotter = Network_plotter(junccoords)
        juncpoints = [
            JuncP(1, [2,6,20]),
            JuncP(2, [1,3,4,5]),
            JuncP(3, [2]),
            JuncP(4, [2]),
            JuncP(5, [2]),
            JuncP(6, [1,7,8,9,14]),
            JuncP(7, [6]),
            JuncP(8, [6]),
            JuncP(9, [6,10]),
            JuncP(10,[9,11,12,13]),
            JuncP(11,[10]),
            JuncP(12,[10]),
            JuncP(13,[10]),
            JuncP(14,[6,23]),
            JuncP(15,[21,16,17]),
            JuncP(16,[15,22]),
            JuncP(17,[15,18,19]),
            JuncP(18,[17]),
            JuncP(19,[17]),
            JuncP(20,[1,21]),
            JuncP(21,[20,15]),
            JuncP(22,[16]),
            JuncP(23,[14])
        ]
        distances = get_distances(juncpoints, plotter, false)
        domain = Network(id, length(juncpoints), juncpoints, 2.0*maximum(distances), Inf, distances)

    elseif (id == "triangle")

        sidesize = 1.0
        distances = [
            0.0         sidesize    sidesize;
            sidesize    0.0         sidesize;
            sidesize    sidesize    0.0
        ]
        domain = Network(id, 3, [JuncP(1,[2,3]), JuncP(2,[1,3]), JuncP(3,[1,2])], 1.5*sidesize, 1.5*sidesize, distances)
        plotter = Network_plotter([[sidesize/2, 0.0], [0.0, sqrt(3)/2*sidesize], [-sidesize/2, 0.0]])

    elseif (id == "circle")

        npoints = 7
        distances = zeros(npoints,npoints)
        for ii in 1:npoints
            for jj in 1:npoints
                distances[ii,jj] = min(abs(ii-jj), npoints-jj+ii, npoints-ii+jj)
            end
        end
        domain = Network(id, npoints, [JuncP(ii,[((ii-2+npoints) % npoints)+1,(ii % npoints)+1]) for ii in 1:npoints], 0.5*npoints, 0.5*npoints, distances)
        plotter = Network_plotter(vec([[cos(ii * 2π/npoints),sin(ii * 2π/npoints)] for ii in 0:(npoints-1)]))

    elseif (id == "intricate")
        juncpoints = [
            JuncP(1,[2,7,8,11]),    # Paris 
            JuncP(2,[1,3,6]),       # Lyon 
            JuncP(3,[2,4,5]),       # Avignon
            JuncP(4,[3,16]),        # Marseille
            JuncP(5,[3,15]),        # Barcelone
            JuncP(6,[2]),           # Turin
            JuncP(7,[1]),           # Strasbourg
            JuncP(8,[1,9,10]),      # Lille
            JuncP(9,[8]),           # Bruxelles
            JuncP(10,[8]),          # Londres
            JuncP(11,[1,12,13]),    # Le Mans
            JuncP(12,[11]),         # Rennes
            JuncP(13,[11,14]),      # Bordeaux
            JuncP(14,[13]),         # Dax 
            JuncP(15,[5]),          # Madrid (Puerta de Atocha)
            JuncP(16,[4]),          # Nice
        ]
        junccoords = [
            [ 0.0,  0.0],  # Paris
            [ 1.54,-3.05], # Lyon 
            [ 1.05,-5.1],  # Avignon
            [ 1.83,-5.85], # Marseille
            [-0.65,-7.8],  # Barcelone
            [ 4.5, -3.5],  # Turin
            [ 3.65,-0.3],  # Strasbourg
            [ 0.53, 2.1],  # Lille
            [ 1.8,  2.3],  # Bruxelles
            [-1.7,  4.0],  # Londres
            [-1.9, -0.8],  # Le Mans
            [-3.05,-0.25], # Rennes
            [-2.55,-3.55], # Bordeaux
            [-3.15,-4.95], # Dax 
            [-2.35,-7.9],  # Madrid (Puerta de Atocha)
            [ 3.4,-5.65],  # Nice
        ]
        plotter = Network_plotter(junccoords)
        distances = get_distances(juncpoints, plotter)
        diameter = maximum(distances)
        domain = Network(id, length(juncpoints), juncpoints, diameter, 2*diameter, distances)

    else
        throw(ErrorException("unknown domain id '$id'"))
    end

    verb(verbose, "... done.")
    return domain, plotter
end

export get_Scalar
function get_Scalar(id::String, domain::Network, verbose::Bool=false)
    verb(verbose, "Creating scalar function $id...")
    
    if (id == "one")  

        res = MathFunction{Scalar}(id, x::NetPoint -> 1.0)
    
    elseif (id == "norm")

        res = MathFunction{Scalar}(id, x::NetPoint -> get_distance(domain, domain.juncpoints[1], x))

    elseif (id == "norm_middle" && domain.id == "segment")

        res = MathFunction{Scalar}(id, x::NetPoint -> get_distance(domain, EdgeP(1,2,0.5,0.5), x))

    elseif (id == "norm_truncated")

        res = MathFunction{Scalar}(id, x::NetPoint -> min(get_distance(domain, domain.juncpoints[1], x),1.0))

    elseif (id == "sinus")

        res = MathFunction{Scalar}(id, x::NetPoint -> sin(get_distance(domain, domain.juncpoints[1], x)))

    elseif ((id == "coffee") && (domain.id == "coffee"))

        the_machine_id = 14
        res = MathFunction{Scalar}(id, x::NetPoint -> get_distance(domain, domain.juncpoints[the_machine_id], x))

    elseif (id == "ferian")

        if (domain.id == "trijunction")
            the_machine_id = 2
        elseif (domain.id == "tworoads") 
            the_machine_id = 6
        else 
            the_machine_id = 1
        end
        
        res = MathFunction{Scalar}(id, x::NetPoint -> get_distance(domain, domain.juncpoints[the_machine_id], x))

    elseif (id == "dist_south" && domain.id == "intricate")

        south = [domain.juncpoints[4], domain.juncpoints[5], domain.juncpoints[13]] # Marseille, Barcelone, Bordeaux
        res = MathFunction{Scalar}(id, x::NetPoint -> minimum([get_distance(domain, s, x) for s in south]))

    elseif (id == "dist_Bordeaux" && domain.id == "intricate")

        res = MathFunction{Scalar}(id, x::NetPoint -> get_distance(domain, domain.juncpoints[13], x))

    elseif (id == "dist_Marseille" && domain.id == "intricate")

        res = MathFunction{Scalar}(id, x::NetPoint -> get_distance(domain, domain.juncpoints[4], x))

    elseif (id == "dist_Barcelone" && domain.id == "intricate")

        res = MathFunction{Scalar}(id, x::NetPoint -> get_distance(domain, domain.juncpoints[5], x))

    else
        throw(ErrorException("unknown scalar id '$id'"))
    end

    verb(verbose, "... done.")
    return res
end

export get_Dynamic
function get_Dynamic(id::String, domain::Network, verbose::Bool=false)
    verb(verbose, "Creating dynamic $id...")

    bdries = [junc for junc in domain.juncpoints if length(junc.neighbours) <= 1]

    if (id == "eikonal")

        res = MathFunction{Dynamic}(id, x::NetPoint -> [Velocity(x, 0.0, 1.0)] ∪ [Velocity(bdry, 1.0, 1.0) for bdry in bdries])
        # res = MathFunction{Dynamic}(id, x::NetPoint -> [Velocity(x, 0.0, 1.0)] ∪ [Velocity(bdry, m*0.2, 1.0) for bdry in bdries for m in 1:5])

    elseif (domain.id == "tripod" && id == "tripod_twoaims")

        res = MathFunction{Dynamic}(id, x::NetPoint -> [Velocity(domain.juncpoints[3], 1.0, 1.0), Velocity(domain.juncpoints[4], 1.0, 1.0)])

    elseif (domain.id == "intricate" && id == "sncf")

        # to Londres (11) with constant speed 
        # to Nice (16) with slower and slower speed as we get closer to [M,N]
        # to Puerta de Atocha (Madrid) (15) with constant speed 
        # to Dax (14) with low speed but gets faster as we get closer to [D,Bₒ]

        res = MathFunction{Dynamic}(id, 
            x::NetPoint -> 
                if (typeof(x) == JuncP && x.id ∈ [13,14]) || (typeof(x) == EdgeP && (14 == x.left || 14 == x.right))
                    # ax [Dax,Bordeaux]
                    return union(
                        [Velocity(domain.juncpoints[10], 1.0, 1.0)], # Londres
                        [Velocity(domain.juncpoints[16], 1.2, 1.0)], # Nice
                        [Velocity(domain.juncpoints[15], 1.0, 1.0)], # Puerta de Atocha (Madrid)
                        [Velocity(domain.juncpoints[14], 1.5, 1.0)], # Dax 
                    )
                elseif (typeof(x) == JuncP && x.id ∈ [4,16]) || (typeof(x) == EdgeP && (16 == x.left || 16 == x.right))
                    # ax [Nice, Marseille]
                    return union(
                        [Velocity(domain.juncpoints[10], 1.0, 1.0)], # Londres
                        [Velocity(domain.juncpoints[16], 0.2, 1.0)], # Nice
                        [Velocity(domain.juncpoints[15], 1.0, 1.0)], # Puerta de Atocha (Madrid)
                        [Velocity(domain.juncpoints[14], 0.9, 1.0)], # Dax 
                    )
                else
                    return union(
                        [Velocity(domain.juncpoints[10], 1.0, 1.0)], # Londres
                        [Velocity(domain.juncpoints[16], 0.2 + min(1.0, get_distance(domain, domain.juncpoints[4], x)/domain.distances[2,4]), 1.0)], # Nice 
                        [Velocity(domain.juncpoints[15], 1.0, 1.0)], # Puerta de Atocha (Madrid)
                        [Velocity(domain.juncpoints[14], 1.5 - 0.6 * min(1.0, get_distance(domain, domain.juncpoints[13], x)/domain.distances[13,11]), 1.0)], # Dax
                    )
                end
        )
            
    else
        throw(ErrorException("unknown dynamic id '$id'"))
    end

    verb(verbose, "... done.")
    return res
end

export get_Analytical
"""
    Return the exact solution as a MathFunction{Scalar}.
    When probe = true, only return the existence of a solution.
"""
function get_Analytical(domain::Network, dynamicid::String, costid::String, T::Float64, probe::Bool=false, verbose::Bool=false)
    verb(verbose, "Computing analytical V... ")

    id = domain.id * "_" * dynamicid * "_" * costid

    if (dynamicid == "eikonal" && costid == "norm") || (domain.id == "tripod" && dynamicid == "tripod_twoaims" && costid == "norm")

        if probe 
            return true 
        else
            res = MathFunction{Scalar}(id, (t::Float64, x::NetPoint) -> max(0.0, get_distance(domain, domain.juncpoints[1], x) - (T-t)))
        end

    elseif (domain.id == "intricate" && dynamicid == "sncf" && costid in ["dist_south","dist_Bordeaux","dist_Marseille","dist_Barcelone"])

        if probe 
            return true 
        else
            # the true value function is the minimum between the value functions associated to the three 
            # objective points : Bordeaux (13), Marseille (4) and Barcelone (5)

            function go_constant_speed(FUNC::Function, t::Float64, x::NetPoint, dir::NetPoint, speed::Float64, aim::NetPoint)
                if get_distance(domain, x, dir) >= speed * (T-t) # we do not pass after dir 
                    return max(0.0, get_distance(domain, x, aim) - speed*(T-t))
                else # reaching dir before the end of the allowed time 
                    return FUNC(t + get_distance(domain, x, dir)/speed, dir) # jumping there and restarting
                end 
            end

            # First : Bordeaux 
            function V13(t::Float64, x::NetPoint)
                if (t >= T)
                    return get_distance(domain, x, domain.juncpoints[13])
                elseif ((x in [domain.juncpoints[13],domain.juncpoints[14]]) || (typeof(x)==EdgeP && 14 in [x.left, x.right]))
                    # Dax (14) - Bordeaux (13) : take the speed 1.2 towards Nice until Bordeaux
                    return max(0.0, get_distance(domain, x, domain.juncpoints[13]) - 1.2*(T-t)) 
                elseif ((x == domain.juncpoints[11]) || (typeof(x)==EdgeP && 13 in [x.left, x.right]))
                    # we are between Le Mans (11) and Bordeaux (13), speedup
                    return max(0.0, exp(0.6*(T-t)/domain.distances[13,11]) * get_distance(domain,x,domain.juncpoints[13]) + (1.0 - exp(0.6*(T-t)/domain.distances[13,11])) * domain.distances[13,11] * 1.5/0.6)
                elseif (x == domain.juncpoints[12]) || (typeof(x) == EdgeP && 12 in [x.left, x.right])
                    # Rennes-Le Mans : we take the speed 1.2 towards Nice until Le Mans (11)
                    return go_constant_speed(V13, t, x, domain.juncpoints[11], 1.2, domain.juncpoints[13])
                elseif (x == domain.juncpoints[1]) || (typeof(x) == EdgeP && 1 in [x.left, x.right] && 11 in [x.left, x.right])
                    # Paris-Le Mans : we take the speed 0.9 towards Dax 
                    return go_constant_speed(V13, t, x, domain.juncpoints[11], 0.9, domain.juncpoints[13])
                elseif (typeof(x) == JuncP && x.id in [7,8,9,10]) || (typeof(x) == EdgeP && (x.left in [7,8,9,10] || x.right in [7,8,9,10]))
                    # Near Lille or Strasbourg: we take the speed 1.2 towards Nice until Paris
                    return go_constant_speed(V13, t, x, domain.juncpoints[1], 1.2, domain.juncpoints[13])
                elseif (typeof(x) == JuncP && x.id == 6) || (typeof(x) == EdgeP && (6 in [x.left,x.right]))
                    # Turin-Lyon : speed 1.2 towards Nice until Lyon
                    return go_constant_speed(V13, t, x, domain.juncpoints[2], 1.2, domain.juncpoints[13])
                elseif (typeof(x) == JuncP && x.id in [2,3,4]) || (typeof(x) == EdgeP && (2 in [x.left,x.right] || 4 in [x.left, x.right]))
                    # Axes Lyon-Paris, Lyon-Avignon and Avignon-Marseille: speed 1 towards Londres until Paris  
                    return go_constant_speed(V13, t, x, domain.juncpoints[1], 1.0, domain.juncpoints[13])
                elseif (typeof(x) == JuncP && x.id == 16) || (typeof(x) == EdgeP && 16 in [x.left,x.right])
                    # Ax Nice-Marseille, take speed 1 towards Puerta de Atocha until Marseille 
                    return go_constant_speed(V13, t, x, domain.juncpoints[4], 1.0, domain.juncpoints[13])
                elseif (typeof(x) == JuncP && x.id == 15) || (typeof(x) == EdgeP && 15 in [x.left,x.right])
                    # Ax Puerta-Barcelone, take speed 1.2 towards Nice until Barcelona 
                    return go_constant_speed(V13, t, x, domain.juncpoints[5], 1.2, domain.juncpoints[13])
                elseif (typeof(x) == JuncP && x.id == 5) || (typeof(x) == EdgeP && 5 in [x.left,x.right])
                    # From Barcelone, take the optimal way to go to Avignon : first follow Nice, then Londres
                    # threshold : 0.2 + min(1.0, d(x,mar)/d(lyon,mar)) = 1.0 quand d(x,mar) = 0.8 d(lyon,mar)
                    if get_distance(domain, x, domain.juncpoints[4]) <= 0.8 * domain.distances[2,4]
                        # follow Londres with speed 1 until Avignon 
                        return go_constant_speed(V13, t, x, domain.juncpoints[3], 1.0, domain.juncpoints[13])
                    else 
                        # follow Nice 
                        if get_distance(domain, x, domain.juncpoints[3]) <= domain.distances[2,3]
                            # slowdown zone. reachtime is the time at which we switch to the dynamic towards Londres
                            # reachtime = t - domain.distances[2,4] * log(domain.distances[2,3] / (get_distance(domain, x, domain.juncpoints[3]) + 0.2 * domain.distances[2,4]))
                            reachtime = t + domain.distances[2,4] * log(0.2 + get_distance(domain, x, domain.juncpoints[4])/domain.distances[2,4])
                            if reachtime < T # we reach the critical point before T 
                                if 0.8*domain.distances[2,4]-domain.distances[3,4] >= T-reachtime # end between the critical point and Avignon
                                    return domain.distances[3,13] + 0.8*domain.distances[2,4]-domain.distances[3,4] - (T-reachtime)
                                else # reaching Avignon before the end of the allowed time 
                                    return V13(reachtime + 0.8*domain.distances[2,4]-domain.distances[3,4], domain.juncpoints[3]) # jumping there and restarting
                                end 
                            else 
                                return domain.distances[3,13] + exp(-(T-t)/domain.distances[2,4]) * get_distance(domain, x, domain.juncpoints[4]) -
                                        0.2 * (1.0 - exp(-(T-t)/domain.distances[2,4])) * domain.distances[2,4] - domain.distances[3,4]                                        
                            end
                        else
                            # constant speed zone : 1.2
                            return go_constant_speed(V13, t, x, EdgeP(3,5,domain.distances[2,3],domain.distances[5,3]-domain.distances[2,3]), 1.2, domain.juncpoints[13])
                        end
                    end
                else 
                    throw(ErrorException("case of point $x not covered in V13"))
                end
            end

            # Second : Marseille 
            function V4(t::Float64, x::NetPoint)
                if (t >= T)
                    return get_distance(domain, x, domain.juncpoints[4])
                elseif ((typeof(x) == JuncP && x.id in [14]) || (typeof(x) == EdgeP && 14 in [x.left,x.right]))
                    # ax [Dax,Bordeaux[; going towards Nice at speed 1.2 until Bordeaux 
                    return go_constant_speed(V4, t, x, domain.juncpoints[13], 1.2, domain.juncpoints[4])
                elseif ((typeof(x) == JuncP && x.id in [12,13]) || (typeof(x) == EdgeP && (12 in [x.left,x.right] || 13 in [x.left,x.right])))
                    # axes [Rennes,Le Mans[ and [Bordeaux, Le Mans[; going towards Nice until Le Mans, speed 1.2
                    return go_constant_speed(V4, t, x, domain.juncpoints[11], 1.2, domain.juncpoints[4])
                elseif ((typeof(x) == JuncP && x.id in [9,10]) || (typeof(x) == EdgeP && (9 in [x.left,x.right] || 10 in [x.left,x.right])))
                    # axes [Londres,Lille[ and [Bruxelles,Lille[; going towards Nice until Lille, speed 1.2
                    return go_constant_speed(V4, t, x, domain.juncpoints[8], 1.2, domain.juncpoints[4])
                elseif ((typeof(x) == JuncP && x.id in [7,8,11]) || (typeof(x) == EdgeP && (7 in [x.left,x.right] || 8 in [x.left,x.right] || 11 in [x.left,x.right])))
                    # axes [Le Mans,Paris[, [Lille,Paris[ and [Strasbourg,Paris[; going towards Nice until Paris, speed 1.2
                    return go_constant_speed(V4, t, x, domain.juncpoints[1], 1.2, domain.juncpoints[4])
                elseif ((typeof(x) == JuncP && x.id in [1,6]) || (typeof(x) == EdgeP && (1 in [x.left,x.right] || 6 in [x.left,x.right])))
                    # axes [Paris,Lyon[ and [Turin,Lyon[; going towards Nice until Lyon, speed 1.2
                    return go_constant_speed(V4, t, x, domain.juncpoints[2], 1.2, domain.juncpoints[4])
                elseif (typeof(x) == JuncP && x.id == 16) || (typeof(x) == EdgeP && 16 in [x.left,x.right])
                    # Ax Nice-Marseille, take speed 1 towards Puerta de Atocha until Marseille 
                    return max(0.0, get_distance(domain, x, domain.juncpoints[4]) - (T-t))
                elseif (typeof(x) == JuncP && x.id == 15) || (typeof(x) == EdgeP && 15 in [x.left,x.right])
                    # Ax Puerta-Barcelone, take speed 1.2 towards Nice until Barcelona 
                    return go_constant_speed(V4, t, x, domain.juncpoints[5], 1.2, domain.juncpoints[4])
                elseif ((typeof(x) == JuncP && x.id in [3,4]) || (typeof(x) == EdgeP && 4 in [x.left,x.right]))
                    # ax [Avignon, Marseille]; going towards Nice, slowdown 
                    return max(0.0, exp(-(T-t)/domain.distances[2,4]) * get_distance(domain, x, domain.juncpoints[4]) - 0.2 * (1.0 - exp(-(T-t)/domain.distances[2,4])) * domain.distances[2,4])
                else
                    # Symmetric situation for the axes [Barcelone, Avignon[ and [Lyon, Avignon[
                    # First a constant speed zone towards Nice (empty in the case of Lyon)
                    # then a slowdown zone 
                    # then a constant speed zone when reaching velocity 1.0, resp. towards Londres and Puerta de Atocha 
                    # threshold : 0.2 + min(1.0, d(x,mar)/d(lyon,mar)) = 1.0 quand d(x,mar) = 0.8 d(lyon,mar)
                    if get_distance(domain, x, domain.juncpoints[4]) <= 0.8 * domain.distances[2,4]
                        # follow Londres/Puerta with speed 1 until Avignon 
                        return go_constant_speed(V4, t, x, domain.juncpoints[3], 1.0, domain.juncpoints[4])
                    else 
                        # follow Nice 
                        if get_distance(domain, x, domain.juncpoints[3]) <= domain.distances[2,3]
                            # slowdown zone. reachtime is the time at which we switch to the 1.0 dynamic
                            reachtime = t + domain.distances[2,4] * log(0.2 + get_distance(domain, x, domain.juncpoints[4])/domain.distances[2,4])
                            if reachtime < T # we reach the critical point before T 
                                if 0.8*domain.distances[2,4]-domain.distances[3,4] >= T-reachtime # end between the critical point and Avignon
                                    return 0.8*domain.distances[2,4] - (T-reachtime)
                                else # reaching Avignon before the end of the allowed time 
                                    return V4(reachtime + 0.8*domain.distances[2,4]-domain.distances[3,4], domain.juncpoints[3]) # jumping there and restarting
                                end 
                            else 
                                return exp(-(T-t)/domain.distances[2,4]) * get_distance(domain, x, domain.juncpoints[4]) -
                                        0.2 * (1.0 - exp(-(T-t)/domain.distances[2,4])) * domain.distances[2,4]
                            end
                        else
                            # constant speed zone : 1.2. Only entered for the ax [Barcelone, Avignon[
                            return go_constant_speed(V4, t, x, EdgeP(3,5,domain.distances[2,3],domain.distances[5,3]-domain.distances[2,3]), 1.2, domain.juncpoints[4])
                        end
                    end
                end
            end

            # Third : Barcelone
            function V5(t::Float64, x::NetPoint)
                if (t >= T)
                    return get_distance(domain, x, domain.juncpoints[5])
                elseif ((typeof(x) == JuncP && x.id in [14]) || (typeof(x) == EdgeP && 14 in [x.left,x.right]))
                    # ax [Dax,Bordeaux[; going towards Nice at speed 1.2 until Bordeaux
                    return go_constant_speed(V5, t, x, domain.juncpoints[13], 1.2, domain.juncpoints[5])
                elseif ((typeof(x) == JuncP && x.id in [12,13]) || (typeof(x) == EdgeP && (12 in [x.left,x.right] || 13 in [x.left,x.right])))
                    # axes [Rennes,Le Mans[ and [Bordeaux, Le Mans[; going towards Nice until Le Mans, speed 1.2
                    return go_constant_speed(V5, t, x, domain.juncpoints[11], 1.2, domain.juncpoints[5])
                elseif ((typeof(x) == JuncP && x.id in [9,10]) || (typeof(x) == EdgeP && (9 in [x.left,x.right] || 10 in [x.left,x.right])))
                    # axes [Londres,Lille[ and [Bruxelles,Lille[; going towards Nice until Lille, speed 1.2
                    return go_constant_speed(V5, t, x, domain.juncpoints[8], 1.2, domain.juncpoints[5])
                elseif ((typeof(x) == JuncP && x.id in [7,8,11]) || (typeof(x) == EdgeP && (7 in [x.left,x.right] || 8 in [x.left,x.right] || 11 in [x.left,x.right])))
                    # axes [Le Mans,Paris[, [Lille,Paris[ and [Strasbourg,Paris[; going towards Nice until Paris, speed 1.2
                    return go_constant_speed(V5, t, x, domain.juncpoints[1], 1.2, domain.juncpoints[5])
                elseif ((typeof(x) == JuncP && x.id in [1,6]) || (typeof(x) == EdgeP && (1 in [x.left,x.right] || 6 in [x.left,x.right])))
                    # axes [Paris,Lyon[ and [Turin,Lyon[; going towards Nice until Lyon, speed 1.2
                    return go_constant_speed(V5, t, x, domain.juncpoints[2], 1.2, domain.juncpoints[5])
                elseif ((x == domain.juncpoints[4]) || (typeof(x) == EdgeP && 4 in [x.left,x.right]))
                    # ax [Marseille, Avignon[; going towards Puerta de Atocha until Avignon, speed 1.0 
                    return go_constant_speed(V5, t, x, domain.juncpoints[3], 1.0, domain.juncpoints[5])
                elseif (typeof(x) == JuncP && x.id == 16) || (typeof(x) == EdgeP && 16 in [x.left,x.right])
                    # Ax Nice-Marseille, take speed 1 towards Puerta de Atocha until Marseille 
                    return go_constant_speed(V5, t, x, domain.juncpoints[4], 1.0, domain.juncpoints[5])
                elseif (typeof(x) == JuncP && x.id == 15) || (typeof(x) == EdgeP && 15 in [x.left,x.right])
                    # Ax Puerta-Barcelone, take speed 1.2 towards Nice until Barcelona 
                    return max(0.0, get_distance(domain, x, domain.juncpoints[5]) - 1.2 * (T-t))
                elseif ((typeof(x) == JuncP && x.id in [3,5]) || (typeof(x) == EdgeP && 5 in [x.left,x.right]))
                    # ax [Avignon, Barcelone]; going towards Puerta de Atocha, speed 1.0 
                    return max(0.0, get_distance(domain, x, domain.juncpoints[5]) - (T-t))
                else
                    # First a slowdown zone, then a constant speed zone when reaching velocity 1.0 towards Puerta de Atocha 
                    # threshold : 0.2 + min(1.0, d(x,mar)/d(lyon,mar)) = 1.0 quand d(x,mar) = 0.8 d(lyon,mar)
                    if get_distance(domain, x, domain.juncpoints[4]) <= 0.8 * domain.distances[2,4]
                        # follow Puerta de Atocha with speed 1 until Avignon 
                        return go_constant_speed(V5, t, x, domain.juncpoints[3], 1.0, domain.juncpoints[5])
                    else 
                        # follow Nice, slowdown zone
                        # reachtime is the time at which we switch to the 1.0 dynamic
                        reachtime = t + domain.distances[2,4] * log(0.2 + get_distance(domain, x, domain.juncpoints[4])/domain.distances[2,4])
                        if reachtime < T # we reach the critical point before T 
                            if 0.8*domain.distances[2,4]-domain.distances[3,4] >= T-reachtime # end between the critical point and Avignon
                                return 0.8*domain.distances[2,4]-domain.distances[3,4]+domain.distances[3,5] - (T-reachtime)
                            else # reaching Avignon before the end of the allowed time 
                                return V5(reachtime + 0.8*domain.distances[2,4]-domain.distances[3,4], domain.juncpoints[3]) # jumping there and restarting
                            end 
                        else 
                            return (exp(-(T-t)/domain.distances[2,4]) * get_distance(domain, x, domain.juncpoints[4]) -
                                    0.2 * (1.0 - exp(-(T-t)/domain.distances[2,4])) * domain.distances[2,4]
                                    - domain.distances[3,4] + domain.distances[3,5])
                        end
                    end
                end
            end

            # res = MathFunction{Scalar}(id, (t::Float64, x::NetPoint) -> V13(t,x))
            # res = MathFunction{Scalar}(id, (t::Float64, x::NetPoint) -> V4(t,x))
            # res = MathFunction{Scalar}(id, (t::Float64, x::NetPoint) -> V5(t,x))
            res = MathFunction{Scalar}(id, (t::Float64, x::NetPoint) -> minimum([V4(t,x),V5(t,x),V13(t,x)]))
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
function get_analytical_yu(domain::Network, dynamic::MathFunction{Dynamic}, costid::String, T::Float64, N::Int64, xx::NetPoint, verbose::Bool=false)
    verb(verbose, "Computing analytical (y,u)... ")

    if (dynamic.id == "eikonal" && costid == "norm")

        y = Vector{NetPoint}(undef, N+1)
        y[end] = domain.juncpoints[1]
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

    elseif (domain.id == "tripod" && dynamic.id == "tripod_twoaims" && costid == "norm")

        y = Vector{NetPoint}(undef, N+1)
        y[1] = xx 
        vv = Velocity(domain.juncpoints[1], 1.0, 1.0)
        for n in 1:N
            y[n+1] = get_exponential(domain, vv, y[n], T/N)
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
