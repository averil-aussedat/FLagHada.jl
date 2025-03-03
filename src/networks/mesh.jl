

function get_mesh(domain::Network, step::Float64, verbose::Bool=false)
    verb(verbose, "Meshing $(domain.id) with step $step...")

    tomesh = trues(domain.njunctions)
    realstep = 0.0
    points = Vector{NetPoint}() # I'm sorry again Pierre for using push!

    juncind = Vector{Int64}(undef, domain.njunctions)
    edgeind = Matrix{Tuple{Int64,Int64}}(undef, domain.njunctions, domain.njunctions)
    [edgeind[i,j] = (0,0) for i in 1:domain.njunctions for j in 1:domain.njunctions]

    for (ipt,pt) in enumerate(domain.juncpoints) 

        push!(points, pt)
        juncind[ipt] = length(points)

        for neigh in pt.neighbours 
            if tomesh[neigh]
                dd = domain.distances[pt.id, neigh]
                nn = floor(Int,dd/step)
                locstep = dd / (nn+1)

                edgeind[ipt,neigh] = (length(points)+1,length(points)+nn)

                for in = 1:nn
                    push!(points, EdgeP(pt.id, neigh, in * locstep, (nn+1-in) * locstep))
                end
                realstep = max(locstep, realstep)
            end
        end

        tomesh[pt.id] = false
    end

    verb(verbose, "... done with $(length(points)) points.")
    return NetMesh(length(points), realstep, points, juncind, edgeind)

end

export search_in_mesh
"""
    Return the index of the closest point of mesh to x
"""
function search_in_mesh(domain::Network, mesh::Mesh, x::NetPoint)
    if typeof(x) == JuncP 
        return mesh.juncind[x.id]
    else
        init, last = mesh.edgeind[min(x.left,x.right),max(x.left,x.right)]
        if (init >= last) # no interior points in the edge 
            return x.left * (x.dist_left <= x.dist_right) + x.right * (x.dist_left > x.dist_right)
        else 
            if 0.5*(x.dist_left + x.dist_right)/(last-init+2) >= x.dist_left * (x.left < x.right) + x.dist_right * (x.right < x.left)
                return mesh.juncind[min(x.left,x.right)]
            elseif x.dist_left + x.dist_right - 0.5*(x.dist_left + x.dist_right)/(last-init+2) <= x.dist_left * (x.left < x.right) + x.dist_right * (x.right < x.left)
                return mesh.juncind[max(x.left,x.right)]
            else 
                return init - 1 + ceil(Int64, (last-init+2) * 
                    (x.dist_left * (x.left < x.right) + x.dist_right * (x.right < x.left) - 0.5*(x.dist_left + x.dist_right)/(last-init+2)) 
                    / (x.dist_left + x.dist_right))
            end
        end
    end 
end

"""
    Return the list of indexes of the mesh 
    between p1 and p2 under the assumption that they share the same edge 
"""
function get_edge_conv(domain::Network, mesh::Mesh, p1::JuncP, p2::JuncP)
    return union(mesh.juncind[p1.id], mesh.juncind[p2.id], map(a -> a[1]:a[2], [mesh.edgeind[min(p1.id,p2.id), max(p1.id, p2.id)]])[1])
end

function get_edge_conv(domain::Network, mesh::Mesh, p1::JuncP, p2::EdgeP)
    if p2.left < p2.right 
        if p1.id == p2.left 
            #  p1 -------- p2 ------- p2.right       ------> meshing direction
            return union(mesh.juncind[p1.id], mesh.edgeind[p2.left,p2.right][1]:search_in_mesh(domain, mesh, p2))
        else 
            # p2.left ------ p2 ----- p2.right = p1  ------> meshing direction
            return union(mesh.juncind[p1.id], search_in_mesh(domain, mesh, p2):mesh.edgeind[p2.left,p2.right][2])
        end
    else 
        if p1.id == p2.left 
            #  p1 -------- p2 ------- p2.right       <------ meshing direction
            return union(mesh.juncind[p1.id], search_in_mesh(domain, mesh, p2):mesh.edgeind[p2.right,p2.left][2])
        else 
            # p2.left ------ p2 ----- p2.right = p1  <------ meshing direction
            return union(mesh.juncind[p1.id], mesh.edgeind[p2.right,p2.left][1]:search_in_mesh(domain, mesh, p2))
        end
    end
end

function get_edge_conv(domain::Network, mesh::Mesh, p1::EdgeP, p2::JuncP)
    return get_edge_conv(domain, mesh, p2, p1)
end

function get_edge_conv(domain::Network, mesh::Mesh, p1::EdgeP, p2::EdgeP)
    return map((a,b) -> min(a,b):max(a,b), search_in_mesh(domain, mesh, p1), search_in_mesh(domain, mesh, p2))
end

export get_mesh_convex
"""
    Return the set of points of the mesh that belong to the convex hull of points.
"""
function get_mesh_convex(domain::Network, mesh::Mesh, points::Vector{<:NetPoint}, aretrimmed::Bool=false)
    if length(points) == 1
        return [search_in_mesh(domain, mesh, points[1])]
    else 
        if aretrimmed
            # keep union([]...) to remove duplicates
            dirs = union([get_direction(domain, points[1], q) for q in points[2:end]]...)
            if length(dirs) == 1
                argmaxdist = argmax([get_distance(domain, points[1], q) for q in points[2:end]])
                if get_distance(domain, points[1], points[1+argmaxdist]) <= get_distance(domain, points[1], domain.juncpoints[dirs[1]])
                    # either all points[2:end] are on the edge [point[1], dirs[1]]
                    return get_edge_conv(domain, mesh, points[1], points[1+argmaxdist])
                else 
                    # or we have to go further 
                    return union(get_edge_conv(domain, mesh, points[1], domain.juncpoints[dirs[1]]), 
                            get_mesh_convex(domain, mesh, union([domain.juncpoints[dirs[1]]], points[2:end]), true)
                        )
                end
            else 
                return union([get_mesh_convex(domain, mesh, union([points[1]], [q for q in points[2:end] if get_direction(domain, points[1], q) == d]), true) for d in dirs]...)
            end
        else
            # theconv = get_convex(domain, points) # trimming
            trimmedpoints = get_convex(domain, points) # trimming
            return get_mesh_convex(domain, mesh, union([mesh.points[search_in_mesh(domain, mesh, p)] for p in trimmedpoints]), true)
        end
    end
end