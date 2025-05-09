using LinearAlgebra, Printf
import JuMP, MathOptInterface, Ipopt, MadNLP


@doc raw"""
    Polyhedron{T}

Implementation of convex polyhedron, used here to compute maximum available torque about a desired axis.
"""
struct Polyhedron{T}
    normals::Matrix{T}
    distances::Vector{T}
    vertices::Matrix{T}
    facets::Vector{Matrix{T}}

    Polyhedron(normals::Matrix{Float64}, distances::Vector{Float64}) = begin
        (n,m) = size(normals)
        if !any(size(normals) .== 3)
            error("Expect a 3xn or nx3 matrix of face normal vectors.")
        end
        if (m != 3) normals = normals' end
        n = 3
        m = length(normals[:,1])
        if length(distances) != m
            error("distances must have same length as number of unit vectors in normals.")
        end

        faces = zeros(size(normals))
        [faces[i,:] = normals[i,:]*distances[i] for i in range(1,m)]
        
        relerr = 1e-3
        L = min(distances...)
        eps = L * relerr

        nconstraints = m
        vertices = Matrix{Float64}(undef, 3, 0)
        for i in 1:nconstraints
            for j in i:nconstraints
                for k in j:nconstraints
                    if i == j || j == k || i == k
                        continue
                    end
                    P = normals[[i,j,k],:]
                    q = distances[[i,j,k]]
                    if rank(P) == 3
                        # check if vertex is solution to LMI for polytope surface
                        #   this check excludes all facet intersections which are external to the polytope.
                        x = P\q
                        if all(normals*x .<= (distances .+ eps))
                            # check if vertex already in list:
                            test_x = [x'*vertex/norm(vertex)/norm(x) .- 1 for vertex in eachcol(vertices)]
                            if any(abs.(test_x) .< eps)
                                continue
                            end
                            vertices = hcat(vertices, x)
                        end
                    end
                end
            end
        end

        # define all facet geometry (primarily for plotting):
        facets = Vector{Matrix{Float64}}(undef, length(normals[:,1]))
        for i in 1:length(normals[:,1])
            facet = Matrix{Float64}(undef, 3, 0)
            angle = Vector{Float64}(undef, 0)
            # idx = Vector{Int64}(undef, 0)
            for j in 1:length(vertices[1,:])
                plane_pt = vertices[:,j][:] - faces[i,:][:]
                if abs(LinearAlgebra.dot(faces[i,:][:], plane_pt)/LinearAlgebra.norm(faces[i,:][:])/LinearAlgebra.norm(plane_pt)) <= eps
                    facet = hcat(facet, vertices[:,j][:])
                    fo = facet[:,1] - faces[i,:][:]
                    fi = facet[:,end] - faces[i,:][:]
                    # this sgn was supposed to flip surface normal on facets, but seems to do nothing.
                    sgn = sign(dot(faces[i,:][:], normals[i,:]))
                    x = fo
                    z = sgn*normals[i,:]
                    # z = normals[i,:]
                    y = LinearAlgebra.cross(z,x)
                    xp = fi'*x
                    yp = fi'*y
                    this_angle = atan(yp,xp)
                    push!(angle, this_angle)
                    # push!(idx, j)
                end
            end
            facet = facet[:,sortperm(angle)]
            # add last vertex to facet again, for drawing.
            facet = hcat(facet, facet[:,1])
            facets[i] = facet
        end
        return new{Float64}(normals, distances, vertices, facets)
    end
end

function project(u::Vector{<:Real}, p::Polyhedron{<:Real})
    # condition the input based on the size of the Polyhedron
    Lmin = min(p.distances...)
    Lmax = max(p.distances...)
    v = u/norm(u)*Lmin
    eps = Lmin/1000

    model = JuMP.Model(MadNLP.Optimizer)
    JuMP.set_silent(model)
    # JuMP.set_attribute(model, "algorithm", :LD_MMA)
    # scalar solution in gamma:
    JuMP.@variable(model, gamma >= eps)
    JuMP.@constraint(model, gamma*p.normals*v <= p.distances)
    JuMP.@objective(model, Max, gamma^2)
    
    # linear solution in w:
    # JuMP.@variable(model, w[1:3] .>= eps)
    # JuMP.@constraint(model, p.normals*w <= p.distances)
    # JuMP.@objective(model, Max, v'*w)

    # linear solution in w (SDP):
    # JuMP.@variable(model, w[1:3])
    # JuMP.@constraint(model, p.normals*w <= p.distances)
    # JuMP.@constraint(model, (w[1]*v[1] + w[2]*v[2] + w[3]*v[3])/sqrt(w[1]^2 + w[2]^2 + w[3]^2)/norm(v) == 1.0)
    # JuMP.@objective(model, Max, w[1]^2 + w[2]^2 + w[3]^2)

    JuMP.optimize!(model)
    if JuMP.is_solved_and_feasible(model)
        return JuMP.value.(gamma)*v
        # return JuMP.value.(w)
    else
        error("failed to project point ", u, " onto polyhedron!")
    end
end

function write(file::String, p::Polyhedron{<:Real})
    suffix = splitext(file)[2]
    res_str = ""
    if lowercase(suffix) == ".obj"
        error(".obj format not supported!")
        res_str = gen_obj(p)
    elseif lowercase(suffix) == ".stl"
        res_str = gen_stl(p)
    end

    open(file, "w") do file
        Base.write(file, res_str)
    end    
end

function gen_obj(p::Polyhedron)::String
    vs = ""
    for (k, facet) in enumerate(p.facets)
        for i in 1:length(facet[1,:])
            vs *= "v "
            for j in 1:3
                vs *= @sprintf "%0.6f " facet[j,i]
            end
            vs *= "\n"
        end
    end
    println(vs)
    return vs
end

function gen_stl(p::Polyhedron)::String
    name = "P"
    vs = "solid "
    vs *= name
    vs *= "\n"
    for (k, facet) in enumerate(p.facets)
        vs *= "facet normal "
        n = cross(facet[:,1])*facet[:,2]
        for i in 1:3
            vs *= @sprintf "%e " n[i]
        end
        vs *= "\n\touter loop\n"
        for j in 1:length(facet[1,:])
            vs *= "\t\tvertex " 
            for i in 1:3
                vs *= @sprintf "%e " facet[i,j]
            end
            vs *= "\n"
        end
        vs *= "\tendloop\n"
        vs *= "endfacet\n"
    end
    vs *= "endsolid "
    vs *= name
    # println(vs)
    return vs
end