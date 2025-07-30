using OpenSOAP
using GLMakie, GeometryBasics
import LinearAlgebra
import JuMP, Ipopt, MathOptInterface

function wheel_momenta()
# todo: move the wheel config into yaml.

    # define wheel array:
    ang = 26.57*pi/180
    w1 = r_euler2(ang)'*[1;0;0]
    w2 = r_euler1(-ang)'*[0;1;0]
    w3 = r_euler2(-ang)'*[-1;0;0]
    w4 = r_euler1(ang)'*[0;-1;0]
    # w5 = [0;0;1]
    nwheels = 4

    W = [w1 w2 w3 w4]       # wheel axes
    h = 2e-3*ones(nwheels)    # momentum at 6000 RPM
    h[3] = h[3]

    # check all pairs of wheels:
    #   define half-spaces (planes) constraining momentum:
    normals = zeros(nwheels*(nwheels - 1),3)
    distances = zeros(nwheels*(nwheels - 1))
    faces = zeros(nwheels*(nwheels - 1),3)
    v = 1
    for i in 1:nwheels
        for j in 1:nwheels
            if i == j
                continue
            end
            # note: these `n` define the normals to bounding planes 
            n = cross(W[:,i])*W[:,j]
            n = n/LinearAlgebra.norm(n)
            
            sgns = [sign(W[:,k]'*n) for k in 1:nwheels]

            distances[v] = (W*(h.*sgns))'*n
            normals[v,:] = n
            faces[v,:] = distances[v].*normals[v,:] 
            v += 1
        end
    end

    L = min(distances...)
    D = max(distances...)
    eps = L/1000
    
    # find all vertices of bounding polyhedron:
    nconstraints = v - 1
    vertices = Matrix{Float64}(undef, 3, 0)
    for i in 1:nconstraints
        for j in i:nconstraints
            for k in j:nconstraints
                if i == j || j == k || i == k
                    continue
                end
                P = normals[[i,j,k],:]
                q = distances[[i,j,k]]
                if LinearAlgebra.rank(P) == 3
                    # check if vertex is solution to LMI for polytope surface
                    #   this check excludes all facet intersections which are external to the polytope.
                    x = P\q
                    if all(normals*x .<= (distances .+ eps))
                        # check if vertex already in list:
                        test_x = [x'*vertex/LinearAlgebra.norm(vertex)/LinearAlgebra.norm(x) .- 1 for vertex in eachcol(vertices)]
                        if any(abs.(test_x) .< eps)
                            continue
                        end
                        vertices = hcat(vertices, x)
                    end
                end
            end
        end
    end

    println("vertices: ", length(vertices[1,:]))
    println("faces: ", length(normals[:,1]))

    # define all facet geometry (primarily for plotting):
    facets = Vector{Matrix{Float64}}(undef, length(normals[:,1]))
    facet_list = Matrix{Float64}(undef, 4, 0)
    for i in 1:length(normals[:,1])
        facet = Matrix{Float64}(undef,3,0)
        angle = Vector{Float64}(undef, 0)
        idx = Vector{Float64}(undef, 0)
        for j in 1:length(vertices[1,:])
            plane_pt = vertices[:,j][:] - faces[i,:][:]
            if abs(LinearAlgebra.dot(faces[i,:][:], plane_pt)/LinearAlgebra.norm(faces[i,:][:])/LinearAlgebra.norm(plane_pt)) <= eps
                facet = hcat(facet, vertices[:,j][:])
                fo = facet[:,1] - faces[i,:][:]
                fi = facet[:,end] - faces[i,:][:]
                x = fo
                z = normals[i,:]
                y = LinearAlgebra.cross(z,x)
                xp = fi'*x
                yp = fi'*y
                this_angle = atan(yp,xp)
                push!(angle, this_angle)
                push!(idx, j)
            end
        end
        facet = facet[:,sortperm(angle)]
        # add last vertex to facet again, for drawing.
        facet = hcat(facet, facet[:,1])
        facets[i] = facet
        
        # println("size of facet: ", size(facet))
        # println("length of idx: ", size(idx))
        facet_list = hcat(facet_list, idx)
    end

    # solve QP for nearest point to polytope:
    max_constraint = distances .+ eps
    min_constraint = distances .- eps
    npts = 1000  # number of points to solve (project onto hull)
    sources = rand(3,npts)
    pts = zeros(3,npts)
    
    # rotate each vector into random orientation (to cover the sphere)
    [sources[:,k] = r_random()*sources[:,k]/LinearAlgebra.norm(sources[:,k]) for k in 1:npts]

    # scale each vector to the inscribed sphere of the polyhedron:
    sources = 1.0 .* sources .* min(distances...)

    println("solving...")
    @time begin
        for k in 1:npts
            model = JuMP.Model(Ipopt.Optimizer)
            JuMP.set_silent(model)
            # scalar version:
            JuMP.@variable(model, gamma >= eps)
            JuMP.@constraint(model, gamma*normals*sources[:,k] <= distances)
            JuMP.@objective(model, Max, gamma^2)
            JuMP.optimize!(model)
            if JuMP.is_solved_and_feasible(model)
                pts[:,k] = JuMP.value.(gamma)*sources[:,k]
            else
                println("failed to solve $k !")
            end
        end
    end

    # plotting:
    GLMakie.activate!()

    ambient_brightness = 0.8
    al = AmbientLight(RGBf(ambient_brightness))
    fig = Figure(size=(800,600))
    ax = LScene(
        fig[1,1], 
        show_axis=false, 
        # height=400,
        # width=400,
        scenekw=(
            lights=[al],
            clear=true
        )
    )
    scatter!(
        ax,
        faces[:,1],
        faces[:,2],
        faces[:,3],
        color=:black
    )

    drawaxes = false
    if drawaxes
        colors = [:red, :green, :blue]
        axes = LinearAlgebra.diagm([1;1;1])
        for i in 1:3
            arrows!(
                ax,
                [Point3f(0)],
                [Vec3f(axes[1,i],axes[2,i],axes[3,i])],
                linewidth=0.02,
                arrowsize=Vec3f(0.02, 0.02, 0.1),
                color=colors[i]
            )
        end
    end
    # draw origin:
    scatter!(
        ax,
        0,0,0,
        color=:blue,
        marker='+',
        markersize=30
    )

    scatter!(
        ax,
        sources[1,:],
        sources[2,:],
        sources[3,:],
        color=:lightgreen,
        alpha=0.3
    )
    scatter!(
        ax,
        pts[1,:],
        pts[2,:],
        pts[3,:],
        color=:blue,
        markersize=3
    )
    scatter!(
        ax,
        vertices[1,:],
        vertices[2,:],
        vertices[3,:],
        color=:red
    )

    for (i,facet) in enumerate(facets)
        lines!(
            ax,
            facet[1,:],
            facet[2,:],
            facet[3,:],
            # colorrange=(1,length(facets)),
            # color=i,
            # colormap=:magma
            color=:black
        )
        mesh!(
            ax,
            facet[1,1:3],
            facet[2,1:3],
            facet[3,1:3],
            color=:black,
            alpha=0.3
        )
        mesh!(
            ax,
            facet[1,3:5],
            facet[2,3:5],
            facet[3,3:5],
            color=:black,
            alpha=0.3
        )
    end

    cam3d!(
        ax,
        upvector=Vec3f(0,0,1)
    )
    # cam3d!(
    #     ax_polytope,
    #     upvector=Vec3f(0,0,1),
    #     lookat=Vec3f(pts[1,1],pts[2,1],pts[3,1])
    # )
    display(fig)
end

function check_polyhedron()
    sim = load_mission("config/mission.yaml")

    P_t = sim.mission.spacecraft.attitude.wheels.torque_env
    P_m = sim.mission.spacecraft.attitude.wheels.momentum_env

    npts = 1000  # number of points to solve (project onto hull)
    sources = rand(3,npts)
    pts = zeros(3,npts)
    
    # rotate each vector into random orientation (to cover the sphere)
    for i in 1:npts
        sources[:,i] = r_random()*sources[:,i]/LinearAlgebra.norm(sources[:,i]) .* min(P_t.distances...) * 2.0
    end
    # sources = sources .* min(P.distances...) .* 2

    @time for i in 1:npts
        pts[:,i] = project(sources[:,i], P_t)
    end

    GLMakie.activate!()

    ambient_brightness = 0.8
    al = AmbientLight(RGBf(ambient_brightness))
    pl = PointLight(RGBf(ambient_brightness), Point3f(10.0,0,0))
    fig = Figure(size=(800,600))
    ax = LScene(
        fig[1,1], 
        show_axis=false, 
        # height=400,
        # width=400,
        scenekw=(
            lights=[al, pl],
            clear=true
        )
    )

    polyplot!(P_t)
    # polyplot!(P_m)

    OpenSOAP.write("test.stl", P_t)

    scatter!(
        ax,
        pts[1,:],
        pts[2,:],
        pts[3,:],
        color=:red,
        markersize=3
    )

    cam3d!(
        ax,
        upvector=Vec3f(0,0,1)
    )

    display(fig)
end

function slew()
    sim = load_mission("config/slew_test.yaml")
    
    C_test = [
        0.644037   0.192398  -0.740405;
        0.697663   0.249319   0.671645;
        0.313821  -0.949117   0.0263418
    ]

    # define a maneuver
    maneuver = Maneuver(
        sim.tspan[end] - (sim.tspan[end] - sim.tspan[1])/6,
        (sim.tspan[end] - sim.tspan[1])/2,
        r_random()
        # C_test
    )

    println("Maneuver: ")
    println("\tSim start: ", sim.tspan[1])
    println("\tManeuver start: ", maneuver.tf - maneuver.dt)
    println("\tManeuver end: ", maneuver.tf)
    println("\tSim end: ", sim.tspan[end])
    initial = State{Float64}(sim)

    println("Integrating...")
    soln = integrate_system(dynamics_attitude!, initial, sim.tspan, sim.dt, maneuver, sim)

    GLMakie.activate!()
    fig = Figure(size=(860,860))
    set_theme!(theme_light())
    plot_state!(fig, maneuver, soln, sim)

    display(fig)
    # call integrate_system()...
    # back out attitude trajectory and plot.

end

function vectors_frames()
    GLMakie.activate!()
    fig = Figure(size=(400,400))
    set_theme!(theme_light())

    ax = LScene(
        fig[1, 1],
        show_axis=false
    )
        
    basis = diagm([1.0,1.0,1.0])
    cols = [:red, :green, :blue]
    lscale = 0.5
    wscale = 0.02
    C1 = r_random()
    C2 = r_random()
    C3 = basis

    r1 = [1,0,1]
    r2 = [0,1,1]
    r3 = zeros(3)

    Cs = (C1, C2, C3)
    rs = (r1, r2, r3)
    for d in zip(Cs, rs)
        C = d[1]
        r = d[2]

        tails = [Point3f(r) for k in 1:3]
        heads = [Vec3f(lscale*C*basis[:,k]) for k in 1:3]
        for k in 1:3
            arrows!(
                tails,
                heads,
                color=cols,
                linewidth = wscale,
                arrowsize = Vec3f(wscale,wscale,wscale)
            )
        end
    end

    display(fig)
end