using GLMakie, Colors, GeometryBasics, LinearAlgebra, Rotations
using SatelliteToolboxBase, SatelliteToolboxTransformations
using FileIO

function lookat_orbit!(ax::Makie.LScene, t_jd_s, soln::Dict)
    dist_scale = 4
    point = lift(t_jd_s) do t_jd_s
        (_,i) = findmin(x->abs(x - t_jd_s), soln["time"])
        orbit_normal = LinearAlgebra.cross(soln["state"][1:3, i], soln["state"][4:6, i])
        return Vec3d(orbit_normal/norm(orbit_normal) * norm(soln["state"][1:3, i]) * dist_scale)
    end

    cam3d!(
        ax, 
        projectiontype="Orthographic", 
        lookat=Vec3d(0), 
        upvector=Vec3d(0,0,1), 
        eyeposition=point
    )
end

function hud_conop!(ax::Makie.LScene, t_jd_s, target_histories, targets, soln, texture_ecef, eops)
    r_E = EARTH_EQUATORIAL_RADIUS
    r_P = EARTH_POLAR_RADIUS
    
    r0m = norm(soln["state"][1:3,1])
    orbit_period = 2*π*sqrt(r0m^3/SatelliteToolboxBase.GM_EARTH)

    sun_targets = [target for target in targets if typeof(target) === SunTarget]
    ground_targets = [target for target in targets if typeof(target) === GroundTarget]
    
    dist_scale = 2

    (_,i_start) = findmin(x->abs(x - t_jd_s), soln["time"])
    (_,i_end) = findmin(x->abs(x - t_jd_s - orbit_period), soln["time"])
    span = i_start:i_end
    
    C_IF_r = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, eops)
    C_IF = Matrix(C_IF_r)

    gs_trace = []
    for target in ground_targets
        trace = Vector{Point3f}(undef, length(span))
        for (ti, si) in enumerate(span)
            trace[ti] = Point3f(position_eci(target, soln["time"][si]))
        end
        push!(gs_trace, trace)
    end
        

    states = soln["state"][21, span[1]:span[end]]
    start_state = span[1]
    orbit_normal = LinearAlgebra.cross(soln["state"][1:3, span[1]], soln["state"][4:6, span[1]])
    eye_pos = Vec3d(orbit_normal/norm(orbit_normal) * norm(soln["state"][1:3, i_start]) * dist_scale)
    start_r = soln["state"][1:3, span[1]]

    δstate_i = findall(diff(states) .!= 0)
    push!(δstate_i, length(span))

    span_per_slice = Vector{Tuple{Int64, Int64}}()
    for i in eachindex(δstate_i)
        if i > 1
            push!(span_per_slice, (δstate_i[i - 1], δstate_i[i] - 1))
        else
            push!(span_per_slice, (1, δstate_i[i] - 1))
        end
    end

    # println("length of time:\t\t\t", length(soln["time"]))
    # println("absolute span:\t\t\t", span)
    # println("length of absolute span:\t", length(span))
    # println("δstate_i:\t\t\t", δstate_i)
    # println("span_per_slice:\t\t\t", span_per_slice)
    
    # making earth surface texture:
    nθ = length(texture_ecef[:,1])
    nφ = 100
    θ = range(0, stop=2π, length=nθ)
    φ = range(0, stop=π, length=nφ)
    
    mesh_X = zeros(nθ, nφ)
    mesh_Y = zeros(nθ, nφ)
    mesh_Z = zeros(nθ, nφ)

    local_r = r_E*cos(35*pi/180)
    local_z = r_P*sin(35*pi/180)
    p35 = [Point3f(C_IF*[local_r * cos(ang); local_r * sin(ang); local_z]) for ang in θ]
    n35 = [Point3f(C_IF*[local_r * cos(ang); local_r * sin(ang); -local_z]) for ang in θ]
    equat = [Point3f(C_IF*[r_E * cos(ang); r_E * sin(ang); 0]) for ang in θ]
    for i in 1:nθ
        for j in 1:nφ
            # a point in the mesh
            p = C_IF*[r_E*cos(θ[i])*sin(φ[j]); r_E*sin(θ[i])*sin(φ[j]); r_P*cos(φ[j])]
            mesh_X[i,j] = p[1]
            mesh_Y[i,j] = p[2]
            mesh_Z[i,j] = p[3]
        end
    end

    # build slice meshes
    slice_mesh = Dict()
    r_inner = norm(start_r)*1.1
    bandwidth = norm(start_r)*3

    for (k, this_span) in enumerate(span_per_slice)
        # println("slice: ", this_span)
        this_span_range = this_span[1]:this_span[2]

        nθ = length(this_span_range)
        nr = 2
        # rs = LinRange(k*r_inner, k*(r_inner + bandwidth), nr)
        rs = LinRange(r_inner, (r_inner + bandwidth), nr)
        X = zeros(nθ, nr)
        Y = zeros(nθ, nr)
        Z = zeros(nθ, nr)

        c = Int64(states[this_span_range[1] + 1])

        for i in 1:nθ
            this_r = norm(soln["state"][1:3, i])
            for j in 1:nr
                X[i,j] = soln["state"][1, span[this_span_range[i]]]/this_r*rs[j]
                Y[i,j] = soln["state"][2, span[this_span_range[i]]]/this_r*rs[j]
                Z[i,j] = soln["state"][3, span[this_span_range[i]]]/this_r*rs[j]
            end
        end
        slice_mesh[k] = (X, Y, Z, c)
    end

    # draw earth
    surface!(
        ax,
        mesh_X,
        mesh_Y,
        mesh_Z,
        color=texture_ecef,
        diffuse=0.8
    )
    # orbit
    taillength = Int64(floor(length(span) / 2)) + 1
    taillength = length(span) - 100
    local_range = range(length(span) - taillength, length(span))
    n = length(local_range)
    lines!(
        ax,
        soln["state"][1,span[local_range]],
        soln["state"][2,span[local_range]],
        soln["state"][3,span[local_range]],
        # color=1:length(soln["state"][1, span[1:taillength]]),
        # colormap=Reverse(:inferno)
        color=[RGBAf(1, 1, 0, i/n) for i in 1:n]
    )
    # spacecraft position
    scatter!(
        ax,
        Point3f(soln["state"][1:3,i_end]),
        color=:yellow,
        markersize=15
    )
    # groundstations
    for trace in gs_trace
        # lines!(
        #     ax,
        #     trace,
        #     color=:lightgreen
        # )
        scatter!(
            ax,
            trace[end],
            color=:lightgreen,
            markersize=20
        )
    end
    # scatter!(
    #     ax,
    #     gs_pts,
    #     color=:lightgreen,
    #     markersize=20
    # )
    # 35º latitude lines:
    lines!(
        ax,
        p35,
        color=:white
    )
    lines!(
        ax,
        n35,
        color=:white,
        label="35º latitude",
    )
    lines!(
        ax,
        equat,
        color=:white,
        linestyle=:dashdot
    )
    text!(
        ax,
        p35[1 + Int64(round(length(p35)/2))],
        text="+35º",
        color=:white,
        depth_shift=1,
        overdraw=true
    )
    text!(
        ax,
        n35[1 + Int64(round(length(p35)/2))],
        text="-35º",
        color=:white,
        depth_shift=1,
        overdraw=true
    )

    colors = [:blue, :orange, :lightgreen, :violet]
    names = ["safe", "power", "downlink", "science"]
    for slice in keys(slice_mesh)
        surface!(
            ax,
            slice_mesh[slice][1],
            slice_mesh[slice][2],
            slice_mesh[slice][3],
            color=fill((colors[slice_mesh[slice][4]], 0.8), 20, 2),
            invert_normals=true,
            fxaa=true,
        )
    end
    # place camera:
    cam3d!(
        ax, 
        projectiontype="Orthographic", 
        lookat=Vec3d(0), 
        upvector=Vec3d(0,0,1), 
        eyeposition=eye_pos
    )
    zoom!(
        ax.scene,
        cameracontrols(ax.scene),
        1
    )
    update_cam!(ax.scene, cameracontrols(ax.scene))
end

function plot_meta!(ax::Makie.LScene, t_jd_s, enable, target_histories::Dict, soln::Dict, eops)
    r_E = EARTH_EQUATORIAL_RADIUS
    r0m = norm(soln["state"][1:3,1])
    r_inner = 1.2*r0m
    bandwidth = 0.2*r0m
    
    orbit_period = 2*π*sqrt(r0m^3/SatelliteToolboxBase.GM_EARTH)
    
    dist_scale = 2

    meshes = lift(t_jd_s) do t_jd_s
        (_,i_start) = findmin(x->abs(x - t_jd_s), soln["time"])
        (_,i_end) = findmin(x->abs(x - t_jd_s - orbit_period), soln["time"])
        span = i_start:i_end
        
        C_IF_r = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, eops)
        C_IF = Matrix(C_IF_r)

        states = soln["state"][21, span[1]:span[end]]
        start_state = span[1]
        orbit_normal = LinearAlgebra.cross(soln["state"][1:3, span[1]], soln["state"][4:6, span[1]])
        eye_pos = Vec3d(orbit_normal/norm(orbit_normal) * norm(soln["state"][1:3, i_start]) * dist_scale)
        start_r = soln["state"][1:3, span[1]]

        δstate_i = findall(diff(states) .!= 0)
        push!(δstate_i, length(span))

        span_per_slice = Vector{Tuple{Int64, Int64}}()
        for i in eachindex(δstate_i)
            if i > 1
                push!(span_per_slice, (δstate_i[i - 1], δstate_i[i] - 1))
            else
                push!(span_per_slice, (1, δstate_i[i] - 1))
            end
        end
        
        # build slice meshes
        slice_mesh = Dict()
        r_inner = norm(start_r)*1.1
        bandwidth = norm(start_r)*3

        for (k, this_span) in enumerate(span_per_slice)
            # println("slice: ", this_span)
            this_span_range = this_span[1]:this_span[2]

            nθ = length(this_span_range)
            nr = 2
            # rs = LinRange(k*r_inner, k*(r_inner + bandwidth), nr)
            rs = LinRange(r_inner, (r_inner + bandwidth), nr)
            X = zeros(nθ, nr)
            Y = zeros(nθ, nr)
            Z = zeros(nθ, nr)

            c = Int64(states[this_span_range[1] + 1])

            for i in 1:nθ
                this_r = norm(soln["state"][1:3, i])
                for j in 1:nr
                    X[i,j] = soln["state"][1, span[this_span_range[i]]]/this_r*rs[j]
                    Y[i,j] = soln["state"][2, span[this_span_range[i]]]/this_r*rs[j]
                    Z[i,j] = soln["state"][3, span[this_span_range[i]]]/this_r*rs[j]
                end
            end
            slice_mesh[k] = (X, Y, Z, c)
        end
        return slice_mesh
    end

    pre_mesh = Dict{Int64, Observable{Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Int64}}}()
    # plot_mesh_x = Dict{Int64, Observable{Matrix{Float64}}}()
    # plot_mesh_y = Dict{Int64, Observable{Matrix{Float64}}}()
    # plot_mesh_z = Dict{Int64, Observable{Matrix{Float64}}}()
    # plot_mesh_c = Dict{Int64, Observable{Int64}}()

    # outer index: slice
    # inner index: coord (x,y,z)
    plot_mesh_x = Vector{Observable{Matrix{Float64}}}()
    plot_mesh_y = Vector{Observable{Matrix{Float64}}}()
    plot_mesh_z = Vector{Observable{Matrix{Float64}}}()

    k = 1
    more_slice = true
    while more_slice
        temp = lift(meshes) do meshes
            if k in keys(meshes)
                k += 1
                return meshes[k - 1][1]
            else
                more_slice = false
                println("no key!")
            end
        end
        if more_slice
            push!(plot_mesh_x, temp) 
        end
    end
    k = 1
    more_slice = true
    while more_slice
        temp = lift(meshes) do meshes
            if k in keys(meshes)
                k += 1
                return meshes[k - 1][2]
            else
                more_slice = false
                println("no key!")
            end
        end
        if more_slice
            push!(plot_mesh_y, temp) 
        end
    end
    k = 1
    more_slice = true
    while more_slice
        temp = lift(meshes) do meshes
            if k in keys(meshes)
                k += 1
                return meshes[k - 1][3]
            else
                more_slice = false
                println("no key!")
            end
        end
        if more_slice
            push!(plot_mesh_z, temp) 
        end
    end

    for j in 1:4
        surface!(
            ax,
            plot_mesh_x[j],
            plot_mesh_y[j],
            plot_mesh_z[j]
        )
    end
    # n_mesh = lift(meshes) do meshes
    #     return (meshes, length(meshes))
    # end

    # colors = [:blue, :orange, :green, :purple]

    # for j in 1:(n_mesh[][2])
    #     surface!(
    #         ax,
    #         meshes[j][1],
    #         meshes[j][2],
    #         meshes[j][3],
    #         color=fill((colors[meshes[j][4]], 0.9), 20, 2)
    #     )
    # end

    # result = lift(meshes) do mesh
    #     for j in eachindex(mesh)
    #         surface!(
    #             ax,
    #             mesh[j][1],
    #             mesh[j][2],
    #             mesh[j][3],
    #             color=fill((colors[mesh[j][4]], 0.9), 20, 2)
    #         )
    #     end
    # end
        
    # for j in 1:n_mesh[]
    #     surface!(
    #         ax,
    #         meshes[j][1],
    #         meshes[j][2],
    #         meshes[j][3],
    #         # color=:lightgreen,
    #         # color=fill((colors[meshes[4][j]],0.8),20,2)
    #     )
    # end
end

function plot_mode!(ax::Makie.Axis, t_jd_s, target_histories::Dict, soln::Dict)
    playhead = lift(t_jd_s) do t_jd_s
        t_0 = soln["time"][1]
        return t_jd_s - t_0
    end
    
    names = ["safe", "power", "downlink", "science"]
    n_modes = Int64(max(soln["state"][21,:]...) - min(soln["state"][21,:]...))
    m_hist = zeros(Int64, 1+n_modes, length(soln["state"][21,:]))
    for i = 1:length(soln["time"])
        m_hist[Int64(soln["state"][21,i]),i] = 1
    end
    t_0 = soln["time"][1]

    for i in 1:n_modes+1
        change = diff(m_hist[i,:])
        starts = findall(i->(i > 0), change)
        stops = findall(i->(i < 0), change)
        if m_hist[i,1] == 1
            # starts = [1; starts]
        end
        if stops[1] == 1
            stops = stops[2:end]
        end
        
        if length(starts) == 0 && length(stops) == 0
            println("found no contacts for key!")
            continue
        end

        plotval = [Rect(soln["time"][starts[j]] - t_0, 0, soln["time"][stops[j] - 1] - soln["time"][starts[j]], 1) for j in 1:min(length(starts), length(stops))]

        poly!(
            ax,
            plotval,
            alpha=1.0,
            strokewidth=0.1,
            # strokecolor=:white,
            label=names[i]
        )
    end
    
    vlines!(
        ax,
        playhead,
        color=:yellow,
        alpha=0.7,
        linewidth=1
    )
end

function plot_power!(ax::Makie.Axis, t_jd_s, target_histories::Dict, soln::Dict)
    playhead = lift(t_jd_s) do t_jd_s
        t_0 = soln["time"][1]
        return t_jd_s - t_0
    end
    
    e_hist = soln["state"][19,:]
    
    lw = 2.0
    lines!(
        ax, 
        soln["time"] .- soln["time"][1],
        e_hist/3600,
        color=:yellow,
        linewidth=lw
    )

    vlines!(
        ax,
        playhead,
        color=:yellow,
        alpha=0.7,
        linewidth=1
    )
end

function plot_data!(ax::Makie.Axis, t_jd_s, target_histories::Dict, soln::Dict)
    playhead = lift(t_jd_s) do t_jd_s
        t_0 = soln["time"][1]
        return t_jd_s - t_0
    end
    
    e_hist = soln["state"][20,:]
    
    lw = 2.0
    lines!(
        ax, 
        soln["time"] .- soln["time"][1],
        e_hist/8e6,
        color=:cyan,
        linewidth=lw
    )

    vlines!(
        ax,
        playhead,
        color=:yellow,
        alpha=0.7,
        linewidth=1
    )
end

function plot_visibilities!(ax::Makie.Axis, t_jd_s, target_histories::Dict, soln::Dict)

    playhead = lift(t_jd_s) do t_jd_s
        t_0 = soln["time"][1]
        return t_jd_s - t_0
    end

    li = 1
    for (key, val) in target_histories
        t_0 = soln["time"][1]

        # println(key)
        change = diff(val)
        starts = findall(i->(i > 0), change)
        stops = findall(i->(i < 0), change)
        if val[1] == 1
            starts = [1; starts]
        end 
        if length(starts) == 0 && length(stops) == 0
            println("found no contacts for key!")
            continue
        end
        # println(starts)
        # println(stops)
        plotval = [Rect(soln["time"][starts[i]] - t_0, li - 1, soln["time"][stops[i] - 1] - soln["time"][starts[i]], 1) for i in 1:min(length(starts), length(stops))]

        poly!(
            ax,
            plotval,
            alpha=1.0,
            strokewidth=0.1,
            strokecolor=:white,
            label=key
        )

        li += 1
    end
    vlines!(
        ax,
        playhead,
        color=:yellow,
        alpha=0.7,
        linewidth=1
    )
end

function plot_detail!(ax::Makie.Axis, t_jd_s, soln::Dict)

    x_hist = soln["state"][7,:]
    y_hist = soln["state"][8,:]
    z_hist = soln["state"][9,:]

    playhead = lift(t_jd_s) do t_jd_s
        t_0 = soln["time"][1]
        return t_jd_s - t_0
    end

    lw = 1.0
    lines!(
        ax, 
        soln["time"] .- soln["time"][1],
        x_hist,
        color=:red,
        linewidth=lw
    )
    lines!(
        ax, 
        soln["time"] .- soln["time"][1],
        y_hist,
        color=:green,
        linewidth=lw
    )
    lines!(
        ax, 
        soln["time"] .- soln["time"][1],
        z_hist,
        color=:blue,
        linewidth=lw
    )
    vlines!(
        ax,
        playhead,
        color=:yellow,
        alpha=0.7,
        linewidth=1
    )
end

function plot_spacecraft!(ax::Makie.LScene, t_jd_s, tail_length::Integer, soln::Dict)

    this_i = lift(t_jd_s) do t_jd_s
        (_,i) = findmin(x->abs(x - t_jd_s), soln["time"])
        return i
    end
    pos = lift(this_i) do this_i
        return Point3f(soln["state"][1:3,this_i])
    end
    pos_hist = lift(this_i) do this_i
        first_i = this_i - tail_length
        start_i = max(first_i, 1)
        return Point3f[(soln["state"][1:3,j]) for j in start_i:this_i]
    end
    
    scatter!(
        ax,
        pos,
        color=:yellow,
        markersize=10
    )
    lines!(
        ax,
        pos_hist,
        color=:yellow,
        linewidth=0.5,
        fxaa=true
    )
    
end

function plot_targets!(ax::Makie.LScene, targets::Vector{<:AbstractTarget}, t_jd_s, soln, eops)
    r_E = EARTH_EQUATORIAL_RADIUS
    arrow_scale = 0.01
    axis_scale = 2

    sun_targets = [target for target in targets if typeof(target) === SunTarget]
    ground_targets = [target for target in targets if typeof(target) === GroundTarget]

    gs_pts = lift(t_jd_s) do t_jd_s
        # this is a hack---need to check that the target actually is in ECEF.
        return [Point3f(position_eci(target, t_jd_s)) for target in ground_targets]
    end

    # find a much less hacky way to do this:
    sun_pos = lift(t_jd_s) do t_jd_s
        sun_I = position_eci(sun_targets[1], t_jd_s)
        sun_proj = axis_scale*r_E*sun_I/norm(sun_I)
        return [Vec3f(sun_proj)]
    end
    sun_root = lift(t_jd_s) do t_jd_s
        return [Point3f([0;0;0])]
    end

    circles = lift(t_jd_s) do t_jd_s
        C_IF_r = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, eops)
        C_IF = Matrix(C_IF_r)

        # assume circular orbit:
        r_S = norm(soln["state"][1:3,1])
        

        nθ = 24
        θ = range(0, stop=2π, length=nθ)
        X = zeros(nθ, length(ground_targets))
        Y = zeros(nθ, length(ground_targets))
        Z = zeros(nθ, length(ground_targets))
        for k in eachindex(ground_targets)
            # first draw cone z-up
            # then transform to align with zenith-up at target.position
            target = ground_targets[k]
            α = π - target.cone
            β = π - target.cone - π/2
            γ = asin(r_E/r_S*sin(α))
            r_circle = r_S*cos(γ + β)   # radius of gs cone's circle projected onto sphere of orbit
            δ = r_circle*tan(β)

            pos_F = Vector(position_ecef(target, t_jd_s))
            C_FT = r_min_arc([0;0;1], pos_F / norm(pos_F))

            c_F = pos_F + (1+δ)*pos_F/norm(pos_F)  # center of gs cone's circle projected onto sphere of orbit

            for i in 1:nθ
                p_T = r_circle .* [
                    cos(θ[i]);
                    sin(θ[i]);
                    0
                ]
                p = C_IF*(C_FT*p_T) + C_IF*c_F
                # p = C_IF*(C_FT*p_T) + position_eci(target, t_jd_s)
                X[i,k] = p[1]
                Y[i,k] = p[2]
                Z[i,k] = p[3]
            end
        end
        return (X, Y, Z)
    end

    # mesh for whole cone
    # cone_meshes = lift(t_jd_s) do t_jd_s
    #     C_IF_r = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, eops)
    #     C_IF = Matrix(C_IF_r)

    #     height = 0.05

    #     nθ = 20
    #     nζ = 2
    #     θ = range(0, stop=2π, length=nθ)
    #     ζ = range(0, stop=r_E*height, length=nζ)
        
    #     mesh_X = zeros(nθ, nζ, length(ground_targets))
    #     mesh_Y = zeros(nθ, nζ, length(ground_targets))
    #     mesh_Z = zeros(nθ, nζ, length(ground_targets))
    #     # for (k, target) in ground_targets
    #     for k in eachindex(ground_targets)
    #         # first draw cone z-up
    #         # then transform to align with zenith-up at target.position
    #         target = ground_targets[k]
    #         b = tan(target.cone)
    #         pos_F = Vector(position_ecef(target, t_jd_s))
    #         C_FT = r_min_arc([0;0;1], pos_F / norm(pos_F))

    #         # note: currently, this ignores FixedFrameTarget.direction
    #         for i in 1:nθ
    #             for j in 1:nζ
    #                 # a point in the mesh
    #                 p_T = [
    #                     b*cos(θ[i])*ζ[j]; 
    #                     b*sin(θ[i])*ζ[j]; 
    #                     ζ[j]
    #                 ]
    #                 p = C_IF*(C_FT*p_T) + position_eci(target, t_jd_s)
    #                 mesh_X[i,j,k] = p[1]
    #                 mesh_Y[i,j,k] = p[2]
    #                 mesh_Z[i,j,k] = p[3]
    #             end
    #         end
    #     end
    #     return (mesh_X, mesh_Y, mesh_Z)
    # end

    # meshs = Matrix{Observable{Matrix{Float64}}}(undef, 3,length(targets))
    # for i in 1:3
    #     for j in 1:length(ground_targets)
    #         meshs[i,j] = lift(cone_meshes) do globe_mesh
    #             return globe_mesh[i][:,:,j]
    #         end
    #     end
    # end

    gs_lines = Matrix{Observable{Vector{Float64}}}(undef, 3, length(targets))
    for i in 1:3
        for k in 1:length(ground_targets)
            gs_lines[i,k] = lift(circles) do cone_circle
                return cone_circle[i][:,k]
            end
        end
    end
    
    scatter!(
        ax,
        gs_pts,
        color=:lightgreen,
        markersize=10
    )

    for i in 1:length(ground_targets)
        # surface!(
        #     ax,
        #     meshs[1,i],
        #     meshs[2,i],
        #     meshs[3,i],
        #     # color=:lightgreen,
        #     color=fill((:lightgreen,0.33),20,2)
        # )
        lines!(
            ax,
            gs_lines[1,i],
            gs_lines[2,i],
            gs_lines[3,i],
            color=:lightgreen
        )
    end

    # for i in 1:length(sun_targets)
    #     arrows!(
    #         ax,
    #         sun_root,
    #         sun_pos,
    #         color=:yellow,
    #         linewidth=arrow_scale*r_E,
    #         arrowsize=Vec3f(arrow_scale*r_E, arrow_scale*r_E, arrow_scale*3*r_E)
    #     )
    #     break
    # end
    for i in 1:length(sun_targets)
        scatter!(
            ax,
            sun_pos,
            color=RGBf(243/255, 241/255, 218/255),
            marker='\u2609',
            markersize=20
        )
        break
    end
end

function plot_earth!(ax::Makie.LScene, t_jd_s, iers_eops, texture_ecef)
    r_E = EARTH_EQUATORIAL_RADIUS
    r_P = EARTH_POLAR_RADIUS
    
    C_IF = lift(t_jd_s) do t_jd_s
        C_IF_r = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, iers_eops)
        C_IF_m = Matrix(C_IF_r)
        # eci vector = C_IF * ecef vector 
        return C_IF_m
    end

    # making surface texture:
    # nθ = 3600
    nθ = length(texture_ecef[:,1])
    nφ = 100
    θ = range(0, stop=2π, length=nθ)
    φ = range(0, stop=π, length=nφ)

    local_r = r_E*cos(35*pi/180)
    local_z = r_P*sin(35*pi/180)

    p35 = lift(C_IF) do C_IF
        return [Point3f(C_IF*[local_r * cos(ang); local_r * sin(ang); local_z]) for ang in θ]    
    end
    n35 = lift(C_IF) do C_IF
       return [Point3f(C_IF*[local_r * cos(ang); local_r * sin(ang); -local_z]) for ang in θ]
    end
    equat = lift(C_IF) do C_IF
        return [Point3f(C_IF*[r_E * cos(ang); r_E * sin(ang); 0]) for ang in θ]
    end
    
    
    globe_mesh = lift(C_IF) do C_IF
        mesh_X = zeros(nθ, nφ)
        mesh_Y = zeros(nθ, nφ)
        mesh_Z = zeros(nθ, nφ)
        for i in 1:nθ
            for j in 1:nφ
                # a point in the mesh
                p = C_IF*[r_E*cos(θ[i])*sin(φ[j]); r_E*sin(θ[i])*sin(φ[j]); r_P*cos(φ[j])]
                mesh_X[i,j] = p[1]
                mesh_Y[i,j] = p[2]
                mesh_Z[i,j] = p[3]
            end
        end
        return (mesh_X, mesh_Y, mesh_Z)
    end
    X = lift(globe_mesh) do globe_mesh
        return globe_mesh[1]
    end
    Y = lift(globe_mesh) do globe_mesh
        return globe_mesh[2]
    end
    Z = lift(globe_mesh) do globe_mesh
        return globe_mesh[3]
    end

    surface!(
        ax,
        X,
        Y,
        Z,
        color=texture_ecef,
        diffuse=0.8
    )

    lines!(
        ax,
        p35,
        color=:white
    )
    lines!(
        ax,
        n35,
        color=:white
    )
    lines!(
        ax,
        equat,
        color=:white,
        linestyle=:dashdot
    )
end

function plot_magnetic_field!(ax::Makie.LScene, t_jd_s, soln::Dict, eops)
    r_E = EARTH_EQUATORIAL_RADIUS

    m = MagneticTarget("mag", eops)

    tail_length = 600
    step = 10
    
    n_t = length(soln["time"])
    b = zeros(3, n_t)

    this_i = lift(t_jd_s) do t_jd_s
        (_,i) = findmin(x->abs(x - t_jd_s), soln["time"])
        return i
    end
    pos = lift(this_i) do this_i
        return Point3f(soln["state"][1:3,this_i])
    end
    # pos_hist = lift(this_i) do this_i
    #     first_i = this_i - tail_length
    #     start_i = max(first_i, 1)
    #     return Point3f[(soln["state"][1:3,j]) for j in start_i:this_i]
    # end
    
    C_IF = lift(t_jd_s) do t_jd_s
        C_IF_r = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, eops)
        C_IF_m = Matrix(C_IF_r)
        # eci vector = C_IF * ecef vector 
        return C_IF_m
    end
    
    lim = 20
    lats = LinRange(-pi/2+0.1, pi/2-0.1, lim)
    lons = LinRange(-pi+0.1, pi-0.1, 2*lim)
    alts = LinRange(2*r_E, 8*r_E, 2)

    pos_I = obs_last_orbit_eci(t_jd_s, soln)

    b_all = lift(t_jd_s, pos_I) do t_jd_s, pos_I
        C_IF_r = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, eops)
        C_IF_m = Matrix(C_IF_r)

        pos_ecef = Point3f[C_IF_m'*pos for pos in pos_I]
        
        igrf_eci = Point3f[C_IF_m*(position_ecef(m, Vector(pos), t_jd_s)) for pos in pos_ecef]
        tail = Point3f[C_IF_m*pos for pos in pos_ecef]

        norms = [norm(b) for b in igrf_eci]
        head = Point3f[tip/max(norms...) for tip in igrf_eci]
        lens = [(b'*[0;0;1]) for b in head]
        lens = [p[3] for p in tail]
        return (tail, head, lens)
    end

    b_tail_all = lift(b_all) do b_all
        return b_all[1]
    end

    b_head_all = lift(b_all) do b_all
        return b_all[2]
    end

    b_lens_all = lift(b_all) do b_all
        return b_all[3]
    end

    width_scale = 0.01
    length_scale = 0.1
    arrows!(
        ax,
        b_tail_all,
        b_head_all,
        color=b_lens_all,
        colormap=:cool,
        linewidth=r_E*width_scale,
        lengthscale=r_E*length_scale,
        arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale),
        # align=:center
        # diffuse=0.8,
        # backlight=1.0
    )
end

function plot_frames!(ax::Makie.LScene, t_jd_s, soln::Dict, eops)
    r_E = EARTH_EQUATORIAL_RADIUS
    r_P = EARTH_POLAR_RADIUS

    body_arrow_scale = 0.01
    body_axis_scale = 0.1

    tail = zeros(3,3)
    body_tip = diagm([1;1;1])*r_E*body_axis_scale
    
    ecef_arrow_scale = 0.025
    ecef_axis_scale = 1.33

    tail = zeros(3,3)
    ecef_tip = diagm([1;1;1])*r_E*ecef_axis_scale
    colors = [:red,:green,:blue]

    this_i = lift(t_jd_s) do t_jd_s
        (_,i) = findmin(x->abs(x - t_jd_s), soln["time"])
        return i
    end
    pos = lift(this_i) do this_i
        return Point3f(soln["state"][1:3,this_i])
    end

    C_BI = lift(this_i) do this_i
        C_BI_m = Matrix(reshape(soln["state"][10:18, this_i], (3,3)))
        return C_BI_m
    end
    body_axes_heads = lift(C_BI) do C_BI # todo: is this wrong? should be heads = [Vec3f(soln["state"][1:3] + C_BI'*body_tip[:,i]) for i in 1:3]?
        return [Vec3f(C_BI*body_tip[:,i]) for i in 1:3]
    end
    body_axes_tails = lift(pos) do pos
        return [Point3f(pos[:,i]) for i in 1:3]
    end

    C_IF = lift(t_jd_s) do t_jd_s
        C_IF_r = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, eops)
        C_IF_m = Matrix(C_IF_r)
        # eci vector = C_IF * ecef vector 
        return C_IF_m
    end

    ecef_axes_heads = lift(C_IF) do C_IF
        return [Vec3f(C_IF*ecef_tip[:,i]) for i in 1:3]
    end

    ecef_axes_tails = lift(t_jd_s) do t_jd_s
        return [Point3f(tail[:,i]) for i in 1:3]
    end

    # drawing earth axes:
    arrows!(
        ax,
        ecef_axes_tails,
        ecef_axes_heads,
        color=colors,
        linewidth=ecef_arrow_scale*r_E,
        arrowsize=Vec3f(ecef_arrow_scale*r_E, ecef_arrow_scale*r_E, ecef_arrow_scale*3*r_E)
    )

    # drawing spacecraft body frame
    arrows!(
        ax,
        body_axes_tails,
        body_axes_heads,
        color=colors,
        linewidth=body_arrow_scale*r_E,
        arrowsize=Vec3f(body_arrow_scale*r_E, body_arrow_scale*r_E, body_arrow_scale*3*r_E)
    )
end

function load_eops()
    return SatelliteToolboxTransformations.fetch_iers_eop()
end

function load_earth_texture_to_ecef(path::String)
    texture = load(path)
    texture = texture'
    W_uv = length(texture[:,1])
    texture = circshift(texture, (round(W_uv)/2, 0))
    return texture
end

function obs_last_orbit_eci(t_jd_s::Observable{Float64}, soln::Dict)
    pos_eci = lift(t_jd_s) do t_jd_s
        step = 5
        
        (_,this_i) = findmin(x->abs(x - t_jd_s), soln["time"])
        this_pos = soln["state"][1:3,this_i]
        period = 2*pi*sqrt(norm(this_pos)^3/SatelliteToolboxBase.GM_EARTH)
        last_orbit = t_jd_s - period
        (_,last_i) = findmin(x->abs(x - last_orbit), soln["time"])

        steps = last_i:step:this_i
        pos_eci = Point3f[soln["state"][1:3,j] for j in steps]
        return pos_eci
    end
    return pos_eci
end

# function plot!(ax::Makie.LScene, mission::Mission)

# end

# Makie.convert_arguments(T::Type{<:Scatter}, pp::Polyhedron) = Makie.convert_arguments(T, pp.vertices[1,:], pp.vertices[2,:],pp.vertices[3,:],)

function plot_state!(fig::Makie.Figure, maneuver::Maneuver, soln::Dict, params)
    torque_p = Makie.Axis(fig[1,1], xlabel="Time [s]", ylabel="Torque [Nm]")
    angrate_p = Makie.Axis(fig[2,1], xlabel="Time [s]", ylabel="Ang. Rate [deg/s]")
    atterr_p = Makie.Axis(fig[3,1], xlabel="Time [s]", ylabel="Attitude Error [deg]")
    ypr_p = Makie.Axis(fig[4,1], xlabel="Time [s]", ylabel="Yaw/Pitch/Roll [deg]")
    axis_p = Makie.Axis(fig[5,1], xlabel="Time [s]", ylabel=L"\mathbf{a}")
    ang_p = Makie.Axis(fig[6,1], xlabel="Time [s]", ylabel=L"\theta \ \mathrm{[deg]}")
    attcheck_p = Makie.Axis(fig[7,1], xlabel="Time [s]", ylabel="DCM Valid")
    # ang2_p = Makie.Axis(fig[1,1])
    # ang3_p = Makie.Axis(fig[1,1])

    torques = backout_free_torque(soln, params)
    angrate = zeros(3, length(soln["time"]))
    dcm_qual = zeros(length(soln["time"]))
    atterr = zeros(3,length(soln["time"]))
    ypr = zeros(3, length(soln["time"]))
    axis = zeros(3, length(soln["time"]))
    ang = zeros(length(soln["time"]))
    t_axis = zeros(3, length(soln["time"]))
    t_ang = zeros(length(soln["time"]))

    for ti in eachindex(soln["time"])
        angrate[:,ti] = soln["state"][ti].angular_velocity
        dcm_qual[ti] = det(soln["state"][ti].attitude) + LinearAlgebra.tr(soln["state"][ti].attitude'*soln["state"][ti].attitude) - 4
        atterr[:,ti] = ang321(maneuver.C) .- ang321(soln["state"][ti].attitude)
        ypr[:,ti] = ang321(soln["state"][ti].attitude)
        axis[:,ti], ang[ti] = axisangle(soln["state"][ti].attitude)
        
        t_M = RotMatrix{3}(soln["state"][ti].attitude)
        t_axis[:,ti] = rotation_axis(t_M)
        t_ang[ti] = rotation_angle(t_M)

    end

    plot_ref = false

    msize = 3
    for i in 1:3
        scatter!(
            torque_p,
            soln["time"] .- soln["time"][1],
            torques[i,:],
            markersize=msize
        )
        scatter!(
            angrate_p,
            soln["time"] .- soln["time"][1],
            angrate[i,:] .* 180/pi,
            markersize=msize
        )
        scatter!(
            atterr_p,
            soln["time"] .- soln["time"][1],
            atterr[i,:] .* 180/pi,
            markersize=msize
        )
        scatter!(
            ypr_p,
            soln["time"] .- soln["time"][1],
            ypr[i,:] .* 180/pi,
            markersize=msize
        )
        scatter!(
            axis_p,
            soln["time"] .- soln["time"][1],
            axis[i,:],
            markersize=msize
        )

        if plot_ref
            lines!(
                axis_p,
                soln["time"] .- soln["time"][1],
                t_axis[i,:],
                linewidth=1
            )
        end
    end
    scatter!(
        ang_p,
        soln["time"] .- soln["time"][1],
        ang .* 180/pi,
        markersize=msize
    )
    if plot_ref
        lines!(
            ang_p,
            soln["time"] .- soln["time"][1],
            t_ang .* 180/pi,
            linewidth=1
        )
    end
    hlines!(
        ypr_p,
        ang321(maneuver.C) .* 180/pi,
        color=:black,
        linestyle=:dash
    )
    lines!(
        attcheck_p,
        soln["time"] .- soln["time"][1],
        dcm_qual
    )

    linkxaxes!(torque_p, angrate_p, atterr_p, ypr_p, axis_p, ang_p, attcheck_p)
end

@recipe(PolyPlot, shape) do scene
    Theme(
        plot_color = :red
    )
end

Makie.args_preferred_axis(::Type{<:PolyPlot}, shape) =  Makie.LScene

function Makie.plot!(polyplot::PolyPlot{<:Tuple{Polyhedron}})
    println(polyplot)
    shape = polyplot[1][]
    scale = min(shape.distances...)
    scatter!(polyplot, 
        shape.vertices[1,:],shape.vertices[2,:],shape.vertices[3,:],
        color=:blue
    )
    
    nscale = 0.01
    for (i, facet) in enumerate(shape.facets)
        # centroid of facet, don't double-count repeated first vertex
        c = mean(facet[:,1:end-1], dims=2)[:]
        # normal to facet, scaled to shape size
        n = cross(facet[:,1] - c)*(facet[:,2] - c)
        n = scale*n/LinearAlgebra.norm(n)

        lines!(polyplot, 
            facet[1,:],facet[2,:],facet[3,:],
            color=:blue
        )

        arrows!(polyplot,
            [Point3f(c)],
            [Vec3f(n)],
            color=:black,
            linewidth=nscale*scale,
            arrowsize=Vec3f(nscale*scale, nscale*scale, 3*nscale*scale)
        )
    end
    polyplot
end