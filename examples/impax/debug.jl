using OpenSOAP, LinearAlgebra
import SatelliteToolboxTransformations
using GLMakie

@doc raw"""
    main()

This is a boilerplate example for running an orbit simulation, checking target visibility, and writing outputs.
"""
function main()

    println("loading transformations...")
    eops = SatelliteToolboxTransformations.fetch_iers_eop()

    println("loading parameters...")
    sim = load_mission(joinpath("config", "mission.yaml"))
    try # create the "cases" folder to store results, if you don't have it already
        mkdir("cases")
    catch e
    end

    times, states, target_visibilities, target_choice, reference_directions, reference_attitude = simulate(sim, Dict{String,String}())

    x_I = [state.position[1] for state in states]
    y_I = [state.position[2] for state in states]
    z_I = [state.position[3] for state in states]
    mode = [state.mode for state in states]

    GLMakie.activate!(title="OpenSOAP")
    fig = Figure(size=(900, 900)) # create a new figure
    display(fig) # show the Figure on the screen

    # plotting
    axes_colors = [:red, :green, :blue]
    al = AmbientLight(RGBf(0.4,0.4,0.4))
    ax = LScene(
        fig[1:3,1], 
        show_axis = false, 
        scenekw = (
            # lights=[dl, al], 
            lights = [al],
            # backgroundcolor=:black, 
            clear=true
        ),
        tellheight = false
    )
    
    stepw = 100
    r_I = Point3f[states[k].position for k in 1:stepw:length(times)]
    
    x_B = Vec3f[states[k].attitude[:,1] for k in 1:stepw:length(times)]
    y_B = Vec3f[states[k].attitude[:,2] for k in 1:stepw:length(times)]
    z_B = Vec3f[states[k].attitude[:,3] for k in 1:stepw:length(times)]
    
    ref_I = Vec3f[reference_directions[:,k] / norm(reference_directions[:,k]) for k in 1:stepw:length(times)]
    
    ecef_arrow_scale = 0.005
    lengthscale = 0.1
    r_E = 6371e3
    
    lines!(ax, r_I, color=:black)
    arrows!(ax, r_I, x_B, 
        color=:red, 
        linewidth=ecef_arrow_scale*r_E,
        arrowsize=Vec3f(ecef_arrow_scale*r_E, ecef_arrow_scale*r_E, ecef_arrow_scale*3*r_E),
        lengthscale=r_E*lengthscale
    )
    arrows!(ax, r_I, y_B, 
        color=:green, 
        linewidth=ecef_arrow_scale*r_E,
        arrowsize=Vec3f(ecef_arrow_scale*r_E, ecef_arrow_scale*r_E, ecef_arrow_scale*3*r_E),
        lengthscale=r_E*lengthscale
    )
    arrows!(ax, r_I, z_B, 
        color=:blue, 
        linewidth=ecef_arrow_scale*r_E,
        arrowsize=Vec3f(ecef_arrow_scale*r_E, ecef_arrow_scale*r_E, ecef_arrow_scale*3*r_E),
        lengthscale=r_E*lengthscale
    )
    arrows!(ax, r_I, ref_I, 
        color=:yellow,
        linewidth=ecef_arrow_scale*r_E*0.5,
        arrowsize=Vec3f(ecef_arrow_scale*r_E, ecef_arrow_scale*r_E, ecef_arrow_scale*3*r_E),
        lengthscale=r_E*lengthscale*2
    )
    scatter!(ax, [r_I[1]], color=:black, markersize=15)
    
    return fig
end

function check_interp()
    C_BI1 = diagm([1.0;1;1])
    # C_BI2 = r_euler3(pi/4)
    C_BI2 = r_random()
    n = 5
    seq = rotinterp(C_BI1, C_BI2, n)
    
    GLMakie.activate!(title="OpenSOAP")
    fig = Figure(size=(900, 900)) # create a new figure
    display(fig) # show the Figure on the screen

    # plotting
    axes_colors = [:red, :green, :blue]
    al = AmbientLight(RGBf(0.4,0.4,0.4))
    ax = LScene(
        fig[1:3,1], 
        show_axis = true, 
        scenekw = (
            # lights=[dl, al], 
            lights = [al],
            # backgroundcolor=:black, 
            clear=true
        ),
        tellheight = false
    ) 
    
    arrows!(ax, [Point3f(0.0)], [Vec3f(C_BI1[1,1], C_BI1[2,1], C_BI1[3,1])], color=:red)
    arrows!(ax, [Point3f(0.0)], [Vec3f(C_BI1[1,2], C_BI1[2,2], C_BI1[3,2])], color=:green)
    arrows!(ax, [Point3f(0.0)], [Vec3f(C_BI1[1,3], C_BI1[2,3], C_BI1[3,3])], color=:blue)
    arrows!(ax, [Point3f(0.0)], [Vec3f(C_BI2[1,1], C_BI2[2,1], C_BI2[3,1])], color=:red)
    arrows!(ax, [Point3f(0.0)], [Vec3f(C_BI2[1,2], C_BI2[2,2], C_BI2[3,2])], color=:green)
    arrows!(ax, [Point3f(0.0)], [Vec3f(C_BI2[1,3], C_BI2[2,3], C_BI2[3,3])], color=:blue)
    
    for k in 1:n
        arrows!(ax, [Point3f(0.0)], [Vec3f(seq[1,1,k], seq[2,1,k], seq[3,1,k])], color=:magenta)
        arrows!(ax, [Point3f(0.0)], [Vec3f(seq[1,2,k], seq[2,2,k], seq[3,2,k])], color=:yellow)
        arrows!(ax, [Point3f(0.0)], [Vec3f(seq[1,3,k], seq[2,3,k], seq[3,3,k])], color=:cyan)
    end
end