using OpenSOAP, LinearAlgebra
import SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using GLMakie

@doc raw"""
    groundstations()

This is a boilerplate example for running an orbit simulation, checking target visibility, and displaying visibilities for each target.
"""
function groundstations()
    println("loading transformations...")
    eops = SatelliteToolboxTransformations.fetch_iers_eop()

    println("loading parameters")
    sim = load_mission(joinpath("config","mission.yaml"))

    try # create the "cases" folder to store results, if you don't have it already
        mkdir("cases")
    catch e
    end

    println("integrating...")
    soln = integrate_system(dynamics_orbit!, sim.initstate, sim.tspan, sim.dt, sim)
    time = soln["time"]

    # soln contains the time history of the simulation and the state vector at each point in time. 
    # - soln["time"] is a vector of each time step, stored as seconds but scaled to Julia Day.
    #   - to recover Julian Day from soln["time"], do soln["time"]/24/3600
    # - soln["state"] is a matrix. Each column is the state vector for the corresponding timestep in soln["time"]. Each row is:
    #   - soln["state"][1:3] is the spacecraft position, expressed in the Earth-Centered, Inertial frame.
    #   - soln["state"][4:6] is the spacecraft velocity, expressed in the Earth-Centered, Inertial frame.
    #   - soln["state"][7:9] is the spacecraft angular velocity of the Body frame relative to the Earth-Centered, Inertial frame.
    #   - soln["state"][10:18] is the spacecraft attitude of the Body frame relative to the Earth-Centered, Inertial frame.
    #   - soln["state"][19] is the battery level.
    #   - soln["state"][20] is the disk usage.
    #   - soln["state"][21] is a number corresponding to mode.

    # Below, I separate each state element into its own variable. 
    # The suffix I indicates Earth-Centered, Inertial frame. The suffix B indicates Body frame.
    position_I = soln["state"][1:3,:]
    velocity_I = soln["state"][4:6,:]
    angular_velocity_B = soln["state"][7:9,:]
    attitude_BI = [reshape(soln["state"][10:18, k], 3,3) for k in 1:length(time)]
    battery = soln["state"][19, :]
    storage = soln["state"][20, :]
    mode = soln["state"][21, :]

    # compute when each target is visible from spacecraft:
    println("computing visibilities...")
    target_list = sim.mission.targets
    visibilities = Dict()
    for target in target_list
        visibilities[target.name] = visibility_history(target, soln)
    end

    # by default, the Sun is included in the target list. Filter it out here so we just worry about groundstations:
    ground_targets = [target for target in target_list if typeof(target) === GroundTarget]
    
    println("saving mission statistics...")
    mission_stats(soln, visibilities, sim)

    # make a "day of mission" time variable for clearer plotting:
    t_day = (time .- time[1]) / 3600 / 24



    # add your own analysis code here...


    
    # an example of plotting visibility for each target:
    #   Note: the Sun is included in the targets list here. 
    #   Also the target list is very large with all groundstations included, which makes this plot very hard to read without zooming.
    GLMakie.activate!(title="OpenSOAP")
    fig = Figure(size=(600,600)) # create a new figure
    display(fig) # show the Figure on the screen

    # make a new axis inside the figure to plot visibility of targets:
    vis_ax = Axis(fig[1,1], xlabel="mission time [days]", ylabel="target eclipse") # create an Axis to plot into

    # plot the visibility (1 if visible, 0 if blocked) for each target:
    for target in ground_targets
        lines!(vis_ax, t_day, visibilities[target.name], linewidth=0.5, label=target.name)   # lineplot of battery over time
    end
    fig[1:3,2] = Legend(fig, vis_ax, "Target", framevisible = false, labelsize=10.0, patchsize=(20.0f0, 10.0f0), halign=:left)
    
    # 3D plot of Earth and groundstations
    # lighting for the plot:
    al = AmbientLight(RGBf(243/255, 241/255, 230/255))
    ax_globe = LScene(
        fig[2:3,1], 
        show_axis=false, 
        scenekw=(
            # lights=[dl, al], 
            lights = [al],
            # backgroundcolor=:black, 
            clear=true
        )
    )
    globe_plot_time = time[1]
    do_gs_shading = true
    texture = load_earth_texture_to_ecef(joinpath("assets","map_bw.png"))
    plot_earth_static!(ax_globe, globe_plot_time, eops, texture)
    plot_targets_static!(ax_globe, ground_targets, do_gs_shading, globe_plot_time, soln, eops)

end