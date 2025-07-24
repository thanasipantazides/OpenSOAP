using OpenSOAP, LinearAlgebra  
import SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using GLMakie 
using Statistics  # For calculating mean

#====================================================================================
This function simulates the orbit of a satellite, calculates how often it's visible 
to specific ground stations, and plots the results! :)
=====================================================================================#

function groundstations()
    println("loading transformations...")
    # Load Earth orientation parameters (needed for accurate coordinate transformations)
    eops = SatelliteToolboxTransformations.fetch_iers_eop()

    println("loading parameters")
    # Load the mission configuration from the YAML file
    sim = load_mission(joinpath("config","mission.yaml"))

    # Try to make a folder to store output, ignore error if it already exists
    try
        mkdir("cases")
    catch e
    end

    # Simulate the satellite's orbit and state over time
    soln = integrate_system(dynamics_orbit!, sim.initstate, sim.tspan, sim.dt, sim)
    time = soln["time"]  # Time steps of the simulation

    # Extract different parts of the satellite's state from the simulation
    position_I = soln["state"][1:3,:]  # Position in inertial frame
    velocity_I = soln["state"][4:6,:]  # Velocity in inertial frame
    angular_velocity_B = soln["state"][7:9,:]  # Angular velocity in body frame
    attitude_BI = [reshape(soln["state"][10:18, k], 3,3) for k in 1:length(time)]  # Attitude matrix at each timestep
    battery = soln["state"][19, :]  # Battery level
    storage = soln["state"][20, :]  # Data storage level
    mode = soln["state"][21, :]  # Operational mode

    println("computing visibilities...")
    # Goes through the list of all targets defined in the mission
    target_list = sim.mission.targets
    visibilities = Dict()  #dictionary to store visibility data for each target

    # Calculates when the satellite can see target
    for target in target_list
        visibilities[target.name] = visibility_history(target, soln)
    end

    # This only keeps ground targets
    ground_targets = [target for target in target_list if typeof(target) === GroundTarget]

    # Converts time to days for plotting!
    t_day = (time .- time[1]) / 3600 / 24
    Δt = sim.dt  # Time step size in seconds

    # Calculates the total simulation time and number of orbits completed
    total_mission_time_sec = time[end] - time[1]
    orbit_period_sec = 90 * 60  # Assuming a 90-minute orbit
    num_orbits = total_mission_time_sec / orbit_period_sec

    # Print out a summary for each ground target
    println("\n==== Normalized Ground Target Visibility Summary ====")

    for target in ground_targets
        vis = visibilities[target.name]

        # Analyzes when the satellite was in view of the groundstation
        pass_durations = Float64[]  # Store durations of each pass
        in_pass = false
        start_idx = 0

        #loop of a visibility array: (1 = visible, 0 = not visible)
        for (i, val) in enumerate(vis)
            if val == 1 && !in_pass
                # Start of a new pass
                start_idx = i
                in_pass = true
            elseif val == 0 && in_pass
                # End of a pass = store duration!
                push!(pass_durations, (i - start_idx) * Δt)
                in_pass = false
            end
        end
        if in_pass
            # If still in pass at end of data, close it!
            push!(pass_durations, (length(vis) - start_idx + 1) * Δt)
        end

        # This converts pass durations from seconds to minutes
        durations_min = pass_durations ./ 60
        total_visible_time_hr = sum(durations_min) / 60  # Converts total visible time to hours
        passes = length(durations_min)

        # Prints out the stats for each ground target
        println("Target: $(target.name)")
        println("  Total visible time: $(round(total_visible_time_hr, digits=2)) hours")
        println("  Number of passes: $(passes)")
        if passes > 0
            println("  Min pass duration: $(round(minimum(durations_min), digits=2)) min")
            println("  Mean pass duration: $(round(mean(durations_min), digits=2)) min")
            println("  Max pass duration: $(round(maximum(durations_min), digits=2)) min")
            println("  Mean visible time per orbit: $(round(total_visible_time_hr / num_orbits, digits=2)) hr/orbit")
            println("  Mean passes per orbit: $(round(passes / num_orbits, digits=2)) passes/orbit\n")
        else
            println("  No passes detected.\n")
        end
    end

    #======== Visualization Section =============
    =#

    # Makes a 2D figure window
    GLMakie.activate!(title="OpenSOAP")
    fig = Figure(size=(600,600))
    display(fig)

    # Visibility timeline for each target
    vis_ax = Axis(fig[1,1], xlabel="mission time [days]", ylabel="target eclipse")
    for target in ground_targets
        lines!(vis_ax, t_day, visibilities[target.name], linewidth=0.5, label=target.name)
    end

    #This adds a legend to explain which line represents which target
    fig[1:3,2] = Legend(fig, vis_ax, "Target", framevisible = false, labelsize=10.0, patchsize=(20.0f0, 10.0f0), halign=:left)

    # This is the 3D Earth plot with visible groundstation locations
    al = AmbientLight(RGBf(243/255, 241/255, 230/255))  # Makes it look fancy! :)
    ax_globe = LScene(
        fig[2:3,1], 
        show_axis=false, 
        scenekw=(lights = [al], clear=true)
    )

    globe_plot_time = time[1]  # This picks a single time to visualize groundstation locations
    do_gs_shading = true  # Whether to shade the groundstations
    texture = load_earth_texture_to_ecef(joinpath("assets","map_bw.png"))  # This loads Earth's map
    plot_earth_static!(ax_globe, globe_plot_time, eops, texture)  # Plots Earths surface
    plot_targets_static!(ax_globe, ground_targets, do_gs_shading, globe_plot_time, soln, eops)  # This plots station locations

end
