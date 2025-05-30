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
    
    println("saving mission statistics...")
    mission_stats(soln, visibilities, sim)

    # add your own analysis code here...

    # uncomment for an example plot of battery level:
    GLMakie.activate!(title="OpenSOAP")
    fig = Figure(size=(600,400)) # create a new figure
    display(fig) # show the Figure on the screen

    ax = Axis(fig[1,1], xlabel="time [s]", ylabel="battery level [J]") # create an Axis to plot into
    lines!(ax, time, battery, linewidth=0.5)   # lineplot of battery over time
end