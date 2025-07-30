using OpenSOAP, LinearAlgebra
import SatelliteToolboxTransformations
using GLMakie

using Profile, PProf

@doc raw"""
    main()

This is a boilerplate example for running an orbit simulation, checking target visibility, and writing outputs.
"""
function main()

    Profile.clear()
    println("loading transformations...")
    eops = SatelliteToolboxTransformations.fetch_iers_eop()

    println("loading parameters...")
    sim = load_mission(joinpath("config", "mission.yaml"))
    try # create the "cases" folder to store results, if you don't have it already
        mkdir("cases")
    catch e
    end

    println("trying new code...")
    @profile times, states, target_visibilities, target_choice = simulate(sim, Dict{String,String}())

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
    position_I = soln["state"][1:3, :]
    velocity_I = soln["state"][4:6, :]
    angular_velocity_B = soln["state"][7:9, :]
    attitude_BI = [reshape(soln["state"][10:18, k], 3, 3) for k in 1:length(time)]
    battery = soln["state"][19, :]
    storage = soln["state"][20, :]
    # mode = soln["state"][21, :]

    x_I = [state.position[1] for state in states]
    y_I = [state.position[2] for state in states]
    z_I = [state.position[3] for state in states]
    mode = [state.mode for state in states]

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
    fig = Figure(size=(1400, 900)) # create a new figure
    display(fig) # show the Figure on the screen

    plot_moc!(fig, times, states, target_visibilities, target_choice, zeros(3,3,1), sim)

    # out = [(target_visibilities[t,:], target) for (t, target) in enumerate(sim.mission.targets) if isa(target, SunTarget)]
    # sun_vis = out[1][1]
    # sun = out[1][2]

    # sun_I = position_eci(sun, times[1])
    # sun_I = 6371e3 * 1.25 * sun_I / norm(sun_I)

    # times_d = (times .- times[1])/3600/24
    # target_idx = [t for t in eachindex(sim.mission.targets)]
    # target_labels = [target.name for target in sim.mission.targets]

    # ax1 = Axis(fig[1, 1], xlabel="time [days]", ylabel="Mode") # create an Axis to plot into
    # lines!(ax1, times_d, mode, linewidth=0.5)   # lineplot of battery over time
    # ax2 = Axis(fig[2,1], xlabel="time [days]", ylabel="target visible", yticks = (target_idx, target_labels))

    # # heatmap!(ax2, Resampler(target_visibilities'), colormap=cgrad(["white","black"],2))
    # plot_visibilities!(ax2, times_d, states, target_visibilities, target_choice, sim)
    # linkxaxes!(ax1,ax2)
    # ax3 = LScene(fig[3,1])
    # lines!(ax3, x_I, y_I, z_I)
    # scatter!(ax3, sun_I, color=:yellow)

    # pprof()
end
