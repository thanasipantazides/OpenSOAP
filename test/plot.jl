using GLMakie, LinearAlgebra
import CairoMakie
import SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using OpenSOAP
using Filesystem

@doc raw"""
    setup_parameters()

Convenience function to build all necessary simulation data, and create a `LEOSimulation` structure to pass to the integrator.
"""
function setup_parameters()::LEOSimulation

    leo_sim = load_mission(joinpath("config","mission.yaml"))
    try 
        mkdir("cases")
    catch e
    end

    return leo_sim
end

@doc raw"""
    run_orbit(sim::LEOSimulation)

Perform numerical integration of simulation problem defined in `LEOSimulation` structure, and return `soln::Dict` with solution data (keys `"time"`, a time history; and `"state"`, a state vector history).
"""
function run_orbit(sim::LEOSimulation)
    soln = integrate_system(dynamics_orbit!, sim.initstate, sim.tspan, sim.dt, sim)
    return soln
end

function save_power(soln)
    CairoMakie.activate!(type="png")
    set_theme!(theme_light())
    fig = Figure()
    ax1 = Axis(
        fig[1,1],
        title = "Battery capacity history",
        xlabel = "Time [days]",
        ylabel = "Battery capacity [Wh]"
    )
    ax2 = Axis(
        fig[2,1],
        title = "Mode history",
        xlabel = "Time [days]",
        ylabel = "Mode"
    )
    lines!(ax1, (soln["time"] .- soln["time"][1])/24/3600, soln["state"][19,:]/3600)
    lines!(ax2, (soln["time"] .- soln["time"][1])/24/3600, soln["state"][21,:])
    save(joinpath("cases","battery_dod.pdf"), fig)
end

function plot_main()
    println("loading transformations...")
    eops = SatelliteToolboxTransformations.fetch_iers_eop()
    # compute orbit
    sim = setup_parameters()
    println("integrating...")
    soln = run_orbit(sim)
    n_t = length(soln["time"])

    # compute when each target is visible from spacecraft:
    println("computing visibilities...")
    target_list = sim.mission.targets
    visibilities = Dict()
    for target in target_list
        visibilities[target.name] = visibility_history(target, soln)
    end
    
    println("saving mission statistics...")
    mission_stats(soln, visibilities, sim)

    save_power(soln)

    # plotting
    println("plotting...")
    GLMakie.activate!(title="OpenSOAP")
    
    fig = Figure(size=(1400,840))
    display(fig)
    
    # make observable slider:
    time_sliders = SliderGrid(fig[6, 2], (label = "time", format = "{:d}", range = 1:1:n_t, startvalue = 1))
    slider_observables = [s.value for s in time_sliders.sliders]
    
    t_jd_s = lift(slider_observables[1]) do t_i
        return soln["time"][t_i]
    end

    play_button = Button(fig[6,5], label="play")
    back_button = Button(fig[6,1], label="<<")
    fwrd_button = Button(fig[6,3], label=">>")

    # lighting position of sun:
    sun_light = lift(t_jd_s) do t_jd_s
        pos_eci = SatelliteToolboxCelestialBodies.sun_position_mod(t_jd_s/3600/24)
        return Vec3f(pos_eci)
    end

    

    # set up lighting:
    # pl = PointLight(Point3f([4*6371e3;0;0]), RGBf(20, 20, 20))
    # dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), Vec3f(-1, 0, 0))
    dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), sun_light)
    # dl = PointLight(RGBf(243/255, 241/255, 218/255), sun_light)
    al = AmbientLight(RGBf(0.3, 0.3, 0.3))

    # start main scene:
    ax = LScene(
        fig[1:5,1:3], 
        show_axis=false, 
        scenekw=(lights=[dl, al], 
        backgroundcolor=:black, 
        clear=true)
    )
    # populate auxiliary axes:
    # detail_ax = Axis(
    #     fig[5,4], 
    #     backgroundcolor=:black, 
    #     limits=(0, soln["time"][end] - soln["time"][1], -0.2, 0.2), 
    #     title="Angular rates in body frame", 
    #     xlabel="Time [s]", 
    #     ylabel="Angular rate [rad/s]"
    # )
    visible_ax = Axis(
        fig[2,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, length(visibilities)), 
        # title="Target visibility", 
        xlabel="Time [s]", 
        ylabel="Visible?"
    )
    mode_ax = Axis(
        fig[1,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, 1),
        # title="Mode", 
        xlabel="Time [s]", 
        ylabel="Mode"
    )
    power_ax = Axis(
        fig[3,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, sim.mission.spacecraft.power.capacity/3600), 
        # title="Power", 
        xlabel="Time [s]", 
        ylabel="Battery capacity [Wh]"
    )
    data_ax = Axis(
        fig[4,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, sim.mission.spacecraft.data.capacity/8e6),
        # title="Data", 
        xlabel="Time [s]", 
        ylabel="Data storage [MB]"
    )    
    
    # load Earth texture (these are all equirectangular projection/plate carreé):
    texture = load_earth_texture_to_ecef(joinpath("assets","map_diffuse.png"))
    # texture = load_earth_texture_to_ecef(joinpath("assets","map_bathy.png"))
    # load_earth_texture_to_ecef(joinpath("assets","map_veggie.png"))
    
    set_theme!(theme_dark())
    
    plot_earth!(ax, t_jd_s, eops, texture)
    plot_spacecraft!(ax, t_jd_s, 10000, soln)
    plot_targets!(ax, target_list, t_jd_s, soln, eops)
    plot_frames!(ax, t_jd_s, soln, eops)
    plot_magnetic_field!(ax, t_jd_s, soln, eops)
    # plot_detail!(detail_ax, t_jd_s, soln)
    plot_visibilities!(visible_ax, t_jd_s, visibilities, soln)
    plot_power!(power_ax, t_jd_s, visibilities, soln)
    plot_data!(data_ax, t_jd_s, visibilities, soln)
    plot_mode!(mode_ax, t_jd_s, visibilities, soln)

    fig[2:5,5] = Legend(fig, visible_ax, "Visibility", framevisible = false, labelsize=10.0, patchsize=(20.0f0, 10.0f0), valign=:top, halign=:left)
    fig[1,5] = Legend(fig, mode_ax, "Mode", framevisible = false, labelsize=10.0, patchsize=(20.0f0, 10.0f0), halign=:left)
    linkxaxes!(visible_ax, power_ax, data_ax, mode_ax)
    # linkxaxes!(detail_ax, visible_ax, power_ax, data_ax, mode_ax)
    hidexdecorations!(visible_ax, grid=false, ticks=false)
    hideydecorations!(visible_ax, label=false)
    hidexdecorations!(mode_ax, grid=false, ticks=false)
    hideydecorations!(mode_ax, label=false)
    hidexdecorations!(power_ax, grid=false, ticks=false)

    # callback for back arrow
    on(back_button.clicks, priority=2) do n
        slider_pos = slider_observables[1][] - 1
        Makie.set_close_to!(time_sliders.sliders[1], slider_pos)
    end

    # callback for forward arrow
    on(fwrd_button.clicks, priority=2) do n
        slider_pos = slider_observables[1][] + 1
        Makie.set_close_to!(time_sliders.sliders[1], slider_pos)
    end
    
    # callback for play button
    on(play_button.clicks, priority=1) do n
        if play_button.label == "play"
            
        else
            play_button.label = "stop"
            frame_rate = 60
            step = Int64(round(n_t / 100)) + 1
            @async for i in 1:step:n_t
                Makie.set_close_to!(time_sliders.sliders[1], i)
                sleep(1/frame_rate)
            end
            play_button.label = "play"
        end
        Consume(true)
    end
end