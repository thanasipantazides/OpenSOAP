using GLMakie, LinearAlgebra, GeometryBasics
import CairoMakie
import SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using OpenSOAP, Printf

function polar()
    println("loading transformations...")
    eops = SatelliteToolboxTransformations.fetch_iers_eop()
    # compute orbit
    sim = load_mission(joinpath("config", "mission.yaml"))
    
    println("integrating...")
    soln = integrate_system(dynamics_orbit!, sim.initstate, sim.tspan, sim.dt, sim)
    n_t = length(soln["time"])

    println("computing visibilities...")
    target_list = sim.mission.targets
    visibilities = Dict()
    for target in target_list
        visibilities[target.name] = visibility_history(target, soln)
    end

    stats = mission_stats(soln, visibilities, sim)
    time_frac = [stats["state_mean_frac_safe"]; stats["state_mean_frac_power"]; stats["state_mean_frac_downlink"]; stats["state_mean_frac_science"]]
    colors = [:blue; :orange; :lightgreen; :violet]

    GLMakie.activate!()

    set_theme!(theme_dark())
    # t_jd_s = soln["time"][4960]
    t_jd_s = soln["time"][4935]
    # t_jd_s = soln["time"][400]

    (_,i_start) = findmin(x->abs(x - t_jd_s), soln["time"])
    orbit_normal = LinearAlgebra.cross(soln["state"][1:3, i_start], soln["state"][4:6, i_start])
    eye_pos = Vec3d(orbit_normal/norm(orbit_normal) * norm(soln["state"][1:3, i_start]) * 2)

    texture = load_earth_texture_to_ecef(joinpath("assets", "map_diffuse.png"))
    # lighting position of sun:
    pos_eci = SatelliteToolboxCelestialBodies.sun_position_mod(t_jd_s/3600/24)
    sun_light = Vec3f(pos_eci)*1e4

    # set up lighting:
    dl = PointLight(-Point3f(sun_light), RGBf(243/255, 241/255, 218/255))
    el = PointLight(10*Point3f(eye_pos), RGBf(1,1,1))
    # dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), sun_light)
    # al = AmbientLight(RGBf(0.1, 0.1, 0.1))
    # al = AmbientLight(RGBf(0.2, 0.2, 0.2))
    al = AmbientLight(RGBf(0.1, 0.1, 0.1))

    fig = Figure(size=(1100,800))
    ax = LScene(
        fig[1:2,1:2], 
        show_axis=false, 
        # height=400,
        # width=400,
        scenekw=(
            lights=[dl, al, el], 
            backgroundcolor=:black, 
            clear=true
        )
    )
    pie_ax = Axis(fig[1,3], backgroundcolor=:black, aspect=AxisAspect(1), title="Mean fraction in mode", height=200, width=200)
    # pie_ax = Axis(fig[1,3], backgroundcolor=:black, axis=(autolimitaspect = 1, ), title="Fraction in mode", height=200, width=200)
    hidedecorations!(pie_ax)

    labels = [ @sprintf "%0.2f" time for time in time_frac]

    r1 = 6
    for (i, time) in enumerate(time_frac)
        δθ = 2*pi*time
        θ = 2*pi*sum(time_frac[1:i]) - δθ/2
        text!(
            pie_ax, 
            labels[i], 
            position=r1.*(cos(θ), sin(θ)), 
            rotation=0, 
            fontsize=16, 
            space=:data, 
            align=(:center, :center)
        )
    end

    hud_conop!(ax, t_jd_s, visibilities, target_list, soln, texture, eops)
    
    # fig[1,3] = Legend(fig, ax, "Mode", framevisible = false, labelsize=10.0, patchsize=(20.0f0, 10.0f0), halign=:left)
    # manually build legend: 
    safe_patch = [PolyElement(color = fill((:blue, 1.0), 20, 2), strokecolor = :black, strokewidth = 1)]
    power_patch = [PolyElement(color = fill((:orange, 1.0), 20, 2), strokecolor = :black, strokewidth = 1)]
    downlink_patch = [PolyElement(color = fill((:lightgreen, 1.0), 20, 2), strokecolor = :black, strokewidth = 1)]
    science_patch = [PolyElement(color = fill((:violet, 1.0), 20, 2), strokecolor = :black, strokewidth = 1)]
    
    Legend(
        fig[2,3],
        [safe_patch, power_patch, downlink_patch, science_patch],
        ["safe", "power", "downlink", "science"],
        "Modes",
        patchsize=(35,35), rowgap=10
    )

    pie!(pie_ax, time_frac, color=colors, radius=4, inner_radius=2)
    pie_scale = 8
    ylims!(pie_ax, -pie_scale, pie_scale)
    xlims!(pie_ax, -pie_scale, pie_scale)
    
    save(joinpath("cases","hud_conop.png"), fig, update=false, px_per_unit=10.0)

    # display(fig, update=false)
    
    # focus_button = Button(fig[3,1], label="look")
    # on(focus_button.clicks, priority=1) do n
    #     cam3d!(
    #         ax, 
    #         projectiontype="Orthographic", 
    #         lookat=Vec3d(0), 
    #         upvector=Vec3d(0,0,1), 
    #         eyeposition=eye_pos
    #     )
    # end

    
end

function main()
    println("loading transformations...")
    eops = SatelliteToolboxTransformations.fetch_iers_eop()
    # compute orbit
    sim = load_mission(joinpath("config","mission.yaml"))
    
    println("integrating...")
    soln = integrate_system(dynamics_orbit!, sim.initstate, sim.tspan, sim.dt, sim)
    n_t = length(soln["time"])

    println("computing visibilities...")
    target_list = sim.mission.targets
    visibilities = Dict()
    for target in target_list
        visibilities[target.name] = visibility_history(target, soln)
    end
    
    println("saving mission statistics...")
    mission_stats(soln, visibilities, sim)

    # plotting
    println("plotting...")
    GLMakie.activate!(title="OpenSOAP")
    
    fig = Figure(size=(1400,840))
    display(fig)
    
    # make observable slider:
    time_sliders = SliderGrid(fig[6, 2], (label = "time", format = "{:d}", range = 1:1:n_t, startvalue = 500))
    slider_observables = [s.value for s in time_sliders.sliders]
    
    t_jd_s = lift(slider_observables[1]) do t_i
        return soln["time"][t_i]
    end

    play_button = Button(fig[6,5], label="play")
    back_button = Button(fig[6,1], label="◀")
    fwrd_button = Button(fig[6,3], label="▶")

    # lighting position of sun:
    sun_light = lift(t_jd_s) do t_jd_s
        pos_eci = SatelliteToolboxCelestialBodies.sun_position_mod(t_jd_s/3600/24)
        return Vec3f(pos_eci)
    end

    # set up lighting:
    # pl = PointLight(Point3f([4*6371e3;0;0]), RGBf(20, 20, 20))
    # dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), Vec3f(-1, 0, 0))
    dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), sun_light)
    al = AmbientLight(RGBf(0.1, 0.1, 0.1))

    # start main scene:
    ax = LScene(
        fig[1:5,1:3], 
        show_axis=false, 
        scenekw=(lights=[dl, al], 
        backgroundcolor=:black, 
        clear=true)
    )
    # populate auxiliary axes:
    detail_ax = Axis(
        fig[5,4], 
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], -0.2, 0.2), 
        title="Angular rates in body frame", 
        xlabel="Time [s]", 
        ylabel="Angular rate [rad/s]"
    )
    visible_ax = Axis(
        fig[2,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, length(visibilities)), 
        xlabel="Time [s]", 
        ylabel="Visible?"
    )
    mode_ax = Axis(
        fig[1,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, 1),
        xlabel="Time [s]", 
        ylabel="Mode"
    )
    power_ax = Axis(
        fig[3,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, sim.mission.spacecraft.power.capacity/3600), 
        xlabel="Time [s]", 
        ylabel="Battery capacity [Wh]"
    )
    data_ax = Axis(
        fig[4,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, sim.mission.spacecraft.data.capacity/8e6),
        xlabel="Time [s]", 
        ylabel="Data storage [MB]"
    )    
    
    # load Earth texture (these are all equirectangular projection/plate carreé):
    texture = load_earth_texture_to_ecef(joinpath("assets","map_diffuse.png"))
    # texture = load_earth_texture_to_ecef("assets/map_bathy.png")
    # texture = load_earth_texture_to_ecef("assets/map_veggie.jpeg")
    # texture = load_earth_texture_to_ecef("assets/map_pol.png")
    
    set_theme!(theme_dark())
    
    plot_earth!(ax, t_jd_s, eops, texture)
    plot_meta!(ax, t_jd_s, true, visibilities, soln, eops)
    # plot_frames!(ax, t_jd_s, visibilities, soln, eops)
    plot_spacecraft!(ax, t_jd_s, 10000, soln)
    # plot_targets!(ax, target_list, t_jd_s, soln, eops)
    plot_detail!(detail_ax, t_jd_s, soln)
    plot_visibilities!(visible_ax, t_jd_s, visibilities, soln)
    plot_power!(power_ax, t_jd_s, visibilities, soln)
    plot_data!(data_ax, t_jd_s, visibilities, soln)
    plot_mode!(mode_ax, t_jd_s, visibilities, soln)

    lookat_orbit!(ax, t_jd_s, soln)

    fig[2:5,5] = Legend(fig, visible_ax, "Visibility", framevisible = false, labelsize=10.0, patchsize=(20.0f0, 10.0f0), valign=:top, halign=:left)
    fig[1,5] = Legend(fig, mode_ax, "Mode", framevisible = false, labelsize=10.0, patchsize=(20.0f0, 10.0f0), halign=:left)
    linkxaxes!(detail_ax, visible_ax, power_ax, data_ax, mode_ax)
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