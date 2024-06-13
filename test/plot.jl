using GLMakie, LinearAlgebra
import SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using openSOAR

function make_solar_panels()::Vector{SolarPanel}
    pxf = SolarPanel([1.0;0;0], 0.295, 0.005*4)
    pxb = SolarPanel([-1.0;0;0], 0.295, 0.005*4)
    pyf = SolarPanel([0.0;1;0], 0.295, 0.005*6)
    pyb = SolarPanel([0.0;-1;0], 0.295, 0.005*6)
    pzf = SolarPanel([0.0;0;1], 0.295, 0.005*12)
    pzb = SolarPanel([0.0;0;-1], 0.295, 0.005*12)
    return [pxf; pxb; pyf; pyb; pzf; pzb]
end

function setup_parameters()::LEOSimulation
    
    earth_data = EarthProperties(3.986e14, 1.081874e-4, SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS, 1361)

    start_time_jd = SatelliteToolboxTransformations.date_to_jd(2027, 11, 28, 0, 25, 0)
    start_time_s = start_time_jd * 3600 * 24
    duration_s = 3600*100
    tspan = [start_time_s, start_time_s + duration_s]
    dt_s = 1

    inc = 70*pi/180
    r0 = [earth_data.r+550e3; 0; 0]
    v0m = sqrt(earth_data.mu/norm(r0))
    v0 = 1.0*[0; v0m*cos(inc); v0m*sin(inc)]
    w0 = [1;1;1]*1e-2
    C_BI0 = diagm([1;1;1])
    E0 = 40*3600

    x0 = [r0;v0;w0;vec(C_BI0);E0]

    mass_data = MassProperties(
        10.0, # kg
        diagm([2;1;4])*1e-3
    )
    power_data = PowerProperties(
        84*60*60,   # Whr to J
        10,         # W
        make_solar_panels()
    )
    spacecraft_data = SpacecraftProperties(
        "impax",
        power_data,
        # data_data,
        mass_data
    )
    mission_data = Mission(
        "impax",
        spacecraft_data
    )

    leo_sim = LEOSimulation(
        earth=earth_data,
        mission=mission_data,
        tspan=tspan,
        dt=dt_s,
        initstate=x0
    )

    return leo_sim
end

function run_orbit(tspan::Vector{<:Real}, dt::Real)

    m = 10 # kg
    I = diagm([2;1;4])*1e-3
    p_system = SystemData(m,I)
    p_earth = EarthData(3.986e14, 1.081874e-4, 6371e3, 1361)
    p_sim = SimulationParameters(p_earth, p_system)
    
    inc = 70*pi/180
    x0 = [p_sim.earth.r+550e3; 0; 0]
    v0m = sqrt(p_sim.earth.mu/norm(x0))
    v0 = 1.0*[0; v0m*cos(inc); v0m*sin(inc)]
    w0 = [1;1;1]*1e-2
    C_BI0 = I

    x0 = [x0;v0;w0;vec(C_BI0)]
    
    soln = integrate_system(dynamics_orbit!, x0, tspan, dt, p_sim)
    return soln
end

function run_orbit(sim::LEOSimulation)
    soln = integrate_system(dynamics_orbit!, sim.initstate, sim.tspan, sim.dt, sim)
    return soln
end

function gen_targets(targets::Vector{FrameFixedTarget}, times_s::Vector{<:Real}, start_jd::Real)

end

function plot_main()
    # compute orbit
    println("integrating...")
    sim = setup_parameters()
    soln = run_orbit(sim)
    n_t = length(soln["time"])

    # compute targets
    println("making targets...")
    eops = SatelliteToolboxTransformations.fetch_iers_eop()

    sun = FrameFixedTarget(
        "sun",
        3,
        [0;0;0],
        [0;0;0],
        pi,
        :ICRS,
        eops
    )
    
    gmt_loc = SatelliteToolboxTransformations.geodetic_to_ecef(51.4934*π/180, 0.0, 3)
    mpls_loc = SatelliteToolboxTransformations.geodetic_to_ecef(44.9778*π/180, -93.265*π/180, 3)
    gmt_loc = Vector(gmt_loc[:])
    mpls_loc = Vector(mpls_loc[:])
    angt = -pi/2
    C_target = [
        cos(angt) -sin(angt) 0;
        sin(angt) cos(angt) 0;
        0 0 1
    ]
    targecef1 = FrameFixedTarget(
        "station1",
        1,
        gmt_loc,
        [1;0;0],
        pi/3,
        :ECEF,
        eops
    )
    
    targetecef2 = FrameFixedTarget(
        "station2", 
        2, 
        mpls_loc,
        [1;0;0],
        pi/6,
        :ECEF,
        eops
    )

    target_list = [sun; targecef1; targetecef2]

    # plotting
    println("plotting...")
    GLMakie.activate!()
    
    fig = Figure(size=(1400,840))
    display(fig)
    
    # make observable slider:
    time_sliders = SliderGrid(fig[5, 1:3], (label = "time", format = "{:d}", range = 1:1:n_t, startvalue = 1))
    slider_observables = [s.value for s in time_sliders.sliders]
    
    t_jd_s = lift(slider_observables[1]) do t_i
        return soln["time"][t_i]
    end

    play_button = Button(fig[5,5], label="play")

    # lighting position of sun:
    sun_light = lift(t_jd_s) do t_jd_s
        pos_eci = SatelliteToolboxCelestialBodies.sun_position_mod(t_jd_s/3600/24)
        return Vec3f(pos_eci)
    end

    # compute when each target is visible from spacecraft:
    visibilities = Dict()
    for target in target_list
        visibilities[target.name] = visibility_history(target, soln)
    end

    # set up lighting:
    # pl = PointLight(Point3f([4*6371e3;0;0]), RGBf(20, 20, 20))
    # dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), Vec3f(-1, 0, 0))
    dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), sun_light)
    al = AmbientLight(RGBf(0.3, 0.3, 0.3))

    # start main scene:
    ax = LScene(
        fig[1:4,1:3], 
        show_axis=false, 
        scenekw=(lights=[dl, al], 
        backgroundcolor=:black, 
        clear=true)
    )
    # populate auxiliary axes:
    detail_ax = Axis(
        fig[1,4], 
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
        title="Target visibility", 
        xlabel="Time [s]", 
        ylabel="Visible?"
    )
    power_ax = Axis(
        fig[3,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, sim.mission.spacecraft.power.capacity/3600), 
        title="Power", 
        xlabel="Time [s]", 
        ylabel="Battery capacity [Wh]"
    )
    data_ax = Axis(
        fig[4,4],
        backgroundcolor=:black, 
        limits=(0, soln["time"][end] - soln["time"][1], 0, 10), # sim.mission.spacecraft.data.capacity/8e6
        title="Data [PLACEHOLDER]", 
        xlabel="Time [s]", 
        ylabel="Storage [MB]"
    )
    
    # load Earth texture (these are all equirectangular projection/plate carreé):
    texture = load_earth_texture_to_ecef("assets/map_diffuse.png")
    # texture = load_earth_texture_to_ecef("assets/map_bathy.png")
    # texture = load_earth_texture_to_ecef("assets/map_veggie.jpeg")
    
    set_theme!(theme_dark())
    
    plot_earth!(ax, t_jd_s, eops, texture)
    plot_spacecraft!(ax, t_jd_s, 10000, soln)
    plot_targets!(ax, target_list, t_jd_s, eops)
    plot_detail!(detail_ax, t_jd_s, soln)
    plot_visibilities!(visible_ax, t_jd_s, visibilities, soln)
    plot_power!(power_ax, t_jd_s, visibilities, soln)

    fig[2,5] = Legend(fig, visible_ax, "Visibility", framevisible = false)
    linkxaxes!(detail_ax, visible_ax, power_ax)
   
    on(play_button.clicks, priority=1) do n
        if play_button.label == "play"
            
        else
            play_button.label = "stop"
            frame_rate = 60
            @async for i in 1:1000:n_t
                Makie.set_close_to!(time_sliders.sliders[1], i)
                sleep(1/frame_rate)
            end
            play_button.label = "play"
        end
        Consume(true)
    end
end