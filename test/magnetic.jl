using GLMakie
using LinearAlgebra
using OpenSOAP


function plot_magnetic_field()
    println("loading transformations...")
    eops = SatelliteToolboxTransformations.fetch_iers_eop()
    m = MagneticTarget("mag", eops)
    
    sim = setup_parameters()
    println("integrating...")
    soln = run_orbit(sim)
    n_t = length(soln["time"])
    print(n_t)

    target_list = sim.mission.targets

    GLMakie.activate!(title="OpenSOAP")
    texture = load_earth_texture_to_ecef("assets/map_diffuse.png")
    
    fig = Figure(size=(1400,840))
    display(fig)
    
    # make observable slider:
    time_sliders = SliderGrid(fig[6, 2], (label = "time", format = "{:d}", range = 1:1:n_t, startvalue = 1))
    slider_observables = [s.value for s in time_sliders.sliders]
    
    t_jd_s = lift(slider_observables[1]) do t_i
        return soln["time"][t_i]
    end

    sun_light = lift(t_jd_s) do t_jd_s
        pos_eci = SatelliteToolboxCelestialBodies.sun_position_mod(t_jd_s/3600/24)
        return Vec3f(pos_eci)
    end

    dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), sun_light)
    al = AmbientLight(RGBf(0.3, 0.3, 0.3))

    # start main scene:
    ax = LScene(
        fig[1:5,1:3], 
        show_axis=false, 
        scenekw=(lights=[dl, al], 
        backgroundcolor=:black, 
        clear=true)
    )
    set_theme!(theme_dark())
    
    plot_earth!(ax, t_jd_s, eops, texture)
    plot_spacecraft!(ax, t_jd_s, 10000, soln)
    plot_targets!(ax, target_list, t_jd_s, soln, eops)
    plot_magnetic_field!(ax, t_jd_s, soln, eops)
end