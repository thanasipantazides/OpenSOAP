import YAML, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies, SatelliteToolboxBase

function load_mission(file::String)
    if !isfile(file)
        error("invalid input file!")
    end
    ext = file[findlast(==('.'), file)+1:end]
    # if ext != "yaml" || ext != "yml"
    #     error("invalid file type: ." * ext)
    # end

    mission_data = YAML.load_file(file)
    gs_file = mission_data["groundstations"]
    sp_file = mission_data["spacecraft"]
    sim_file = mission_data["simulation"]
    sim_spec = load_simulation(sim_file)
    gs_spec = load_groundstations(gs_file, GroundTarget)
    bus_spec = load_spacecraft(sp_file)
    
    mission_data = Mission(
        "impax",
        bus_spec,
        gs_spec
    )

    earth_data = EarthProperties(3.986e14, 1.081874e-4, SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS, 1361)

    tspan = [sim_spec["start_time"]; sim_spec["start_time"] + sim_spec["duration"]]
    leo_sim = LEOSimulation(
        earth=earth_data,
        mission=mission_data,
        tspan=tspan,
        dt=sim_spec["dt"],
        initstate=sim_spec["initial"]
    )
    return leo_sim
end

function load_spacecraft(file::String)
    bus_data = YAML.load_file(file)
    inertia = [bus_data["structure"]["inertia"][1] bus_data["structure"]["inertia"][2] bus_data["structure"]["inertia"][3]]
    mass_data = MassProperties(
        bus_data["structure"]["mass"],
        inertia
    )

    panels = Vector{SolarPanel}(undef, 0)
    for panel in bus_data["power"]["panels"]
        panels = [panels; SolarPanel(panel["normal"], panel["efficiency"], panel["area"])]
    end
    power_data = PowerProperties(
        bus_data["power"]["capacity"],
        bus_data["power"]["base_consumption"],
        panels
    )
    data_data = DataProperties(
        bus_data["data"]["capacity"],
        bus_data["data"]["base_production"],
        bus_data["data"]["base_transmit"]
    )

    wheels = Matrix{Float64}(undef, 3, 0)
    for i in bus_data["attitude"]["wheel_axes"]
        if length(i) == 3
            wheels = hcat(wheels, i)
        else
            wheels = vcat(wheels, i)
        end
    end
    if length(wheels[:,1]) != 3
        wheels = wheels'
    end
    wheels = convert(Matrix{Float64}, wheels)
    attitude_data = AttitudeProperties(
        ReactionWheelProperties(
            wheels,
            bus_data["attitude"]["momenta"],
            bus_data["attitude"]["torques"]
        )
    )
    spacecraft_data = SpacecraftProperties(
        bus_data["name"],
        power_data,
        data_data,
        mass_data,
        attitude_data
    )
    return spacecraft_data
end

function load_simulation(file::String)
    sim_data = YAML.load_file(file)
    time_dict = sim_data["time"]
    start_jd = SatelliteToolboxTransformations.date_to_jd(
        time_dict["start"]["year"],
        time_dict["start"]["month"],
        time_dict["start"]["day"],
        time_dict["start"]["hour"],
        time_dict["start"]["minute"],
        time_dict["start"]["second"],
    )
    start_jd_s = start_jd*24*3600
    duration_s = time_dict["duration"]
    dt_s = time_dict["dt"]

    init_r = sim_data["initial"]["position"]
    init_v = sim_data["initial"]["velocity"]
    init_E = sim_data["initial"]["battery"]
    init_S = sim_data["initial"]["storage"]
    init_ω = sim_data["initial"]["ang_velocity"]
    init_C = diagm([1;1;1])
    init_m = 1.0
    x0 = Float64.([init_r; init_v; init_ω...; init_C...; init_E; init_S; init_m])
    result = Dict("start_time"=>start_jd_s, "duration"=>duration_s, "dt"=>dt_s, "initial"=>x0)
    return result
end

function load_groundstations(file::String, type::Type)

    gs_data = YAML.load_file(file)
    println("\nloading ", length(gs_data), " target records...")

    if type === GroundTarget
        eops = SatelliteToolboxTransformations.fetch_iers_eop()
        result = Vector{GroundTarget}(undef,length(gs_data))
        for i in 1:length(result)
            this_gs = gs_data[i]
            result[i] = GroundTarget(
                this_gs["name"],
                [this_gs["latitude"]*π/180; this_gs["longitude"]*π/180; this_gs["altitude"]],
                [0;0;1],
                (90 - this_gs["elevation"])*π/180,
                eops
            )
        end

        sun = SunTarget(
            "Sun",
            eops
        )
        return [sun; result...]
    end
end
