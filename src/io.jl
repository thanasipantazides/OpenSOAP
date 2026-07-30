import Dates
import JSON
import SatelliteToolboxCelestialBodies,
    SatelliteToolboxTransformations, SatelliteToolboxBase, SatelliteToolboxGeomagneticField
import Makie
import DateFormats
using GeometryBasics

function lookup_target_config(
    name::AbstractString,
    target_configs::Vector{TargetConfig},
)::Union{Nothing,IDType}
    res = findfirst(t -> t.name==name, target_configs)
    if isnothing(res)
        return nothing
    else
        return target_configs[res].id
    end
end

function lookup_target(s::T) where {T<:AbstractString}
    if lowercase(s) == "ephemeris.sun"
        return SunState
    elseif lowercase(s) == "igrf"
        # @warn("IGRF tracking unimplemented, skipping")
        return MagneticFieldState
    elseif lowercase(s) == "lla"
        return GroundState
    elseif lowercase(s) == "lla.constraint"
        return LLAConstraint
    else
        @warn("target lookup failed for target type $(s)")
        return nothing
    end
end

function make_mode(
    json,
    mode_name::AbstractString,
    start_time::Dates.DateTime,
    id_registry::Set{IDType},
)
    mode = json["spacecraft"]["modes"][mode_name]
    tstates, tconfigs = fetch_targets_for_mode(json, mode_name, start_time, id_registry)
    constraints = fetch_constraints_for_mode(json, mode_name, start_time, id_registry)

    tstate_ids = [t.id for t in tstates]
    tconfig_ids = [t.id for t in tconfigs]
    constraint_ids = [c.id for c in constraints]

    m = ModeConfig(
        next_id!(id_registry),
        InlineStrings.String63(mode_name),
        tconfig_ids,
        constraint_ids,
        mode["priority"],
        Makie.RGBAf(mode["color"] ./ 255 ...),
        mode["power_consumption"],
        mode["data_production"],
        Vec3d(mode["direction_Body"]),
    )
    return m, tstates, tconfigs, constraints
    # populate the rest of ModeConfig directly from `mode`
end

function fetch_targets_for_mode(
    json,
    mode_name::AbstractString,
    start_time::Dates.DateTime,
    id_registry::Set{IDType},
)
    mode = json["spacecraft"]["modes"][mode_name]

    target_states = Vector{AbstractTarget}()
    target_configs = Vector{AbstractConfig}()

    if "target" in keys(mode) && length(mode["target"]) > 0
        for tstr in mode["target"]
            if tstr != ""
                targetj = find_by_name(json["spacecraft"]["targets"], tstr)

                if !isnothing(targetj)
                    t, c = make_target(targetj, start_time, id_registry)
                    if !isnothing(t)
                        push!(target_states, t)
                        push!(target_configs, c)
                    end
                end
            end
        end
    end
    return target_states, target_configs
end

function fetch_constraints_for_mode(
    json,
    mode_name::AbstractString,
    start_time::Dates.DateTime,
    id_registry::Set{IDType},
)
    mode = json["spacecraft"]["modes"][mode_name]

    constraints = Vector{AbstractConstraint}()
    if "constraint" in keys(mode) && length(mode["constraint"]) > 0
        for cstr in mode["constraint"]
            if cstr != ""
                constraintj = find_by_name(json["spacecraft"]["constraints"], cstr)

                if !isnothing(constraintj)
                    c = make_constraint(constraintj, start_time, id_registry)
                    if !isnothing(c)
                        push!(constraints, c)
                    end
                end
            end
        end
    end
    return constraints
end

"""
    find_by_name(json, name::AbstractString)

Read a JSON array and return the object containing a `key == name`.
"""
function find_by_name(json, name::AbstractString; key = "name")
    for field in json
        if key in keys(field) && field[key] == name
            return field
        end
    end
    return nothing
end

function default_mode_list()
    modelist = [
        ModeConfig(
            0x0101,
            "idle",
            [],
            [],
            3,
            Makie.RGBAf(42/255, 133/255, 255/255),
            5.0,
            1e3,
            Vec3d(1.0, 0.0, 0.0),
        ),
        ModeConfig(
            0x0102,
            "charging",
            [],
            [],
            2,
            Makie.RGBAf(255/255, 201/255, 74/255),
            5.0,
            1e6,
            Vec3d(-1.0, 0.0, 0.0),
        ),
        ModeConfig(
            0x0103,
            "downlink",
            [],
            [],
            1,
            Makie.RGBAf(157/255, 226/255, 107/255),
            23.0,
            10e3,
            Vec3d(0.0, 1.0, 0.0),
        ),
        # ModeConfig(0x0104, "science",   Nothing,    Makie.RGBAf(206/255, 155/255, 255/255), 20.0, 10e3, Vec3d(1.0,0.0,0.0))
    ]
    return modelist
end

"""
    load_jsonc(path::AbstractString)

Load a JSON file with comments, then parse using JSON.jl.
"""
function load_jsonc(path::AbstractString)

    raw = read(path, String)
    remove = UnitRange{Int64}[]
    nodouble = nostart = false
    k = 1
    while k < length(raw)
        # get start of next comment.
        kstart = findnext("//", raw, k)
        if isnothing(kstart) # try /* */ style comments instead.
            kstart = findnext("/*", raw, k)
            if !isnothing(kstart)
                kend = findnext("*/", raw, kstart[1])
                if isnothing(kend) # file ends in this comment, no newline
                    kend = length(raw)
                end
                push!(remove, UnitRange{Int64}(kstart[1], max(kend)))
                k = kend+1
            else # no more comments
                break
            end
        else
            kend = findnext('\n', raw, kstart[2])
            if isnothing(kend) # file ends in this comment, no newline
                kend = length(raw)
            end
            # remove from the start of the comment to kend
            push!(remove, UnitRange{Int64}(kstart[1], max(kend)-1))
            k = kend+1
        end
    end

    # get complement of the ranges to remove to see what to keep
    keep = UnitRange{Int64}[]
    last = UnitRange{Int64}(-1:0)
    for (k, r) in enumerate(remove)
        if k == 1
            push!(keep, 1:(r.start-1))
        else
            push!(keep, (last.stop+1):(r.start-1))
        end
        last = r
    end
    # N.B. the sizeof(str) != length(str) due to Unicode. Watch out!
    if last.stop+1 < sizeof(raw)
        push!(keep, (last.stop+1):(sizeof(raw)))
    end

    # keep the contents in a new string, to parse:
    sout = ""
    for k in keep
        sout *= raw[k]
    end

    return JSON.parse(sout)
end

function next_id!(id_set::Set{IDType})
    next_id = 1 + max(id_set...)
    push!(id_set, next_id)
    return next_id
end

function make_target(json, start_time::Dates.DateTime, id_registry::Set{IDType})
    name = InlineStrings.InlineString63(json["name"])
    data_source = json["direction"]["source"]
    target_type = lookup_target(data_source)
    if target_type === SunState
        state = SunState(
            next_id!(id_registry),
            0.0, # elapsed time
            1,
            SatelliteToolboxCelestialBodies.sun_position_mod(start_time),
            false,
            false,
        )
        # config = TargetConfig(
        #     next_id!(id_registry),

        # )
        config = TargetConfig(next_id!(id_registry), json["name"], state.id, 0.0, pi, pi/2)
        return (state, config)
    elseif target_type === GroundState
        state = GroundState(
            next_id!(id_registry),
            0.0,
            2,
            pi/180*Point3d(json["direction"]["value"]...),
            SatelliteToolboxTransformations.r_ecef_to_eci(
                SatelliteToolboxTransformations.ITRF(),
                SatelliteToolboxTransformations.J2000(),
                SatelliteToolboxBase.date_to_jd(start_time),
                eop,
            )*SatelliteToolboxTransformations.geodetic_to_ecef(
                pi/180*Point3d(json["direction"]["value"]...),
            ),
            false,
            false,
        )
        config = TargetConfig(
            next_id!(id_registry),
            json["name"],
            state.id,
            json["data_consumption"],
            pi/2 - json["elevation_mask"]*pi/180,
            pi,
        )
        return (state, config)
    elseif target_type === MagneticFieldState

        # create a fake initial position to get a reasonable starting magnetic field value out.
        init_ECI = [6371e3, 0.0, 0.0]
        C_ECEF_ECI = SatelliteToolboxTransformations.r_eci_to_ecef(
            SatelliteToolboxTransformations.J2000(),
            SatelliteToolboxTransformations.ITRF(),
            SatelliteToolboxBase.date_to_jd(start_time),
            eop,
        )
        init_lla = SatelliteToolboxTransformations.ecef_to_geodetic(C_ECEF_ECI*init_ECI)

        state = MagneticFieldState(
            next_id!(id_registry),
            0.0,
            Vec3d(
                C_ECEF_ECI'*SatelliteToolboxTransformations.ned_to_ecef(
                    SatelliteToolboxGeomagneticField.igrf(
                        DateFormats.yeardecimal(start_time),
                        6371e3,
                        0.0,
                        0.0,
                        Val(:geodetic);
                        max_degree = clamp(json["direction"]["model_order"], 1, 13),
                        show_warnings = false,
                    ),
                    init_lla...,
                    translate = false,
                ),
            ),
            false,
            false,
            false,
        )

        config=MagneticFieldConfig(
            next_id!(id_registry),
            json["name"],
            state.id,
            lookup_igrf_normalization(json["direction"]["normalization"]),
            clamp(json["direction"]["model_order"], 1, 13),
        )
        return (state, config)
    else
        # (state, config)
        return (nothing, nothing)
    end
end

function make_constraint(json, start_time::Dates.DateTime, id_registry::Set{IDType})
    # return all "constraint": objects. Expect to be passed a complete target.

    # a target with no constraints
    if !("source" in keys(json))
        return nothing
    end
    data_source = json["source"]
    constraint_type = lookup_target(data_source)
    println("Got constraint type: $constraint_type")

    if constraint_type === LLAConstraint
        lat = [isa(val, Real) ? val : parse(Float64, val) for val in json["value"]["lat"]]
        lon = [isa(val, Real) ? val : parse(Float64, val) for val in json["value"]["lon"]]
        alt = [isa(val, Real) ? val : parse(Float64, val) for val in json["value"]["alt"]]

        lat2 = pi/180*Point2d(min(lat...), max(lat...))
        lon2 = pi/180*Point2d(min(lon...), max(lon...))
        alt2 = pi/180*Point2d(min(alt...), max(alt...))

        c = LLAConstraint(next_id!(id_registry), lat2, lon2, alt2)
        return c
    else
        return nothing
    end
end

function lookup_target_state(
    name::AbstractString,
    target_states::Vector{Union{Nothing,AbstractTarget}},
    target_configs::Vector{TargetConfig},
)
    # compare the `name` to the .name value in target config, return the state that matches.

    for state in target_states
        # find the config for this state:
        config = nothing
        for this_config in target_configs
            if !isnothing(this_config) && !isnothing(state)
                if state.id == this_config.dynamic_id
                    config = this_config
                    break
                end
            end
        end
        # config = findfirst(config -> config.dynamic_id == state.id, target_configs)
        if isnothing(config)
            println("no config for this target!")
        elseif name == config.name
            return state
        else
            continue
        end
    end
    return nothing
end

"""
    collapse_panels(efficiencies::Vector{Float64}, areas::Vector{Float64}, normals::Vector{Vector{Float64}})
    
Return the solar panel with maximum area*efficiency product. 

Plan to deprecate after vectorizing `SatelliteConfig` solar panel data.
"""
function collapse_panels(
    efficiencies::Vector{Float64},
    areas::Vector{Float64},
    normals::Vector{Vector{Float64}},
)
    quality = efficiencies .* areas
    m = findmax(quality)[2]
    return (efficiencies[m], areas[m], normals[m])
end

function get_json_mat3d(major_key::AbstractString, val::AbstractVector)
    value = Mat3d(I)
    if lowercase(major_key) == "row"
        value = vcat(val'...)
    elseif lowercase(major_key) == "column"
        value = hcat(val...)
    else
        @warn "undefined key $(major_key)!"
    end
    return Mat3d(value)
end

function make_mode_table(
    modes::Vector{ModeConfig},
    target_configs::Vector{<:AbstractConfig},
)::Matrix{UInt8}
    M = zeros(UInt8, length(target_configs), length(modes))
    for (m, mode) in enumerate(modes)
        for (t, target) in enumerate(target_configs)
            if target.id in mode.target_ids
                M[t, m] = 1
            end
        end
    end
    return M
end

function load_config(d::Dict{String,Any})
    id_registry = Set{IDType}(23)

    siminittime = d["simulation"]["initial"]["time"]
    sim = SimConfig(
        next_id!(id_registry),
        Dates.DateTime(
            siminittime["year"],
            siminittime["month"],
            siminittime["day"],
            siminittime["hour"],
            siminittime["minute"],
            siminittime["second"],
        ),
        d["simulation"]["timestep"],
        0,
        Dict(),
    )

    initial_attitude = get_json_mat3d(
        d["simulation"]["initial"]["attitude_major"],
        d["simulation"]["initial"]["attitude"],
    )

    sat = SatelliteState(
        next_id!(id_registry),      # ID
        0.0,                        # met
        Vec3d(0.0),     # net force
        Vec3d(0.0),     # net moment
        Vec3d(d["simulation"]["initial"]["position"]),
        Vec3d(d["simulation"]["initial"]["velocity"]),
        Vec3d(d["simulation"]["initial"]["angular_velocity"]),
        initial_attitude,           # attitude
        0.0,
        0.0,                   # net power, net data flows
        d["simulation"]["initial"]["battery"],
        d["simulation"]["initial"]["storage"], # battery, data storage
        IDType(0x00),               # mode
        IDType(0x00),               # target  
        false,
        false,                # target visible, pointed
    )

    inertia = get_json_mat3d(
        d["spacecraft"]["mechanical"]["inertia_major"],
        d["spacecraft"]["mechanical"]["inertia"],
    )

    eff = Vector{Float64}()
    area = Vector{Float64}()
    normal = Vector{Vector{Float64}}()
    for item in d["spacecraft"]["power"]["panels"]
        push!(eff, item["efficiency"])
        push!(area, item["area"])
        push!(normal, item["normal"])
    end

    p = collapse_panels(eff, area, normal)

    # todo: don't hard index this:
    antenna_dir = Vec3d(d["spacecraft"]["data"]["antennas"][1]["normal"])

    sat_config = SatelliteConfig(
        next_id!(id_registry),
        d["spacecraft"]["name"],
        sat.id,
        inertia,
        inv(inertia),
        d["spacecraft"]["control"]["max_slew_rate"],
        d["spacecraft"]["power"]["capacity"],
        d["spacecraft"]["data"]["capacity"],
        p[1],
        p[2],
        p[3],
        antenna_dir,
    )

    # target_states = Vector{AbstractTarget}()
    # target_config = Vector{TargetConfig}()
    # target_constraints = Vector{AbstractConstraint}()

    # # first, build list of targets
    # targets = d["spacecraft"]["targets"]
    # for target_data in targets
    #     this_target_state, this_target_config =
    #         make_target(target_data, sim.start_time, id_registry)
    #     constraints = make_constraint(target_data, sim.start_time, id_registry)
    #     if !isnothing(this_target_state)
    #         push!(target_states, this_target_state)
    #         push!(target_config, this_target_config)
    #     end
    #     # if length(constraints) > 0
    #     #     append!(target_constraints, constraints)
    #     # end

    #     # do some big lookup of objects based on `target_data_source`
    #     # Construct an AbstractTarget, then refer to it in TargetConfig
    # end

    # then build list of modes, looking up targets
    modes = ModeConfig[]
    target_states = Set{AbstractState}()
    target_configs = Set{AbstractConfig}()
    constraints = Set{AbstractConstraint}()

    for mode_name in keys(d["spacecraft"]["modes"])

        # make the mode. If a target or target config or constraint needs making, do it.


        # note: target states, target configs, and constraints need to get returned from this too so they don't just disappear
        mode, new_states, new_configs, new_constraints =
            make_mode(d, mode_name, sim.start_time, id_registry)
        push!(modes, mode)
        union!(target_states, new_states)
        union!(target_configs, new_configs)
        union!(constraints, new_constraints)

        # mode_target_names = d["spacecraft"]["modes"][mode_name]["target"]
        # mode_target_ids =
        #     [lookup_target_config(name, target_config) for name in mode_target_names]
        # # println("> for mode $(mode_name): $(mode_target_ids)")
        # # this is ok to do because target_config don't contain `nothing`
        # mode_targets = filter(x -> x.id in mode_target_ids, target_config)

        # # todo: filter target_constraints and insert for this mode.

        # m = ModeConfig(
        #     next_id!(id_registry),
        #     InlineStrings.InlineString63(mode_name),
        #     # SunTarget, # todo: lookup correct target type
        #     # typeof(lookup_target_state(d["spacecraft"]["modes"][mode_name]["target"], target_states, target_config)),
        #     map(x->x.id, mode_targets),
        #     [], # todo: insert constraints here.
        #     d["spacecraft"]["modes"][mode_name]["priority"],
        #     Makie.RGBAf(d["spacecraft"]["modes"][mode_name]["color"] ./ 255 ...),
        #     d["spacecraft"]["modes"][mode_name]["power_consumption"],
        #     d["spacecraft"]["modes"][mode_name]["data_production"],
        #     d["spacecraft"]["modes"][mode_name]["direction_Body"],  # todo: normalize
        # )

        # push!(modes, mode)
        # if length(intersect(mode_target_ids, map(x->x.id, target_config))) != length(mode_target_ids)
        #     @warn "Some targets missing for mode $(m.id) : $(mode_name). They have been dropped."
        # end
    end

    # convert these to vectors now that they are done being populated
    target_configs = Vector([target_configs...])
    target_states = Vector([target_states...])
    constraints = Vector([constraints...])

    targetless = filter(m -> length(m.target_ids) == 0, modes)
    if length(targetless) != 1
        @error "Can only have one mode without any targets listed (idle)!"
    end
    mode_table = make_mode_table(modes, target_configs)

    sat.mode = modes[1].id
    sat.target = target_states[1].id

    return (sim, sat, sat_config, target_states, target_configs, modes, mode_table)
end

# for now, just serve a default mode list.
function load_mode_config(path::AbstractString)
    return default_mode_list()
end


function tabulate(
    x::Union{Vector{<:AbstractTarget},Vector{<:AbstractState},Vector{<:AbstractConfig}},
)
    return Dict(p.id => p for p in x)
end
