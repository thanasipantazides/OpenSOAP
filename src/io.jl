import Dates
import JSON
import SatelliteToolboxCelestialBodies, SatelliteToolboxTransformations, SatelliteToolboxBase
import Makie
using GeometryBasics


function lookup_target(s::AbstractString)::Union{AbstractTarget, Nothing}
    if lowercase(s) == "ephemeris.sun"
        return SunState
    elseif lowercase(s) == "igrf"
        warn("IGRF tracking unimplemented, skipping")
        return nothing
    elseif lowercase(s) == "lla"
        return GroundState
    else
        warn("target lookup failed for ", s)
        return nothing
    end
end

function default_mode_list()
    modelist = [
        ModeConfig(0x0101, "idle",      Nothing,    Makie.RGBAf(42/255, 133/255, 255/255), 5.0, 1e3, Vec3d(1.0,0.0,0.0)),
        ModeConfig(0x0102, "charging",  SunState,   Makie.RGBAf(255/255, 201/255, 74/255), 5.0, 1e6, Vec3d(-1.0,0.0,0.0)),
        ModeConfig(0x0103, "downlink",  GroundState, Makie.RGBAf(157/255, 226/255, 107/255), 23.0, 10e3, Vec3d(0.0,1.0,0.0)),
        # ModeConfig(0x0104, "science",   Nothing,    Makie.RGBAf(206/255, 155/255, 255/255), 20.0, 10e3, Vec3d(1.0,0.0,0.0))
    ]
    return modelist
end

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
    for (k,r) in enumerate(remove)
        if k == 1
            push!(keep, 1:r.start-1)
        else
            push!(keep, last.stop+1:r.start-1)
        end
        last = r
    end
    if last.stop+1 < length(raw)
        push!(keep, last.stop+1:length(raw)+1)
    end
    
    # keep the contents in a new string, to parse:
    sout = ""
    for k in keep
        sout *= raw[k]
    end
    
    return JSON.parse(sout)
end

function next_id(id_set::Set{Int64})
    return 1 + max(id_set...)
end

function make_target(json::Dict{String, Any}, start_time::Dates.DateTime, id_registry::Set{Int64})
    name = InlineStrings.InlineString63(json["name"])
    data_source = json["direction"]["source"]
    target_type = lookup_target(data_source)
    if target_type === SunState
        state = SunState(
            0.0, # elapsed time
            1,
            SatelliteToolboxCelestialBodies.sun_position_mod(start_time),
            false,
            false
        )
        # config = TargetConfig(
        #     next_id(id_registry),
            
        # )
        config = nothing
        return (state, config)
    elseif target_type === GroundState
        state = GroundState(
            next_id(id_registry),
            0.0,
            2,
            Point3d(json["direction"]["value"]...),
            SatelliteToolboxTransformations.r_ecef_to_eci(SatelliteToolboxTransformations.ITRF(), SatelliteToolboxTransformations.J2000(), SatelliteToolboxBase.date_to_jd(start_time), eop)*SatelliteToolboxTransformations.geodetic_to_ecef(Point3d(json["direction"]["value"]...)),
            false,
            false
        )
        config = TargetConfig(
            next_id(id_registry),
            json["name"],
            state.id,
            pi/2 - json["elevation_mask"]*pi/180,
            pi
        )
    else
        return nothing
    end
end

function load_config(d::Dict{String, Any})
    
    # todo: set this from config
    start_time = Dates.DateTime(2026,3,20,1,2,3)
    
    
    id_registry = Set{Int64}()
    d["spacecraft"] # nothing to do with this currently, except assign initial conditions.
    
    # first, build list of targets
    for target_data in d["targets"]
        # target_name = InlineStrings.InlineString63(target["name"])
        # target_data_source = target["direction"]["source"]
        target_state = make_target(target_data)
        # do some big lookup of objects based on `target_data_source`
        # Construct an AbstractTarget, then refer to it in TargetConfig
    end
    # then build list of modes, looking up targets
    
    modes = []
    for mode_name in keys(d["modes"])
        # if d["modes"][mode_name]["target"] == ... some lookup in target list
        ModeConfig(
            next_id(id_registry),
            InlineStrings.InlineString63(mode_name),
            SunTarget, # todo: lookup correct target type
            Makie.RGBAf(d["modes"][mode_name]["color"] ./ 255 ...),
            d["modes"][mode_name]["power_consumption"],
            d["modes"][mode_name]["data_production"],
            d["modes"][mode_name]["direction_Body"],
        )
    end
end

# for now, just serve a default mode list.
function load_mode_config(path::AbstractString)
    # d = load_jsonc(path)
    return default_mode_list()
end