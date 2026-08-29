
# an isbitstype(IDVector) == true wrapper, to enable good serialization of 
# lists of IDs (e.g. for ModeConfig.)
struct IDVector{N} <: AbstractVector{IDType}
    data::SVector{N,IDType}
    length::IDType          # not fully generic: IDType is meant for unique 
    # objects, so only need to store up to typemax(IDType) of them.
end

Base.size(v::IDVector) = (v.length,)
Base.getindex(v::IDVector, k) = v.data[k]
Base.sizeof(v::IDVector{N}) where {N} = N

function IDVector{N}(v::AbstractVector) where {N}
    if N > SOAP_MAX_ID_BUCKET
        error("N is too large")
    end
    n = length(v)
    if n > N
        error("v is too large!")
    end

    return IDVector{N}(SVector{N}(vcat(v, zeros(IDType, N - n))), n)
end

function IDVector(v::AbstractVector)
    if length(v) > SOAP_MAX_ID_BUCKET
        error("N is too large")
    end
    n = length(v)

    return IDVector{SOAP_MAX_ID_BUCKET}(
        SVector{SOAP_MAX_ID_BUCKET}(vcat(v, zeros(IDType, SOAP_MAX_ID_BUCKET - n))),
        n,
    )
end

# an isbitstype(SizedString) == true wrapper, to enable serialization of strings
struct SizedString{N} <: AbstractString
    data::SVector{N,UInt8}
    length::UInt16
end

function SizedString{N}(s::AbstractString) where {N}
    if N > typemax(UInt16)
        throw(ArgumentError("FixedString{N} requires N <= 65535"))
    end
    bytes = codeunits(s)
    n = length(bytes)
    if n > N
        throw(ArgumentError("string $s exceeds FixedString capacity $N"))
    end
    padded = ntuple(i -> i <= n ? bytes[i] : 0x00, N)
    SizedString{N}(padded, UInt8(n))
end

function SizedString(s::AbstractString)
    if length(codeunits(s)) > SOAP_MAX_STRING_LEN
        throw(
            ArgumentError(
                "FixedString requires a length (in codeunits) < $SOAP_MAX_STRING_LEN",
            ),
        )
    end
    bytes = codeunits(s)
    n = length(bytes)
    padded = ntuple(i -> i <= n ? bytes[i] : 0x00, SOAP_MAX_STRING_LEN)
    SizedString{SOAP_MAX_STRING_LEN}(padded, UInt8(n))
end

# AbstractString interface:
Base.ncodeunits(s::SizedString) = Int(s.length)
Base.codeunit(::SizedString) = UInt8
Base.codeunit(s::SizedString, i::Integer) = s.data[i]
Base.isvalid(s::SizedString, i::Integer) = isvalid(String(s), i)
Base.String(s::SizedString) = String(UInt8[s.data[i] for i = 1:s.length])

Base.show(io::IO, s::SizedString) = show(io, String(s))
Base.:(==)(a::SizedString, b::SizedString) = a.length == b.length && a.data == b.data
Base.hash(s::SizedString, h::UInt) = hash(s.data, hash(s.length, h))

capacity(::SizedString{N}) where {N} = N
capacity(::Type{SizedString{N}}) where {N} = N

function Base.iterate(s::SizedString, i::Int = 1)
    i > s.length && return nothing
    c, j = iterate(String(s), i)  # delegate UTF-8 decoding to String
    return c, j
end

function Base.:(==)(a::SizedString, b::AbstractString)
    bb = codeunits(b)
    a.length == length(bb) || return false
    @inbounds for i = 1:a.length
        a.data[i] == bb[i] || return false
    end
    return true
end

function GeometryBasics.Mat3d(x::Float64)
    return Mat3d(fill(x, (3, 3)))
end

# the "*State" structs should be dynamic. 
# All static config variables should be stored elsewhere.

"""
    NetworkMessage

An abstract supertype for types that can be serialized to the network.
"""
abstract type NetworkMessage end

"""
    ControlMessage<:NetworkMessage

An abstract type used for controlling the simulation speed or execution.
"""
abstract type ControlMessage<:NetworkMessage end

"""
    AskMessage<:ControlMessage

A message used to ask the core simulation for a value.
"""
@kwdef mutable struct AskMessage<:ControlMessage
    message::IDType = 0x00
end
"""
    PlayMessage<:ControlMessage

A message used to pause or play the simulation. A value of 0x00 is paused.
"""
@kwdef mutable struct PlayMessage<:ControlMessage
    message::UInt8 = 0x00
end

"""
    RateMessage<:ControlMessage

A message used to speed up or slow down the simulation.
"""
@kwdef mutable struct RateMessage<:ControlMessage
    message::UInt8 = 0x00
end

@kwdef mutable struct QuitMessage<:ControlMessage
    message::UInt8 = 0x00
end

@kwdef mutable struct PerturbationMessage<:ControlMessage
    moment_Body::Vec3d = Vec3d(0.0)
    moment_duration::UInt32 = 0x00 # number of counter steps to sustain
    force_ECI::Vec3d = Vec3d(0.0)
    force_duration::UInt32 = 0x00# number of counter steps to sustain
end

# some static configuration for stuff the mission wants to look at
abstract type AbstractConfig<:NetworkMessage end
# something with a time-dependent state
abstract type AbstractState<:NetworkMessage end
# dynamic configuration of something that the mission wants to look at
abstract type AbstractTarget<:AbstractState end

abstract type AbstractConstraint<:NetworkMessage end

# server-side, the table of modes can live in params.
# But! want a way to construct the mode table via sockets.
struct ModeConfig<:AbstractConfig
    id::IDType
    name::SizedString{SOAP_MAX_STRING_LEN}    # todo: method for reinterpret on .String63 that actually works
    # target_type::DataType
    target_ids::IDVector{SOAP_MAX_ID_BUCKET}      # lookup TargetConfig by ID.
    constraint_ids::IDVector{SOAP_MAX_ID_BUCKET}  # lookup <:AbstractConstraint by ID.
    priority::IDType                # low value => high priority; high value => low priority
    color::Makie.RGBAf
    power_consumption::Float64
    data_production::Float64
    direction_Body::Vec3d           # body vector to point at the target
end

struct TargetConfig<:AbstractConfig
    id::IDType
    name::SizedString{SOAP_MAX_STRING_LEN}
    dynamic_id::IDType
    # color::Makie.RGBAf
    data_consumption::Float64
    position_cone::Float64  # must contain the spacecraft to get the target
    pointing_cone::Float64  # spacecraft sensor must put target inside this cone
end

struct MagneticFieldConfig<:AbstractConfig
    id::IDType
    name::SizedString{SOAP_MAX_STRING_LEN}
    dynamic_id::IDType
    normalization::UInt8
    model_order::UInt8
end

# placeholder
struct EarthConfig<:AbstractConfig
    id::IDType
    name::SizedString{SOAP_MAX_STRING_LEN}
    dynamic_id::IDType
end

@kwdef struct SatelliteConfig<:AbstractConfig
    id::IDType
    name::SizedString{SOAP_MAX_STRING_LEN}
    dynamic_id::IDType
    inertia_Body::Mat3d
    inertia_inv_Body::Mat3d
    surface_area::Float64       # total surface area
    angular_rate_max::Float64
    power_battery_max::Float64
    data_storage_max::Float64
    power_solar_panel_efficiency::Float64
    power_solar_panel_area::Float64
    power_solar_panel_direction_Body::Vec3d
    data_antenna_direction_Body::Vec3d
end

@kwdef mutable struct SimConfig<:AbstractConfig
    id::IDType
    start_time::Dates.DateTime
    time_step::Float64
    step_count::UInt64
    kw::Dict{String,Any}
end

# find the config struct for a corresponding dynamic (state) struct by id
function find_config(
    state::AbstractState,
    config_lookup::Vector{T},
) where {T<:NetworkMessage}
    res = findfirst(
        p -> p isa AbstractConfig && !(p isa ModeConfig) && p.dynamic_id == state.id,
        config_lookup,
    )
    return config_lookup[res]
end

# find the config struct for a corresponding dynamic (state) struct by id
function find_config(
    state::AbstractState,
    config_lookup::IDDict{T},
) where {T<:NetworkMessage}
    res = config_lookup[findfirst(
        p -> p isa AbstractConfig && !(p isa ModeConfig) && p.dynamic_id == state.id,
        config_lookup,
    )]
    return res
end

function filtertype(T::DataType, sim_environment::IDDict{S}) where {S}
    return filter(s -> typeof(s.second) <: T, sim_environment)
end

@kwdef mutable struct PositionState<:AbstractState
    elapsed_time::Float64 = 0.0
    position_ECI::Point3d = Point3d(0.0)
    velocity_ECI::Point3d = Point3d(0.0)
end

@kwdef mutable struct AttitudeState<:AbstractState
    elapsed_time::Float64 = 0.0
    angular_velocity_ECI_Body::Vec3d = Vec3d(0.0)
    attitude_ECI_Body::Mat3d = Mat3d(I(3))
end

# note: with separate Pos and Att states, 
# incident rate on monitor was ~2-4 MBps
# incident sim rate was ~6300:1
mutable struct SatelliteState<:AbstractState
    const id::IDType
    elapsed_time::Float64

    # dynamic properties
    mass::Float64
    net_force_ECI::Vec3d
    net_moment_Body::Vec3d
    position_ECI::Vec3d
    velocity_ECI::Vec3d
    angular_velocity_ECI_Body::Vec3d
    attitude_ECI_Body::Mat3d

    # system properties
    net_power::Float64
    net_data::Float64
    battery_level::Float64
    storage_level::Float64

    # operations properties
    mode::IDType
    target::IDType
    target_visible::UInt8
    target_pointed::UInt8

    # consider adding private reference to config data, targets, etc.
end

function mode_table(modes::Vector{ModeConfig})::Dict{IDType,ModeConfig}
    d = Dict{IDType,ModeConfig}(m.id => m for m in modes)
    return d
end



# empty constructor
SatelliteState() = SatelliteState(
    0x0000,
    0.0,
    0.0,
    Vec3d(0),
    Vec3d(0),
    Vec3d(0),
    Vec3d(0),
    Vec3d(0),
    Mat3d(I),
    0.0,
    0.0,
    0.0,
    0.0,
    0xffff,
    0xffff,
    0x00,
    0x00,
)

@kwdef mutable struct EarthState<:AbstractTarget
    id::IDType = 0x00
    elapsed_time::Float64 = 0.0
    attitude_ECI_ECEF::Mat3d = Mat3d(I(3))
end

@kwdef mutable struct SunState<:AbstractTarget
    const id::IDType = 0x00
    elapsed_time::Float64 = 0.0
    priority::UInt16 = 0x00# todo: make this a Mode priority, not Target priority.
    position_ECI::Point3d = Point3d(0.0)
    visible::Bool = false
    selected::Bool = false
end

@kwdef mutable struct GroundState<:AbstractTarget
    const id::IDType = 0x00
    elapsed_time::Float64 = 0.0

    priority::UInt16 = 0x00 # todo: make this a Mode priority, not Target priority.
    position_LLA::Point3d = Point3d(0.0)
    position_ECI::Point3d = Point3d(0.0)
    visible::Bool = false
    selected::Bool = false
end

@kwdef mutable struct MagneticFieldState<:AbstractTarget
    const id::IDType = 0x00
    elapsed_time::Float64 = 0.0

    direction_ECI::Vec3d = Vec3d(0.0)
    visible::Bool = false
    available::Bool = false
    selected::Bool = false
end

struct LLAConstraint<:AbstractConstraint
    id::IDType
    name::SizedString{SOAP_MAX_STRING_LEN}
    lat::Point2d
    lon::Point2d
    alt::Point2d
end

function reference_direction(target::T) where {T<:Union{SunState,GroundState}}
    return target.position_ECI
end

function reference_direction(target::MagneticFieldState)
    return target.direction_ECI
end

function +(a::PositionState, b::PositionState)
    if a.elapsed_time != b.elapsed_time
        error("can only add States with the same time")
    end
    return PositionState(
        a.elapsed_time,
        a.position_ECI + b.position_ECI,
        a.velocity_ECI + b.velocity_ECI,
    )
end

function *(c::Real, a::PositionState)
    return PositionState(a.elapsed_time, c*a.position_ECI, c*a.velocity_ECI)
end

function +(a::AttitudeState, b::AttitudeState)
    if a.elapsed_time != b.elapsed_time
        error("can only add States with the same time")
    end
    # linear only in angular rate
    return AttitudeState(
        a.elapsed_time,
        a.angular_velocity_ECI_Body + b.angular_velocity_ECI_Body,
        a.attitude_ECI_Body,
    )
end

function *(c::Real, a::AttitudeState)
    # linear only in angular rate
    return AttitudeState(a.elapsed_time, c*a.angular_velocity_ECI_Body, a.attitude_ECI_Body)
end

"""
    mut_struct_eq(ls, rs)::Bool

check for equality of all fields of mutable structs. 

Note that for two mutable structs T and S, T == S is always false.
"""
function mut_struct_eq(ls, rs)::Bool
    return all(getproperty(ls, f) == getproperty(rs, f) for f in propertynames(ls))
end

function lookup_igrf_normalization(value::AbstractString)::UInt8
    d = Dict{String,UInt8}("none" => 0x00, "nadir" => 0x01, "zenith" => 0x02)
    if value in keys(d)
        return d[value]
    else
        @warn "key $value not a valid IGRF normalization"
        return d["none"]
    end
end

function IDDict(
    val::Vector{T},
) where {T<:Union{AbstractState,AbstractConfig,AbstractConstraint}}
    return Dict(v.id => v for v in val)
end

function IDDict(
    val::Set{T},
) where {T<:Union{AbstractState,AbstractConfig,AbstractConstraint}}
    return Dict(v.id => v for v in val)
end

function Base.show(io::IO, confs::Vector{<:AbstractConfig})
    fmt =
        PrettyTables.TextTableFormat(borders = PrettyTables.text_table_borders__borderless)

    els = fill("", (length(confs), 3))
    for (k, c) in enumerate(confs)
        els[k, 1] = "$(c.name)"
        els[k, 2] = "$(string(c.id, base=16, pad=4))"
        els[k, 3] = c isa ModeConfig ? "" : "-> $(string(c.dynamic_id, base=16, pad=4))"
    end
    PrettyTables.pretty_table(
        io,
        els;
        title = "Target configurations",
        column_labels = ["Name", "Conf. ID", "Dyn. ID"],
        backend = :text,
        table_format = fmt,
    )
end

function Base.show(io::IO, confs::IDDict{<:AbstractConfig})
    fmt =
        PrettyTables.TextTableFormat(borders = PrettyTables.text_table_borders__borderless)

    els = fill("", (length(confs), 3))
    for (k, c) in enumerate(collect(keys(confs)))
        els[k, 1] = "$(confs[c].name)"
        els[k, 2] = "$(string(confs[c].id, base=16, pad=4))"
        els[k, 3] =
            c isa ModeConfig ? "" : "-> $(string(confs[c].dynamic_id, base=16, pad=4))"
    end
    PrettyTables.pretty_table(
        io,
        els;
        title = "Target configurations",
        column_labels = ["Name", "Conf. ID", "Dyn. ID"],
        backend = :text,
        table_format = fmt,
    )
end

function Base.show(
    io::IO,
    data::Tuple{
        Vector{ModeConfig},
        IDDict{<:AbstractState},
        IDDict{<:AbstractConfig},
        IDDict{<:AbstractConstraint},
    },
)
    fmt =
        PrettyTables.TextTableFormat(borders = PrettyTables.text_table_borders__borderless)

    modes = data[1]
    target_configs = data[3]
    constraints = data[4]

    all_objs = cat(collect(values(target_configs)), collect(values(constraints)), dims = 1)
    all_names = [c.name for c in all_objs]

    mode_names = [m.name for m in modes]
    push!(mode_names, "target ID")
    push!(mode_names, "type")

    els = fill("---", (length(target_configs) + length(constraints), length(modes) + 2))
    for m in eachindex(modes)
        trow = findall(c -> c.id in modes[m].target_ids, collect(values(target_configs)))
        crow =
            length(target_configs) .+
            findall(c -> c.id in modes[m].constraint_ids, collect(values(constraints)))
        els[trow, m] .= "-x-"
        els[crow, m] .= "-x-"
    end

    for r in eachindex(all_objs)
        if typeof(all_objs[r]) <: AbstractConfig
            els[
                r,
                length(modes)+1,
            ] = "$(string(all_objs[r].id, base=16, pad=4)) -> $(string(all_objs[r].dynamic_id, base=16, pad=4))"
        else
            els[r, length(modes)+1] = ""
        end
        els[r, length(modes)+2] = string(nameof(typeof(all_objs[r])))
    end

    PrettyTables.pretty_table(
        io,
        els;
        title = "Modes",
        column_labels = mode_names,
        row_labels = all_names,
        alignment = :c,
        backend = :text,
        table_format = fmt,
    )
end
