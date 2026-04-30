import Serialization
import Sockets
import InlineStrings
using GeometryBasics
import Dates
import SatelliteToolboxTransformations, SatelliteToolboxBase, SatelliteToolboxCelestialBodies
import Base: +, *
import Makie

global const unixsockname = "/tmp/misik_out.sock"
global const IDType = UInt16
# global const unixsock"tname = "127.0.0.1"

# the "*State" structs should be dynamic. 
# All static config variables should be stored elsewhere.

abstract type NetworkMessage end

abstract type ControlMessage<:NetworkMessage end

# control pausing/playing the core simulation
mutable struct PlayMessage<:ControlMessage
    message::UInt8
end

# control the run rate of the core simulation
mutable struct RateMessage<:ControlMessage
    message::UInt8
end

# some static configuration for stuff the mission wants to look at
abstract type AbstractConfig<:NetworkMessage end
# something with a time-dependent state
abstract type AbstractState<:NetworkMessage end
# dynamic configuration of something that the mission wants to look at
abstract type AbstractTarget<:AbstractState end
    
# server-side, the table of modes can live in params.
# But! want a way to construct the mode table via sockets.
struct ModeConfig<:AbstractConfig
    id::IDType
    name::InlineStrings.String63    # todo: method for reinterpret on .String63 that actually works
    target_type::DataType
    color::Makie.RGBAf
    power_consumption::Float64
    data_production::Float64
    direction_Body::Vec3d           # body vector to point at the target
end

struct TargetConfig<:AbstractConfig
    id::IDType
    name::InlineStrings.String63
    dynamic_id::IDType
    color::Makie.RGBAf
end

# find the config struct for a corresponding dynamic (state) struct by id
function find_config(from::AbstractState, lookup::Vector{T}) where T<:AbstractConfig
    res = find(p -> p.id == from.dynamic_id, lookup)
    return res
end

mutable struct PositionState<:AbstractState
    elapsed_time::Float64
    position_ECI::Point3d
    velocity_ECI::Point3d
end

mutable struct AttitudeState<:AbstractState
    elapsed_time::Float64
    angular_velocity_ECI_Body::Vec3d
    attitude_ECI_Body::Mat3d
end

# note: with separate Pos and Att states, 
# incident rate on monitor was ~2-4 MBps
# incident sim rate was ~6300:1
mutable struct SatelliteState<:AbstractState
    const id::IDType
    elapsed_time::Float64
    
    # dynamic properties
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
    mode::UInt16
    target_ECI::Vec3d
    target_visible::UInt8
    target_pointed::UInt8
    
    # consider adding private reference to config data, targets, etc.
end

function mode_table(modes::Vector{ModeConfig})::Dict{IDType, ModeConfig}
    d = Dict{IDType, ModeConfig}(m.id => m for m in modes)
    return d
end



# empty constructor
SatelliteState() = SatelliteState(0x0000, 0.0, Vec3d(0), Vec3d(0), Vec3d(0), Vec3d(0),Vec3d(0), Mat3d(I), 0.0,0.0,0.0,0.0, 0xffff, Vec3d(0.0), 0x00, 0x00)

mutable struct EarthState<:AbstractTarget
    elapsed_time::Float64
    attitude_ECI_ECEF::Mat3d
end

mutable struct SunState<:AbstractTarget
    elapsed_time::Float64
    priority::UInt16
    position_ECI::Point3d
    visible::Bool
    selected::Bool
end

# todo: rename "GroundState" or something to convey a const LLA value (rather than a generic target)
mutable struct TargetState<:AbstractTarget
    const id::IDType
    elapsed_time::Float64
    
    priority::UInt16
    const position_LLA::Point3d
    position_ECI::Point3d
    visible::Bool
    selected::Bool
end

function +(a::PositionState, b::PositionState)
    if a.elapsed_time != b.elapsed_time
        error("can only add States with the same time")
    end
    return PositionState(
        a.elapsed_time, 
        a.position_ECI + b.position_ECI, 
        a.velocity_ECI + b.velocity_ECI
    )
end

function *(c::Real, a::PositionState)
    return PositionState(
        a.elapsed_time,
        c*a.position_ECI,
        c*a.velocity_ECI
    )
end

function +(a::AttitudeState, b::AttitudeState)
    if a.elapsed_time != b.elapsed_time
        error("can only add States with the same time")
    end
    # linear only in angular rate
    return AttitudeState(
        a.elapsed_time, 
        a.angular_velocity_ECI_Body + b.angular_velocity_ECI_Body, 
        a.attitude_ECI_Body
    )
end

function *(c::Real, a::AttitudeState)
    # linear only in angular rate
    return AttitudeState(
        a.elapsed_time,
        c*a.angular_velocity_ECI_Body,
        a.attitude_ECI_Body
    )
end

function dcm_to_quat(X::AbstractMatrix)
    qr = 0.5 * sqrt(max(0.0, 1.0 + LinearAlgebra.tr(X)))
    qi = 1 / 4 / qr * (X[3, 2] - X[2, 3])
    qj = 1 / 4 / qr * (X[1, 3] - X[3, 1])
    qk = 1 / 4 / qr * (X[2, 1] - X[1, 2])
    return Makie.Quaternion(qi, qj, qk, qr)
end
Makie.Quaternion(X::AbstractMatrix) = begin
    qr = 0.5 * sqrt(1 + LinearAlgebra.tr(X))
    qi = 1 / 4 / qr * (X[3, 2] - X[2, 3])
    qj = 1 / 4 / qr * (X[1, 3] - X[3, 1])
    qk = 1 / 4 / qr * (X[2, 1] - X[1, 2])
    return Makie.Quaternion(qi, qj, qk, qr)
end

function run()

    server = Sockets.listen(unixsockname)
    start_time = Dates.DateTime(2026,3,20,1,2,3)
    dt = 1.0
    pwidth = 12
    counter = UInt64(0)
    
    inertia_B = diagm([5,10,13])*1e-2
    
    params = Dict{String, Any}(
        "I_Body"            => inertia_B,
        "I_Body_inv"        => inv(inertia_B),
        "max_angular_rate"  => 2e-2, # rad/s
        "battery_max"       => 84*3600.0,
        "storage_max"       => 8e9,
        "irradiance"        => 1360.0,
        "solar_panel_area"  => 0.2*0.3*2,
        "solar_panel_efficiency" => 0.29,
        "downlink_rate" => 17e6,
        "do_J2"             => true
    )
    groundstations = [
        TargetState(0x0001, 0.0, 10, pi/180*Vec3d(78.220, 15.55, 1718.0), SatelliteToolboxTransformations.r_ecef_to_eci(SatelliteToolboxTransformations.ITRF(), SatelliteToolboxTransformations.J2000(), SatelliteToolboxBase.date_to_jd(start_time), eop)*SatelliteToolboxTransformations.geodetic_to_ecef((pi/180*Vec3d(78.22, 15.55, 1718.0))...), false, false),
        TargetState(0x0002, 0.0, 9,  pi/180*Vec3d(-53.167, -70.933, 34.0), SatelliteToolboxTransformations.r_ecef_to_eci(SatelliteToolboxTransformations.ITRF(), SatelliteToolboxTransformations.J2000(), SatelliteToolboxBase.date_to_jd(start_time), eop)*SatelliteToolboxTransformations.geodetic_to_ecef((pi/180*Vec3d(-53.167, -70.933, 34.0))...), false, false),
        TargetState(0x0003, 0.0, 8,  pi/180*Vec3d(-33.86777, 151.21, 5.0), SatelliteToolboxTransformations.r_ecef_to_eci(SatelliteToolboxTransformations.ITRF(), SatelliteToolboxTransformations.J2000(), SatelliteToolboxBase.date_to_jd(start_time), eop)*SatelliteToolboxTransformations.geodetic_to_ecef((pi/180*Vec3d(-33.86777, 151.21, 5.0))...), false, false),
    ]
    
    # modelist = [
    #     ModeConfig(0x0101, "idle",      Nothing,    Makie.RGBAf(42/255, 133/255, 255/255), 5.0, 1e3, Vec3d(1.0,0.0,0.0)),
    #     ModeConfig(0x0102, "charging",  SunState,   Makie.RGBAf(255/255, 201/255, 74/255), 5.0, 10e3, Vec3d(-1.0,0.0,0.0)),
    #     ModeConfig(0x0103, "downlink",  TargetState, Makie.RGBAf(157/255, 226/255, 107/255), 23.0, 10e3, Vec3d(1.0,0.0,0.0)),
    #     ModeConfig(0x0104, "science",   Nothing,    Makie.RGBAf(206/255, 155/255, 255/255), 20.0, 10e3, Vec3d(1.0,0.0,0.0))
    # ]
    
    params["modes"] = mode_table(load_mode_config(""))
    
    r0m = 6871e3
    inc0 = 60*pi/180
    v0m = sqrt(3.986e14/r0m)*1.05
    pos_I0 = [r0m;0.0;0.0]
    vel_I0 = [0.0;v0m*cos(inc0);v0m*sin(inc0)]
    pos_state = PositionState(
        0,
        pos_I0,
        vel_I0
    )
    att_state = AttitudeState(
        0,
        rand(3)*1e-3,
        Mat3d(I)
    )
    sat_state = SatelliteState(
        IDType(0x00a0),             # ID
        0,                          # met
        Vec3d(0.0), Vec3d(0.0),     # net force, moment
        pos_I0, vel_I0,             
        rand(3)*1e-4, Mat3d(I),     # angvel, attitude
        0.0, 0.0,                   # net power, net data flows
        0.0, 0.0,                   # batter, data storage
        UInt16(0x00),               # mode
        Vec3d(NaN),                 # target dir 
        false, false                # target visible, pointed
    )
    sun_state = SunState(
        0,
        12,
        SatelliteToolboxCelestialBodies.sun_position_mod(start_time),
        false, false
    )
    earth_state = EarthState(
        0,
        SatelliteToolboxTransformations.r_eci_to_ecef(SatelliteToolboxTransformations.J2000(), SatelliteToolboxTransformations.ITRF(), SatelliteToolboxBase.date_to_jd(start_time), eop)
    )
    
    sock = Sockets.accept(server)
    println("listening...")
    
    play = Ref(UInt8(0x01))
    playrate = Ref(UInt8(0xff))
    packlen = 1
    headbuff = zeros(UInt8, 8)
    buff = zeros(UInt8, packlen)
    
    # handle received commands
    @async while isopen(sock)
        if !eof(sock)
            # read the header
            nbhead = readbytes!(sock, headbuff)
            type, flags, len = behead(headbuff)
            
            # read the rest
            buff = zeros(UInt8, len)
            nbdata = readbytes!(sock, buff)
            cmd = desmsg(type, buff)
            
            if type === PlayMessage
                play[] = cmd.message
                println("> playing: ", play[])
            end
            if type === RateMessage
                playrate[] = cmd.message
                println("> resting: ", playrate[])
            end
        else
            close(sock)
        end
    end
    
    
    time = start_time
    # run simulation
    while true
        if play[] == 0x01
            # update all dynamic systems
            push_sim!(sat_state, sun_state, earth_state, groundstations, dt, time, params)
            # update timestep
            time = start_time + Dates.Millisecond(1000*sat_state.elapsed_time)
            
            # serialize and send data
            sat_data = packetize(sat_state, 0x0000, counter)
            write(sock, sat_data)
            earth_data = packetize(earth_state, 0x0000, counter)
            write(sock, earth_data)
            sun_data = packetize(sun_state, 0x0000, counter)
            write(sock, sun_data)
            
            for target in groundstations
                gs_data = packetize(target, 0x0000, counter)
                write(sock, gs_data)
            end

            # control playback rate
            if playrate[] != 0xff
                if playrate[] == 0x00
                    sleep(0.01)
                elseif counter % playrate[] == 0
                    sleep(0.01)
                end
            end
        else
            sleep(0.1)
        end
        
        counter += 1
        
        # print(" "*worm(counter%pwidth, pwidth)*" \r")
    end
end

