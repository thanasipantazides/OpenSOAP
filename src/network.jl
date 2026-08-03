
# packet structure idea:
#   protocol byte ID                            [1B]
#   version ID                                  [1B]
#   message type [state, config, control, etc]  [2B]
#   flags                                       [2B]                                 
#   payload len                                 [2B]
#   payload...
# 

global const MESSAGE_TYPES = Dict{DataType,UInt16}(
    PositionState => 0x0010,
    AttitudeState => 0x0011,
    SatelliteState => 0x0020,
    MagneticFieldState => 0x0104,
    EarthState => 0x0103,
    SunState => 0x0100,
    GroundState => 0x0400,
    PlayMessage => 0x0c01,
    RateMessage => 0x0c02,
    QuitMessage => 0x0c03,
    PerturbationMessage => 0x0cd0,
)

MESSAGE_LOOKUP = Dict(value => key for (key, value) in MESSAGE_TYPES)

# a little dangerous String serializer
function reinterpret(::Type{UInt8}, s::InlineStrings.InlineString63)
    io = IOBuffer()
    write(io, s)
    v = take!(io)
    return v
end
# a String deserializer
function reinterpret(
    ::Type{InlineStrings.InlineString63},
    u::Vector{UInt8},
)::InlineStrings.InlineString63
    sarr = String(reverse(convert(Vector{Char}, u)))
    knull = findfirst('\0', sarr)
    return join(sarr[1:(knull-1)])
end

# For PositionState type
function ser(counter::UInt64, state::PositionState)::Vector{UInt8}
    c = reinterpret(UInt8, [counter])
    st = reinterpret(UInt8, [state.elapsed_time])
    sp = reinterpret(UInt8, [state.position_ECI])
    sv = reinterpret(UInt8, [state.velocity_ECI])

    res = vcat(c, st, sp, sv, UInt8[0x0a])
    return res
end

# For AttitudeState type
function ser(counter::UInt64, state::AttitudeState)::Vector{UInt8}
    c = reinterpret(UInt8, [counter])
    st = reinterpret(UInt8, [state.elapsed_time])
    sw = reinterpret(UInt8, [state.angular_velocity_ECI_Body])
    sm = reinterpret(UInt8, [state.attitude_ECI_Body])

    res = vcat(c, st, sw, sm, UInt8[0x0a])
    return res
end

function ser(counter::UInt64, state::SatelliteState)::Vector{UInt8}
    c = reinterpret(UInt8, [counter])

    f0 = reinterpret(UInt8, [state.id])
    f1 = reinterpret(UInt8, [state.elapsed_time])
    f2 = reinterpret(UInt8, [state.net_force_ECI])
    f3 = reinterpret(UInt8, [state.net_moment_Body])
    f4 = reinterpret(UInt8, [state.position_ECI])
    f5 = reinterpret(UInt8, [state.velocity_ECI])
    f6 = reinterpret(UInt8, [state.angular_velocity_ECI_Body])
    f7 = reinterpret(UInt8, [state.attitude_ECI_Body])

    f8 = reinterpret(UInt8, [state.net_power])
    f9 = reinterpret(UInt8, [state.net_data])
    f10 = reinterpret(UInt8, [state.battery_level])
    f11 = reinterpret(UInt8, [state.storage_level])
    f12 = reinterpret(UInt8, [state.mode])
    f13 = reinterpret(UInt8, [state.target])
    f14 = reinterpret(UInt8, [state.target_visible])
    f15 = reinterpret(UInt8, [state.target_pointed])

    res = vcat(
        c,
        f0,
        f1,
        f2,
        f3,
        f4,
        f5,
        f6,
        f7,
        f8,
        f9,
        f10,
        f11,
        f12,
        f13,
        f14,
        f15,
        UInt8[0x0a],
    )
    return res
end

# For MagneticFieldState type
function ser(counter::UInt64, state::MagneticFieldState)::Vector{UInt8}
    c = reinterpret(UInt8, [counter])
    id = reinterpret(UInt8, [state.id])
    st = reinterpret(UInt8, [state.elapsed_time])
    sd = reinterpret(UInt8, [state.direction_ECI])
    sv = reinterpret(UInt8, [state.visible])
    sa = reinterpret(UInt8, [state.available])
    ss = reinterpret(UInt8, [state.selected])

    res = vcat(c, id, st, sd, sv, sa, ss, UInt8[0x0a])
    return res
end

# For EarthState type
function ser(counter::UInt64, state::EarthState)::Vector{UInt8}
    c = reinterpret(UInt8, [counter])
    id = reinterpret(UInt8, [state.id])
    st = reinterpret(UInt8, [state.elapsed_time])
    sa = reinterpret(UInt8, [state.attitude_ECI_ECEF])

    res = vcat(c, id, st, sa, UInt8[0x0a])
    return res
end

# For SunState type
function ser(counter::UInt64, state::SunState)::Vector{UInt8}
    c = reinterpret(UInt8, [counter])
    id = reinterpret(UInt8, [state.id])
    st = reinterpret(UInt8, [state.elapsed_time])
    pr = reinterpret(UInt8, [state.priority])
    sa = reinterpret(UInt8, [state.position_ECI])
    vi = reinterpret(UInt8, [state.visible])
    se = reinterpret(UInt8, [state.selected])

    res = vcat(c, id, st, pr, sa, vi, se, UInt8[0x0a])
    return res
end

# For GroundState type
function ser(counter::UInt64, state::GroundState)::Vector{UInt8}
    c = reinterpret(UInt8, [counter])

    id = reinterpret(UInt8, [state.id])
    st = reinterpret(UInt8, [state.elapsed_time])
    pr = reinterpret(UInt8, [state.priority])
    lla = reinterpret(UInt8, [state.position_LLA])
    eci = reinterpret(UInt8, [state.position_ECI])
    vis = reinterpret(UInt8, [state.visible])
    sel = reinterpret(UInt8, [state.selected])

    res = vcat(c, id, st, pr, lla, eci, vis, sel, UInt8[0x0a])
    return res
end

function ser(msg::PlayMessage)::Vector{UInt8}
    return [msg.message]
end
function ser(msg::RateMessage)::Vector{UInt8}
    return [msg.message]
end
function ser(msg::QuitMessage)::Vector{UInt8}
    return [msg.message]
end
function ser(msg::PerturbationMessage)::Vector{UInt8}
    mm = reinterpret(UInt8, [msg.moment_Body])
    md = reinterpret(UInt8, [msg.moment_duration])
    fm = reinterpret(UInt8, [msg.force_ECI])
    fd = reinterpret(UInt8, [msg.force_duration])

    res = vcat(mm, md, fm, fd)
    return res
end

# general deserializer for AbstractStates
function des(::Type{T}, raw::Vector{UInt8})::Tuple{UInt64,T} where {T<:AbstractState}
    countlen = 8
    c = reinterpret(UInt64, raw[1:countlen])

    agg = countlen+1
    fields = []
    for t in T.types
        size = sizeof(t)
        d = reinterpret(t, raw[agg:(agg+size-1)])
        agg += size
        push!(fields, d[1])
    end
    return (c[1], T(fields...))
end

function des(::Type{T}, raw::Vector{UInt8})::T where {T<:ControlMessage}
    agg = 1
    fields = []
    for t in T.types
        size = sizeof(t)
        d = reinterpret(t, raw[agg:(agg+size-1)])
        agg += size
        push!(fields, d[1])
    end
    return T(fields...)
end

function packetize(payload::NetworkMessage, flags::UInt16, counter::UInt64)::Vector{UInt8}
    protocol_id = UInt8(0xab)
    protocol_ver = UInt8(0x00)
    message_type = MESSAGE_TYPES[typeof(payload)]
    if typeof(payload) <: AbstractState
        payload_B = ser(counter, payload)
    elseif typeof(payload) <: ControlMessage
        payload_B = ser(payload)
    else
        throw("Unknown payload type.")
    end
    payload_len = UInt16(length(payload_B))

    res = UInt8[
        protocol_id,
        protocol_ver,
        reverse(reinterpret(UInt8, [message_type]))...,
        reverse(reinterpret(UInt8, [flags]))...,
        reverse(reinterpret(UInt8, [payload_len]))...,
        payload_B...,
    ]
    return res
end

function behead(data::Vector{UInt8})
    if length(data) < 8
        throw(BoundsError("data too short to unpack!"))
    end
    protocol_id = data[1]
    protocol_ver = data[2]
    if protocol_id != 0xab
        throw(DomainError("protocol ID mismatch! Got " * string(protocol_id)))
    end
    if protocol_ver != 0x00
        throw(DomainError("protocol version mismatch!"))
    end

    message_type = UInt16(UInt16(data[3])<<8 | data[4])
    flags = UInt16(UInt16(data[5])<<8 | data[6])
    len = UInt16(UInt16(data[7])<<8 | data[8])

    type = MESSAGE_LOOKUP[message_type]

    return type, flags, len
end

function unpacketize(data::Vector{UInt8}, headless = true)
    if !headless
        if length(data) < 8
            throw(BoundsError("data too short to unpack!"))
        end
        protocol_id = data[1]
        protocol_ver = data[2]
        if protocol_id != 0xab
            throw(DomainError("protocol ID mismatch! Got " * string(protocol_id)))
        end
        if protocol_ver != 0x00
            throw(DomainError("protocol version mismatch!"))
        end

        message_type = UInt16(UInt16(data[3])<<8 | data[4])
        flags = UInt16(UInt16(data[5])<<8 | data[6])
        len = UInt16(UInt16(data[7])<<8 | data[8])
        type = MESSAGE_LOOKUP[message_type]
        payload = data[9:(9+len-1)]
        return des(type, payload)
    end

    # if len > length(data) - 9 - 1
    #     throw(BoundsError("data shorter than claimed!"))
    # end

    # later, return the flags + other stuff as well:

end

mutable struct SocketWrapper{
    T<:Union{Sockets.UDPSocket,Sockets.TCPSocket,Base.PipeEndpoint},
}
    sock::T
    dest::Union{Nothing,Tuple{Sockets.IPv4,UInt16}}
end

function setup_server(addr::String)
    server = Sockets.listen(addr)
    sock = Sockets.accept(server)
    return SocketWrapper{Base.PipeEndpoint}(sock, nothing)
end

function setup_server(addr::Sockets.IPv4, port::Int)
    server = Sockets.listen(addr, port)
    sock = Sockets.accept(server)
    Sockets.nagle(sock, false)
    return SocketWrapper{Sockets.TCPSocket}(sock, nothing)
end

function setup_server(
    addr::Sockets.IPv4,
    port::Int,
    dest_addr::Sockets.IPv4,
    dest_port::Int,
)

    sock = Sockets.UDPSocket()
    bind(sock, addr, port, reuseaddr = true)
    return SocketWrapper{Sockets.UDPSocket}(sock, (dest_addr, dest_port))
end

function setup_client(addr::String)
    sock = Sockets.connect(addr)
    return SocketWrapper{Base.PipeEndpoint}(sock, nothing)
end

function setup_client(addr::Sockets.IPv4, port::Int)
    sock = Sockets.connect(addr, port)
    return SocketWrapper{Sockets.TCPSocket}(sock, nothing)
end

function setup_client(
    addr::Sockets.IPv4,
    port::Int,
    dest_addr::Sockets.IPv4,
    dest_port::Int,
)

    sock = Sockets.UDPSocket()
    bind(sock, addr, port; reuseaddr = true)
    return SocketWrapper{Sockets.UDPSocket}(sock, (dest_addr, dest_port))
end

# I/O methods:
function read_transport(sock::SocketWrapper{Base.PipeEndpoint})
    headbuff = zeros(UInt8, 8)
    nbhead = readbytes!(sock.sock, headbuff)
    type, flags, len = behead(headbuff)

    buff = zeros(UInt8, len)
    nbdata = readbytes!(sock.sock, buff)
    res = des(type, buff)
    if type <: AbstractState
        count = res[1]
        data = res[2]
        return type, len, flags, count, data
    else
        return type, len, flags, res
    end
end

function write_transport(sock::SocketWrapper{Base.PipeEndpoint}, msg::Vector{UInt8})
    write(sock.sock, msg)
end

function read_transport(sock::SocketWrapper{Sockets.TCPSocket})
    headbuff = zeros(UInt8, 8)
    nbhead = readbytes!(sock.sock, headbuff)
    type, flags, len = behead(headbuff)

    buff = zeros(UInt8, len)
    nbdata = readbytes!(sock.sock, buff)
    res = des(type, buff)
    if type <: AbstractState
        count = res[1]
        data = res[2]
        return type, len, flags, count, data
    else
        return type, len, flags, res
    end
end

function write_transport(sock::SocketWrapper{Sockets.TCPSocket}, msg::Vector{UInt8})
    write(sock.sock, msg)
end

function read_transport(sock::SocketWrapper{Sockets.UDPSocket})
    msg = Sockets.recv(sock.sock)
    type, flags, len = behead(msg)
    res = des(type, msg[(8+1):end])
    if type <: AbstractState
        count = res[1]
        data = res[2]
        return type, len, flags, count, data
    else
        return type, len, flags, res
    end

    return flags, count, data
end

function write_transport(sock::SocketWrapper{Sockets.UDPSocket}, msg::Vector{UInt8})
    Sockets.send(sock.sock, sock.dest[1], sock.dest[2], msg)

end
