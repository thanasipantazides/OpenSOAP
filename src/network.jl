
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
    LLAConstraint => 0x0b01,
    SatelliteConfig => 0x0d00,
    GroundConfig => 0x0d01,
    SunConfig => 0x0d02,
    EarthConfig => 0x0d03,
    MagneticFieldConfig => 0x0d04,
    ModeConfig => 0x0d10,
    PlayMessage => 0x0c01,
    RateMessage => 0x0c02,
    QuitMessage => 0x0c03,
    AskMessage => 0x0c04,
    PerturbationMessage => 0x0cd0,
)

MESSAGE_LOOKUP = Dict(value => key for (key, value) in MESSAGE_TYPES)

function ser(network_object::T)::Vector{UInt8} where {T<:NetworkMessage}
    # if !isbitstype(T)
    #     throw(ArgumentError("network_object::$T is not isbits()"))
    # end

    buffer = Vector{UInt8}(undef, sizeof(T))
    ref = Ref(network_object)
    GC.@preserve ref begin
        p = Base.unsafe_convert(Ptr{T}, ref)
        unsafe_copyto!(pointer(buffer), Ptr{UInt8}(p), sizeof(T))
    end
    return buffer
end

function des(::Type{T}, buffer::AbstractVector{UInt8}) where {T<:NetworkMessage}
    ref = nothing
    if isbitstype(T)
        ref = Ref{T}()
    else
        ref = Ref{T}(T())
    end
    GC.@preserve ref begin
        p = Base.unsafe_convert(Ptr{T}, ref)
        # println("p: $p")
        unsafe_copyto!(Ptr{UInt8}(p), pointer(buffer), sizeof(T))
    end
    return ref[]
end

function packetize(payload::NetworkMessage, flags::UInt16)::Vector{UInt8}
    protocol_id = UInt8(0xab)
    protocol_ver = UInt8(SOAP_PROTO_TX_VERSION)
    message_type = MESSAGE_TYPES[typeof(payload)]

    payload_B = ser(payload)
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
    if protocol_ver != SOAP_PROTO_RX_VERSION
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
        if protocol_ver != SOAP_PROTO_RX_VERSION
            throw(DomainError("protocol version mismatch!"))
        end

        message_type = UInt16(UInt16(data[3])<<8 | data[4])
        flags = UInt16(UInt16(data[5])<<8 | data[6])
        len = UInt16(UInt16(data[7])<<8 | data[8])
        type = MESSAGE_LOOKUP[message_type]
        payload = data[9:(9+len)]
        return des(type, payload)
    else
        error("not implemented, sorry!")
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

function setup_server(addr::Sockets.IPv4, port::Integer)
    server = Sockets.listen(addr, port)
    sock = Sockets.accept(server)
    Sockets.nagle(sock, false)
    return SocketWrapper{Sockets.TCPSocket}(sock, nothing)
end

function setup_server(
    addr::Sockets.IPv4,
    port::Integer,
    dest_addr::Sockets.IPv4,
    dest_port::Integer,
)

    sock = Sockets.UDPSocket()
    bind(sock, addr, port, reuseaddr = true)
    return SocketWrapper{Sockets.UDPSocket}(sock, (dest_addr, dest_port))
end

function setup_client(addr::String)
    sock = Sockets.connect(addr)
    return SocketWrapper{Base.PipeEndpoint}(sock, nothing)
end

function setup_client(addr::Sockets.IPv4, port::Integer)
    sock = Sockets.connect(addr, port)
    return SocketWrapper{Sockets.TCPSocket}(sock, nothing)
end

function setup_client(
    addr::Sockets.IPv4,
    port::Integer,
    dest_addr::Sockets.IPv4,
    dest_port::Integer,
)

    sock = Sockets.UDPSocket()
    bind(sock, addr, port; reuseaddr = true)
    return SocketWrapper{Sockets.UDPSocket}(sock, (dest_addr, dest_port))
end

# I/O methods:
function read_transport(
    sock::SocketWrapper{Base.PipeEndpoint};
    buff = zeros(UInt8, SOAP_MAX_BUFF_LEN),
)

    nbhead = readbytes!(sock.sock, buff, 8)
    type, flags, len = behead(buff[1:8])

    nbdata = readbytes!(sock.sock, buff, len)
    data = des(type, buff[1:len])
    return type, len, flags, data
end

function write_transport(sock::SocketWrapper{Base.PipeEndpoint}, msg::Vector{UInt8})
    write(sock.sock, msg)
end

function read_transport(
    sock::SocketWrapper{Sockets.TCPSocket};
    buff = zeros(UInt8, SOAP_MAX_BUFF_LEN),
)
    nbhead = readbytes!(sock.sock, buff, 8)
    type, flags, len = behead(buff[1:8])

    nbdata = readbytes!(sock.sock, buff, len)
    data = des(type, buff[1:len])
    return type, len, flags, data
end

function write_transport(sock::SocketWrapper{Sockets.TCPSocket}, msg::Vector{UInt8})
    write(sock.sock, msg)
end

function read_transport(
    sock::SocketWrapper{Sockets.UDPSocket};
    buff = zeros(UInt8, SOAP_MAX_BUFF_LEN),
)
    msg = Sockets.recv(sock.sock)
    type, flags, len = behead(msg)
    data = des(type, msg[(8+1):end])
    return type, len, flags, data
end

function write_transport(sock::SocketWrapper{Sockets.UDPSocket}, msg::Vector{UInt8})
    Sockets.send(sock.sock, sock.dest[1], sock.dest[2], msg)

end
