module OpenSOAP
import Serialization
import Sockets
import JSON
import InlineStrings
using GeometryBasics
import Dates
import DateFormats
import PrettyTables
import SatelliteToolboxTransformations,
    SatelliteToolboxBase, SatelliteToolboxCelestialBodies, SatelliteToolboxGeomagneticField
import Base: +, *, @kwdef, reinterpret, show
import Makie

global const unixsockname = "/tmp/soap_out.sock"
global const udpsockname = Sockets.@ip_str("127.0.0.1")
global const udpmonport = UInt16(9999)
global const udpcoreport = UInt16(9996)
global const IDType = UInt16
const IDDict = Dict{IDType,T} where {T}
global const eop = SatelliteToolboxTransformations.fetch_iers_eop()

include("types.jl")
export show
include("util.jl")
include("core.jl")
include("visibility.jl")
include("simulate.jl")
include("monitor.jl")
include("network.jl")
include("io.jl")

export worm,
    cross,
    uncross,
    residualSO3,
    residualso3,
    dcm_to_quat,
    r_min_arc,
    projSO3,
    NetworkMessage,
    AbstractState,
    AbstractConfig,
    AbstractTarget,
    ControlMessage,
    PlayMessage,
    RateMessage,
    QuitMessage,
    PerturbationMessage,
    PositionState,
    AttitudeState,
    EarthState,
    SatelliteState,
    SunState,
    GroundState,
    ModeConfig,
    TargetConfig,
    SimConfig,
    find_config,
    mode_table,
    ser,
    des,
    packetize,
    unpacketize,
    behead,
    setup_server,
    setup_client,
    read_transport,
    write_transport,
    step_orbit,
    step_satellite,
    step_earth,
    step!,
    run,
    run_free,
    monitor,
    write_csv,
    load_earth_texture_to_ecef,
    test_texture,
    load_jsonc,
    load_config,
    check_ids

export reinterpret

end # module OpenSOAP
