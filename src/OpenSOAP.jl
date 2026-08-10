module OpenSOAP
import Serialization
import Sockets
import JSON
import RelocatableFolders
import InlineStrings
using GeometryBasics
import Dates
import DateFormats
import PrettyTables
import SatelliteToolboxTransformations,
    SatelliteToolboxBase, SatelliteToolboxCelestialBodies, SatelliteToolboxGeomagneticField
import Base: +, *, @kwdef, reinterpret, show
import Makie

const SOAP_UNIX_SOCK = "/tmp/soap_out.sock"
const SOAP_HOST = Sockets.@ip_str("127.0.0.1")
const SOAP_MON_PORT = UInt16(9999)
const SOAP_CORE_PORT = UInt16(9993)
const SOAP_MAX_BUFF_LEN = 2048
const IDType = UInt16
const IDDict = Dict{IDType,T} where {T}
const asset_path = RelocatableFolders.@path joinpath(@__DIR__, "../assets")
const eop = SatelliteToolboxTransformations.fetch_iers_eop()

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
    load_jsonc,
    load_config,
    check_ids

export reinterpret

BUILD = false
if BUILD
    include("app.jl")
    export julia_main
end

end # module OpenSOAP
