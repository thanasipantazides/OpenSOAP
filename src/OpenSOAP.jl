module OpenSOAP
import Serialization
import Sockets
import JSON
import RelocatableFolders

using StaticArrays

using LinearAlgebra
using GeometryBasics
import Dates
import DateFormats
import PrettyTables

import SatelliteToolboxTransformations
import SatelliteToolboxBase
import SatelliteToolboxCelestialBodies
import SatelliteToolboxGeomagneticField
import SatelliteToolboxAtmosphericModels

import Base:
    +,
    *,
    @kwdef,
    reinterpret,
    show,
    size,
    getindex,
    convert,
    ncodeunits,
    codeunits,
    isvalid,
    iterate,
    hash
import Makie

const SOAP_UNIX_SOCK = "/tmp/soap_out.sock"
const SOAP_HOST = Sockets.@ip_str("127.0.0.1")
const SOAP_MON_PORT = UInt16(9999)
const SOAP_CORE_PORT = UInt16(9993)
const SOAP_REPL_PORT = UInt16(9997)
const SOAP_PROTO_RX_VERSION = UInt8(0x00)
const SOAP_PROTO_TX_VERSION = UInt8(0x00)

# socket input buffer size:
const SOAP_MAX_BUFF_LEN = 2048
# max number of targets or constraints (each up to this limit) that ModeConfig can depend on.
const SOAP_MAX_ID_BUCKET = 64
const SOAP_MAX_STRING_LEN = 255
const IDType = UInt16
const IDDict = Dict{IDType,T} where {T}

const asset_path = RelocatableFolders.@path joinpath(@__DIR__, "../assets")
const eop = SatelliteToolboxTransformations.fetch_iers_eop()

SatelliteToolboxAtmosphericModels.SpaceIndices.init()

include("types.jl")
export SOAP_HOST
export SOAP_MON_PORT
export SOAP_CORE_PORT
export SOAP_MAX_STRING_LEN
export show
export convert
export SizedString
include("util.jl")
include("core.jl")
include("visibility.jl")
include("simulate.jl")
include("monitor.jl")
include("network.jl")
include("io.jl")
include("repl.jl")


export worm,
    cross,
    uncross,
    residualSO3,
    residualso3,
    dcm_to_quat,
    r_min_arc,
    projSO3,
    IDType,
    IDDict,
    IDVector,
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
    SatelliteConfig,
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
    check_ids,
    update

export reinterpret

BUILD = false
if BUILD
    include("app.jl")
    export julia_main
end

end # module OpenSOAP
