module OpenSOAP
import Serialization
import Sockets
import InlineStrings
using GeometryBasics
import Dates
import PrettyTables
import SatelliteToolboxTransformations,
    SatelliteToolboxBase, SatelliteToolboxCelestialBodies
import Base: +, *, @kwdef, reinterpret, show
import Makie

global const unixsockname = "/tmp/soap_out.sock"
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
